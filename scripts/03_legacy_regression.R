# Spec 03: legacy regression (operational ICOS PBLH vs legacy code vs PBLHmod v2).
# Part A: reproduce the operational FLUXES `PBLH` by running the original legacy
#   calc_PBL_height() + proc_PBL_height() on inputs rebuilt from the published L2 files
#   (R/legacy_inputs.R), per processing run (amendment 1), with sensitivity variants.
# This script and tests/testthat/test-legacy_port.R are the only places allowed to
# source legacy/ (CLAUDE.md rule 5; spec 03), into a dedicated environment.
# Outputs: output/tables/03_reproduction_<site>.csv (long: variant, stratum, statistic, value);
#   <derived_dir>/legacy/<site>_legacy_reproduction_<variant>.csv.gz
#
# Usage: Rscript scripts/03_legacy_regression.R [SITE ...]  (default: cfg$test_sites)

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

# --- legacy code: packages must already be installed (never let it install them) ---
# Legacy-only dependencies, declared for renv's dependency scan (never executed).
if (FALSE) {
  library(zoo)
  library(crayon)
  library(bigleaf)
  library(data.table)
}
legacy_pkgs <- c("zoo", "crayon", "bigleaf", "data.table")
missing_pkgs <- legacy_pkgs[!vapply(legacy_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cli::cli_abort(c("Legacy-only packages missing: {.pkg {missing_pkgs}}.",
                   "i" = "Install them with renv (they are in renv.lock); the legacy code must not install packages."))
}
if ("package:dplyr" %in% search()) cli::cli_abort("dplyr is attached: it would mask stats::filter() used by the legacy code.")
legacy_env <- new.env(parent = globalenv())
sys.source(here::here("legacy", "calc_PBL_height.R"), envir = legacy_env)
sys.source(here::here("legacy", "proc_PBL_height.R"), envir = legacy_env)

cfg <- get_cfg()
args <- commandArgs(trailingOnly = TRUE)
site_ids <- if (length(args) > 0) args else unlist(cfg$test_sites)
sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
derived_dir <- resolve_path(cfg$derived_dir, must_exist = FALSE)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(derived_dir, "legacy"), recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 5)

# Run the original legacy pipeline on one processing run. Returns the proc output,
# the calc output (before post-processing) and the fill values of the run.
run_original_run <- function(inp, lat, altit) {
  pass_null <- attr(inp, "pass_null")
  if (!is.na(inp$H[1]) && inp$H[1] > 0) {
    cli::cli_abort("Run starts with H > 0 ({inp$TIMESTAMP[1]}): the legacy code fails there (legacy issue 5).")
  }
  arg <- function(k) if (k %in% pass_null) NULL else inp[[k]]
  calc <- suppressPackageStartupMessages(legacy_env$calc_PBL_height(
    TIMESTAMP = inp$TIMESTAMP, LAT = lat, zm = inp$zm, d = inp$d,
    MOL = inp$MOL, ustar = inp$ustar, H = inp$H, AT = inp$AT,
    AP = arg("AP"), rho = arg("rho"), cp = arg("cp"), altit = altit,
    GF.PBL = TRUE, smooth.PBL = TRUE, plot.atmc = FALSE, plot.PBL = FALSE, save.PBL = FALSE
  ))
  tmp <- tempfile("legacy_png_")
  dir.create(tmp)
  ts_proc <- as.POSIXct(inp$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  proc <- legacy_env$proc_PBL_height(TIMESTAMP = ts_proc, PBL.h = calc, plot.PBL.h.freq = FALSE,
                                     PBL.dir = paste0(tmp, "/"))
  unlink(tmp, recursive = TRUE)
  calc <- as.numeric(calc)
  proc <- as.numeric(proc)
  filled <- is.na(replace(calc, !is.na(calc) & calc < 0, NA))
  day <- as.numeric(format(ts_proc, "%H")) >= 6 & as.numeric(format(ts_proc, "%H")) <= 18
  fills <- c(day = unique(proc[filled & day])[1], night = unique(proc[filled & !day])[1])
  list(calc = calc, proc = proc, filled = filled, fills = fills,
       n_fill = c(day = sum(filled & day), night = sum(filled & !day)))
}

run_original <- function(inp, site_meta) {
  out <- rep(NA_real_, nrow(inp))
  flag_fill <- rep(NA, nrow(inp))
  fills <- list()
  for (r in sort(unique(inp$run_id))) {
    i <- which(inp$run_id == r)
    sub <- inp[i, ]
    attr(sub, "pass_null") <- attr(inp, "pass_null")
    res <- run_original_run(sub, site_meta$lat, site_meta$altitude_m)
    out[i] <- res$proc
    flag_fill[i] <- res$filled
    fills[[length(fills) + 1]] <- data.frame(run_id = r, first = sub$TIMESTAMP[1], last = sub$TIMESTAMP[nrow(sub)],
                                             fill_day = res$fills[["day"]], n_day = res$n_fill[["day"]],
                                             fill_night = res$fills[["night"]], n_night = res$n_fill[["night"]])
  }
  list(pblh = out, flag_fill = flag_fill, fills = do.call(rbind, fills))
}

# Reproduction statistics vs the operational column (rows where it is not NA).
repro_stats <- function(leg, oper, year, regime, variant) {
  ok <- !is.na(oper) & !is.na(leg)
  dif <- leg - oper
  rel <- abs(dif) / oper
  one <- function(sel, stratum) {
    s <- ok & sel
    data.frame(variant = variant, stratum = stratum,
               statistic = c("n", "share_within_0.01m", "median_abs_rel_diff", "p95_abs_rel_diff", "max_abs_rel_diff",
                             "median_abs_diff_m"),
               value = c(sum(s), mean(abs(dif[s]) < 0.01), stats::median(rel[s]),
                         stats::quantile(rel[s], 0.95, names = FALSE), max(rel[s]), stats::median(abs(dif[s]))))
  }
  res <- list(one(rep(TRUE, length(ok)), "all"))
  for (y in sort(unique(year[ok]))) res[[length(res) + 1]] <- one(year == y, paste0("year=", y))
  for (g in sort(unique(regime[ok & !is.na(regime)]))) res[[length(res) + 1]] <- one(regime %in% g, paste0("regime=", g))
  do.call(rbind, res)
}

oper_fills <- c(1654.58, 843.07, 1591.68, 682.86) # operational fill values (STATUS; spec 04 amendment 7)

for (site_id in site_ids) {
  t0 <- Sys.time()
  cli::cli_h1("Legacy regression: {site_id}")
  meta <- sites[sites$site_id == site_id, ]
  d <- read_icos_site(site_id, cfg, sites, sites_tv, exclusions)
  mod_path <- file.path(derived_dir, "pblhmod", paste0(site_id, "_pblhmod_default.csv.gz"))
  if (!file.exists(mod_path)) cli::cli_abort("Run {.code scripts/02_run_pblh_model.R} first ({.path {mod_path}} missing).")
  mod <- readr::read_csv(mod_path, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)

  # --- Part A: reproduction --------------------------------------------------------
  cli::cli_h2("Part A: reproduction of the operational PBLH")
  variants <- list(
    rebuilt = list(split = "202601010000", ambient = "rebuilt"),
    AT_sonic = list(split = "202601010000", ambient = "AT_sonic"),
    AP_altitude = list(split = "202601010000", ambient = "AP_altitude"),
    rho_cp_legacy = list(split = "202601010000", ambient = "rho_cp_legacy"),
    one_run = list(split = NULL, ambient = "rebuilt")
  )
  stats_all <- list()
  fills_all <- list()
  for (v in names(variants)) {
    inp <- legacy_inputs(d, meta, split_local = variants[[v]]$split, ambient = variants[[v]]$ambient)
    tm <- system.time(res <- run_original(inp, meta))
    m <- match(as.numeric(inp$timestamp_start_utc), as.numeric(d$timestamp_start_utc))
    oper <- d$PBLH[m]
    year <- as.integer(substr(d$TIMESTAMP_START[m], 1, 4))
    regime <- mod$regime[match(as.numeric(inp$timestamp_start_utc), as.numeric(mod$timestamp_start_utc))]
    stats_all[[v]] <- repro_stats(res$pblh, oper, year, regime, v)
    fills_all[[v]] <- cbind(variant = v, res$fills)
    readr::write_csv(dplyr::tibble(timestamp_start_utc = inp$timestamp_start_utc, run_id = inp$run_id,
                                   pblh_legacy_m = res$pblh, flag_fill = res$flag_fill, pblh_oper_m = oper),
                     file.path(derived_dir, "legacy", paste0(site_id, "_legacy_reproduction_", v, ".csv.gz")))
    a <- stats_all[[v]][stats_all[[v]]$stratum == "all", ]
    cli::cli_inform("{v}: {nrow(inp)} rows in {length(unique(inp$run_id))} run{?s}, {round(tm[['elapsed']])} s; within 0.01 m {round(100 * a$value[a$statistic == 'share_within_0.01m'], 2)} %, median |rel. diff| {signif(100 * a$value[a$statistic == 'median_abs_rel_diff'], 3)} %, p95 {signif(100 * a$value[a$statistic == 'p95_abs_rel_diff'], 3)} %.")
  }
  stats_df <- do.call(rbind, stats_all)
  fills_df <- do.call(rbind, fills_all)
  fills_long <- data.frame(
    variant = rep(fills_df$variant, 2), stratum = paste0("run=", rep(fills_df$run_id, 2)),
    statistic = rep(c("fill_day_m", "fill_night_m"), each = nrow(fills_df)),
    value = c(fills_df$fill_day, fills_df$fill_night)
  )
  readr::write_csv(rbind(stats_df, fills_long), file.path(tables_dir, paste0("03_reproduction_", site_id, ".csv")))

  cli::cli_h3("Processing runs and fill values (operational: 1654.58 / 843.07 m, 1591.68 / 682.86 m)")
  print_tbl(fills_df)
  cli::cli_h3("Reproduction statistics (rows with operational PBLH)")
  wide <- stats::reshape(stats_df[stats_df$stratum %in% c("all") | grepl("^year=|^regime=", stats_df$stratum), ],
                         idvar = c("variant", "stratum"), timevar = "statistic", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  print_tbl(wide[wide$variant == "rebuilt", ])
  print_tbl(wide[wide$stratum == "all", ])

  # Acceptance (spec 03 section 2, amendment 1), on the rebuilt variant with the two runs.
  a <- stats_df[stats_df$variant == "rebuilt" & stats_df$stratum == "all", ]
  med <- a$value[a$statistic == "median_abs_rel_diff"]
  p95 <- a$value[a$statistic == "p95_abs_rel_diff"]
  fr <- fills_df[fills_df$variant == "rebuilt", ]
  fill_rep <- c(fr$fill_day[1], fr$fill_night[1], fr$fill_day[2], fr$fill_night[2])
  fill_dev <- abs(fill_rep / oper_fills - 1)
  cli::cli_h3("Acceptance (median |rel| < 0.5 %, p95 < 2 %, fills within 0.5 %)")
  print_tbl(data.frame(fill = c("run1 day", "run1 night", "run2 day", "run2 night"), operational = oper_fills,
                       reproduced = round(fill_rep, 2), rel_dev_pct = round(100 * fill_dev, 3)))
  accepted <- med < 0.005 && p95 < 0.02 && all(fill_dev < 0.005)
  if (accepted) {
    cli::cli_alert_success("Part A accepted: median {signif(100 * med, 3)} %, p95 {signif(100 * p95, 3)} %, fills within {signif(100 * max(fill_dev), 3)} %.")
  } else {
    cli::cli_alert_danger("Part A NOT accepted: median {signif(100 * med, 3)} %, p95 {signif(100 * p95, 3)} %, max fill deviation {signif(100 * max(fill_dev), 3)} %. Stop before part C (spec 03 section 2).")
  }
  cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 2))}.")
}
