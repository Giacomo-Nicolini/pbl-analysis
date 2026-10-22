# Spec 03: legacy regression (operational ICOS PBLH vs legacy code vs PBLHmod v2).
# Part A: reproduce the operational FLUXES `PBLH` by running the original legacy code on inputs
#   rebuilt from the published L2 files (R/legacy_inputs.R), per processing run (amendment 1).
#   Reference (amendment 2): legacy/calc_PBL_height_operational.R + legacy/proc_PBL_height.R, the
#   versions sourced by the ETC pipeline (legacy/FFP_analysis_dev.R). Also run: the development
#   version legacy/calc_PBL_height_developing.R (comparison), ambient-variable sensitivity variants,
#   one run over all years, and a diagnostic run on the local ETC cleanset (cfg$cleanset_file,
#   optional, never committed) with a comparison of cleanset and L2 inputs.
# This script and tests/testthat/test-legacy_port.R are the only places allowed to source legacy/
#   (CLAUDE.md rule 5; spec 03), into dedicated environments.
# Outputs: output/tables/03_reproduction_<site>.csv (long: variant, stratum, statistic, value);
#   output/tables/03_cleanset_vs_l2_<site>.csv; <derived_dir>/legacy/<site>_legacy_reproduction_<variant>.csv.gz
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
legacy_envs <- list(operational = new.env(parent = globalenv()), developing = new.env(parent = globalenv()))
sys.source(here::here("legacy", "calc_PBL_height_operational.R"), envir = legacy_envs$operational)
sys.source(here::here("legacy", "calc_PBL_height_developing.R"), envir = legacy_envs$developing)
for (e in legacy_envs) sys.source(here::here("legacy", "proc_PBL_height.R"), envir = e)

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

# Run the original legacy code (one version) on one processing run. Returns the proc output,
# the calc output (before post-processing) and the fill values of the run.
run_original_run <- function(inp, lat, altit, version) {
  env <- legacy_envs[[version]]
  pass_null <- attr(inp, "pass_null")
  if (version == "developing" && !is.na(inp$H[1]) && inp$H[1] > 0) {
    cli::cli_abort("Run starts with H > 0 ({inp$TIMESTAMP[1]}): the developing version fails there (legacy issue 5).")
  }
  arg <- function(k) if (k %in% pass_null) NULL else inp[[k]]
  common <- list(TIMESTAMP = inp$TIMESTAMP, LAT = lat, zm = inp$zm, d = inp$d, MOL = inp$MOL, ustar = inp$ustar,
                 H = inp$H, AT = inp$AT, AP = arg("AP"), rho = arg("rho"), cp = arg("cp"), altit = altit,
                 plot.atmc = FALSE, plot.PBL = FALSE, save.PBL = FALSE)
  if (version == "developing") common <- c(common, list(GF.PBL = TRUE, smooth.PBL = TRUE))
  calc <- suppressPackageStartupMessages(do.call(env$calc_PBL_height, common))
  tmp <- tempfile("legacy_png_")
  dir.create(tmp)
  ts_proc <- as.POSIXct(inp$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  proc <- env$proc_PBL_height(TIMESTAMP = ts_proc, PBL.h = calc, plot.PBL.h.freq = FALSE, PBL.dir = paste0(tmp, "/"))
  unlink(tmp, recursive = TRUE)
  calc <- as.numeric(calc)
  proc <- as.numeric(proc)
  filled <- is.na(replace(calc, !is.na(calc) & calc < 0, NA))
  day <- as.numeric(format(ts_proc, "%H")) >= 6 & as.numeric(format(ts_proc, "%H")) <= 18
  fills <- c(day = unique(proc[filled & day])[1], night = unique(proc[filled & !day])[1])
  list(calc = calc, proc = proc, filled = filled, fills = fills,
       n_fill = c(day = sum(filled & day), night = sum(filled & !day)))
}

run_original <- function(inp, site_meta, version) {
  out <- rep(NA_real_, nrow(inp))
  flag_fill <- rep(NA, nrow(inp))
  fills <- list()
  for (r in sort(unique(inp$run_id))) {
    i <- which(inp$run_id == r)
    sub <- inp[i, ]
    attr(sub, "pass_null") <- attr(inp, "pass_null")
    res <- run_original_run(sub, site_meta$lat, site_meta$altitude_m, version)
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

# Fill counts (day / night, from the hour of TIMESTAMP_END) and operational fill blocks reproduced:
# a block = run of consecutive operational fill rows; reproduced if >= 50 % (or 100 %) of its rows are
# filled by the legacy run.
fill_stats <- function(res, oper, inp, variant) {
  hr <- as.integer(substr(inp$TIMESTAMP, 12, 13))
  day <- hr >= 6 & hr <= 18
  oper_fill <- !is.na(oper) & sprintf("%.2f", oper) %in% sprintf("%.2f", oper_fills)
  r <- rle(oper_fill)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  blk <- which(r$values)
  share <- vapply(blk, function(b) mean(res$flag_fill[starts[b]:ends[b]]), numeric(1))
  data.frame(variant = variant, stratum = "all",
             statistic = c("n_fill_day", "n_fill_night", "n_oper_fill_day", "n_oper_fill_night",
                           "n_oper_fill_blocks", "blocks_reproduced_50pct", "blocks_reproduced_100pct",
                           "share_oper_fill_rows_filled"),
             value = c(sum(res$flag_fill & day), sum(res$flag_fill & !day), sum(oper_fill & day), sum(oper_fill & !day),
                       length(blk), sum(share >= 0.5), sum(share == 1), mean(res$flag_fill[oper_fill])))
}

# Cleanset vs L2 inputs on common rows: missing-value patterns and differences per variable.
compare_inputs <- function(cs_inp, l2_inp, cs_raw_mol, l2_raw_mol) {
  m <- match(as.numeric(cs_inp$timestamp_start_utc), as.numeric(l2_inp$timestamp_start_utc))
  pairs <- list(MO_LENGTH = list(cs_raw_mol, l2_raw_mol[m]), USTAR = list(cs_inp$ustar, l2_inp$ustar[m]),
                H = list(cs_inp$H, l2_inp$H[m]), AT = list(cs_inp$AT, l2_inp$AT[m]),
                AP = list(cs_inp$AP, l2_inp$AP[m]), rho = list(cs_inp$rho, l2_inp$rho[m]),
                cp = list(cs_inp$cp, l2_inp$cp[m]))
  do.call(rbind, lapply(names(pairs), function(k) {
    a <- pairs[[k]][[1]]
    b <- pairs[[k]][[2]]
    both <- !is.na(a) & !is.na(b)
    rel <- abs(a[both] - b[both]) / pmax(abs(a[both]), 1e-12)
    data.frame(variable = k, n = length(a), both = sum(both), na_cleanset_only = sum(is.na(a) & !is.na(b)),
               na_l2_only = sum(!is.na(a) & is.na(b)), na_both = sum(is.na(a) & is.na(b)),
               identical_share = mean(a[both] == b[both]), median_abs_diff = stats::median(abs(a[both] - b[both])),
               p99_rel_diff = stats::quantile(rel, 0.99, names = FALSE), max_rel_diff = max(rel))
  }))
}

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
    rebuilt = list(version = "operational", split = "202601010000", ambient = "rebuilt"),
    AT_sonic = list(version = "operational", split = "202601010000", ambient = "AT_sonic"),
    AP_altitude = list(version = "operational", split = "202601010000", ambient = "AP_altitude"),
    rho_cp_legacy = list(version = "operational", split = "202601010000", ambient = "rho_cp_legacy"),
    one_run = list(version = "operational", split = NULL, ambient = "rebuilt"),
    dev_rebuilt = list(version = "developing", split = "202601010000", ambient = "rebuilt")
  )
  cs_path <- resolve_path(cfg$cleanset_file, site = site_id, must_exist = FALSE)
  cs_inp <- NULL
  if (file.exists(cs_path)) {
    cs <- readr::read_csv(cs_path, col_select = c("TIMESTAMP_START", "TIMESTAMP_END", "MO_LENGTH", "USTAR", "H",
                                                  "AT", "AP", "RHO", "CP"),
                          col_types = readr::cols(TIMESTAMP_START = "c", TIMESTAMP_END = "c", .default = "d"),
                          na = c("-9999", "-9999.0", "NA", ""), progress = FALSE)
    cs_inp <- legacy_inputs_cleanset(cs, meta$utc_offset_h, d)
    cs_raw_mol <- cs$MO_LENGTH[match(as.numeric(cs_inp$timestamp_start_utc),
                                     as.numeric(icos_time_to_utc(cs$TIMESTAMP_START, meta$utc_offset_h)))]
    variants$cleanset <- list(version = "operational", inputs = "cleanset")
    cli::cli_inform("Cleanset: {nrow(cs_inp)} rows, {cs_inp$TIMESTAMP[1]} to {cs_inp$TIMESTAMP[nrow(cs_inp)]} (TIMESTAMP_END, local).")
  } else {
    cli::cli_alert_info("No cleanset file ({.path {cs_path}}): diagnostic run skipped.")
  }

  stats_all <- list()
  fills_all <- list()
  inputs <- list()
  for (v in names(variants)) {
    vv <- variants[[v]]
    inp <- if (identical(vv$inputs, "cleanset")) cs_inp else
      legacy_inputs(d, meta, split_local = vv$split, ambient = vv$ambient, tau_na = isTRUE(vv$tau_na))
    inputs[[v]] <- inp
    tm <- system.time(res <- run_original(inp, meta, vv$version))
    m <- match(as.numeric(inp$timestamp_start_utc), as.numeric(d$timestamp_start_utc))
    oper <- d$PBLH[m]
    year <- as.integer(substr(d$TIMESTAMP_START[m], 1, 4))
    regime <- mod$regime[match(as.numeric(inp$timestamp_start_utc), as.numeric(mod$timestamp_start_utc))]
    stats_all[[v]] <- rbind(repro_stats(res$pblh, oper, year, regime, v), fill_stats(res, oper, inp, v))
    fills_all[[v]] <- cbind(variant = v, res$fills)
    readr::write_csv(dplyr::tibble(timestamp_start_utc = inp$timestamp_start_utc, run_id = inp$run_id,
                                   pblh_legacy_m = res$pblh, flag_fill = res$flag_fill, pblh_oper_m = oper),
                     file.path(derived_dir, "legacy", paste0(site_id, "_legacy_reproduction_", v, ".csv.gz")))
    a <- stats_all[[v]][stats_all[[v]]$stratum == "all", ]
    cli::cli_inform("{v} ({vv$version}): {nrow(inp)} rows in {length(unique(inp$run_id))} run{?s}, {round(tm[['elapsed']])} s; within 0.01 m {round(100 * a$value[a$statistic == 'share_within_0.01m'], 2)} %, median |rel. diff| {signif(100 * a$value[a$statistic == 'median_abs_rel_diff'], 3)} %, p95 {signif(100 * a$value[a$statistic == 'p95_abs_rel_diff'], 3)} %.")
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
  cli::cli_h3("Reproduction statistics, all rows with operational PBLH")
  keep <- c("n", "share_within_0.01m", "median_abs_rel_diff", "p95_abs_rel_diff", "n_fill_day", "n_fill_night",
            "n_oper_fill_day", "n_oper_fill_night", "n_oper_fill_blocks", "blocks_reproduced_50pct",
            "blocks_reproduced_100pct", "share_oper_fill_rows_filled")
  sa <- stats_df[stats_df$stratum == "all" & stats_df$statistic %in% keep, ]
  wa <- stats::reshape(sa, idvar = c("statistic"), timevar = "variant", direction = "wide", drop = "stratum")
  names(wa) <- sub("^value\\.", "", names(wa))
  print_tbl(wa[match(keep, wa$statistic), ])
  cli::cli_h3("Operational version, rebuilt inputs, by year and regime")
  wide <- stats::reshape(stats_df[stats_df$variant == "rebuilt" & grepl("^all$|^year=|^regime=", stats_df$stratum) &
                                    stats_df$statistic %in% keep[1:4], ],
                         idvar = "stratum", timevar = "statistic", direction = "wide", drop = "variant")
  names(wide) <- sub("^value\\.", "", names(wide))
  print_tbl(wide)

  if (!is.null(cs_inp)) {
    cli::cli_h3("Cleanset vs L2 inputs (rebuilt), common rows of run 1")
    l2_raw_mol <- d$MO_LENGTH[match(as.numeric(inputs$rebuilt$timestamp_start_utc), as.numeric(d$timestamp_start_utc))]
    cmp <- compare_inputs(cs_inp, inputs$rebuilt, cs_raw_mol, l2_raw_mol)
    print_tbl(cmp)
    readr::write_csv(cmp, file.path(tables_dir, paste0("03_cleanset_vs_l2_", site_id, ".csv")))
  }

  # Acceptance (spec 03 section 2, amendments 1-2): operational version, rebuilt L2 inputs, two runs.
  a <- stats_df[stats_df$variant == "rebuilt" & stats_df$stratum == "all", ]
  med <- a$value[a$statistic == "median_abs_rel_diff"]
  p95 <- a$value[a$statistic == "p95_abs_rel_diff"]
  fr <- fills_df[fills_df$variant == "rebuilt", ]
  fill_rep <- c(fr$fill_day[1], fr$fill_night[1], fr$fill_day[2], fr$fill_night[2])
  fill_dev <- abs(fill_rep / oper_fills - 1)
  cli::cli_h3("Acceptance (median |rel| < 0.5 %, p95 < 2 %, fills within 0.5 %)")
  fill_tab <- data.frame(fill = c("run1 day", "run1 night", "run2 day", "run2 night"), operational = oper_fills,
                         rebuilt = round(fill_rep, 2), rel_dev_pct = round(100 * fill_dev, 3))
  if (!is.null(cs_inp)) {
    fc <- fills_df[fills_df$variant == "cleanset", ]
    fill_tab$cleanset <- round(c(fc$fill_day[1], fc$fill_night[1], NA, NA), 2)
  }
  print_tbl(fill_tab)
  accepted <- med < 0.005 && p95 < 0.02 && all(fill_dev < 0.005)
  if (accepted) {
    cli::cli_alert_success("Part A accepted: median {signif(100 * med, 3)} %, p95 {signif(100 * p95, 3)} %, fills within {signif(100 * max(fill_dev), 3)} %.")
  } else {
    cli::cli_alert_danger("Part A NOT accepted: median {signif(100 * med, 3)} %, p95 {signif(100 * p95, 3)} %, max fill deviation {signif(100 * max(fill_dev), 3)} %. Stop before part C (spec 03 section 2).")
  }
  cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 2))}.")
}
