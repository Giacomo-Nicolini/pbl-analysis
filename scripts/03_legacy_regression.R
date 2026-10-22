# Spec 03: legacy regression (operational ICOS PBLH vs legacy code vs PBLHmod v2).
# Reference (amendment 2): legacy/calc_PBL_height_operational.R + legacy/proc_PBL_height.R, the
# versions sourced by the ETC pipeline (legacy/FFP_analysis_dev.R). The development version
# (calc_PBL_height_developing.R) is not used in this project (G, 2026-10-22).
# Part A: reproduce the operational FLUXES `PBLH` by running the original code on inputs rebuilt
#   from the published L2 files (R/legacy_inputs.R), per processing run (amendments 1-2), with
#   sensitivity variants and a diagnostic run on the local ETC cleanset (cfg$cleanset_file; optional).
# Part B: the port (R/legacy_port.R, default switches) equals the original on the full series.
# Part C: attribution, from the operational baseline towards PBLHmod v2: one switch at a time and a
#   cumulative path (spec 03 section 5), medians on all valid rows and on common rows (amendment 2); tables and
#   figures (section 6; png only).
# This script and tests/testthat/test-legacy_port.R are the only places allowed to source legacy/
#   (CLAUDE.md rule 5; spec 03), into a dedicated environment.
# Outputs: output/tables/03_reproduction_<site>.csv, 03_cleanset_vs_l2_<site>.csv, 03_variants_<site>.csv
#   (long: variant, stratum, statistic, value); output/figures/03_*_<site>.png/.pdf;
#   <derived_dir>/legacy/<site>_legacy_reproduction_<variant>.csv.gz, <site>_legacy_<variant>.csv.gz
#
# Usage: Rscript scripts/03_legacy_regression.R [SITE ...]  (default: cfg$test_sites)

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

# --- legacy code: packages must already be installed (never let it install them) ---
# Legacy-only dependency, declared for renv's dependency scan (never executed): the operational
# version requires only zoo (legacy/ itself is excluded from the scan by .renvignore).
if (FALSE) {
  library(zoo)
}
legacy_pkgs <- c("zoo")
missing_pkgs <- legacy_pkgs[!vapply(legacy_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cli::cli_abort(c("Legacy-only package missing: {.pkg {missing_pkgs}}.",
                   "i" = "Install it with renv (it is in renv.lock); the legacy code must not install packages."))
}
legacy_env <- new.env(parent = globalenv())
sys.source(here::here("legacy", "calc_PBL_height_operational.R"), envir = legacy_env)
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

# Processing runs of the operational column at BE-Bra (amendments 1-2): split at TIMESTAMP_START
# 2026-01-01 00:00, run 2 ends at TIMESTAMP_END 2026-08-26 00:00 (local standard time).
run_split_local <- "202601010000"
run_end_local <- "202608260000"
oper_fills <- c(1654.58, 843.07, 1591.68, 682.86) # operational fill values (STATUS; spec 04 amendment 7)

# --- original legacy code ---------------------------------------------------------------
run_original_run <- function(inp, lat, altit) {
  pass_null <- attr(inp, "pass_null")
  arg <- function(k) if (k %in% pass_null) NULL else inp[[k]]
  calc <- suppressPackageStartupMessages(legacy_env$calc_PBL_height(
    TIMESTAMP = inp$TIMESTAMP, LAT = lat, zm = inp$zm, d = inp$d, MOL = inp$MOL, ustar = inp$ustar, H = inp$H,
    AT = inp$AT, AP = arg("AP"), rho = arg("rho"), cp = arg("cp"), altit = altit,
    plot.atmc = FALSE, plot.PBL = FALSE, save.PBL = FALSE
  ))
  tmp <- tempfile("legacy_png_")
  dir.create(tmp)
  ts_proc <- as.POSIXct(inp$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  proc <- legacy_env$proc_PBL_height(TIMESTAMP = ts_proc, PBL.h = calc, plot.PBL.h.freq = FALSE, PBL.dir = paste0(tmp, "/"))
  unlink(tmp, recursive = TRUE)
  calc <- as.numeric(calc)
  proc <- as.numeric(proc)
  filled <- is.na(replace(calc, !is.na(calc) & calc < 0, NA))
  day <- as.numeric(format(ts_proc, "%H")) >= 6 & as.numeric(format(ts_proc, "%H")) <= 18
  list(proc = proc, filled = filled, fills = c(day = unique(proc[filled & day])[1], night = unique(proc[filled & !day])[1]),
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

# --- statistics ---------------------------------------------------------------------------
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

fill_stats <- function(flag_fill, oper, inp, variant) {
  hr <- as.integer(substr(inp$TIMESTAMP, 12, 13))
  day <- hr >= 6 & hr <= 18
  oper_fill <- !is.na(oper) & sprintf("%.2f", oper) %in% sprintf("%.2f", oper_fills)
  r <- rle(oper_fill)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  blk <- which(r$values)
  share <- vapply(blk, function(b) mean(flag_fill[starts[b]:ends[b]]), numeric(1))
  data.frame(variant = variant, stratum = "all",
             statistic = c("n_fill_day", "n_fill_night", "n_oper_fill_day", "n_oper_fill_night",
                           "n_oper_fill_blocks", "blocks_reproduced_50pct", "blocks_reproduced_100pct"),
             value = c(sum(flag_fill & day), sum(flag_fill & !day), sum(oper_fill & day), sum(oper_fill & !day),
                       length(blk), sum(share >= 0.5), sum(share == 1)))
}

compare_inputs <- function(cs_inp, l2_inp) {
  m <- match(as.numeric(cs_inp$timestamp_start_utc), as.numeric(l2_inp$timestamp_start_utc))
  pairs <- list(MO_LENGTH = list(cs_inp$MOL_raw, l2_inp$MOL_raw[m]), USTAR = list(cs_inp$ustar, l2_inp$ustar[m]),
                H = list(cs_inp$H, l2_inp$H[m]), AT = list(cs_inp$AT, l2_inp$AT[m]), AP = list(cs_inp$AP, l2_inp$AP[m]),
                rho = list(cs_inp$rho, l2_inp$rho[m]), cp = list(cs_inp$cp, l2_inp$cp[m]))
  do.call(rbind, lapply(names(pairs), function(k) {
    a <- pairs[[k]][[1]]
    b <- pairs[[k]][[2]]
    both <- !is.na(a) & !is.na(b)
    data.frame(variable = k, both = sum(both), na_cleanset_only = sum(is.na(a) & !is.na(b)),
               na_l2_only = sum(!is.na(a) & is.na(b)), median_abs_diff = stats::median(abs(a[both] - b[both])))
  }))
}

# Statistics of one attribution variant (spec 03 section 5.3), long format.
variant_stats <- function(v, h, flag_floor, flag_fill, ref_mod, oper, regime) {
  probs <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
  strata <- list(all = rep(TRUE, length(h)))
  for (g in c("stable", "unstable_diag", "convective")) strata[[g]] <- regime %in% g
  res <- list()
  for (s in names(strata)) {
    sel <- strata[[s]]
    q <- stats::quantile(h[sel], probs, na.rm = TRUE, names = FALSE)
    lr <- function(ref) {
      ok <- sel & !is.na(h) & !is.na(ref) & h > 0 & ref > 0
      x <- log(h[ok] / ref[ok])
      c(stats::median(x), stats::mad(x), sum(ok))
    }
    a <- lr(ref_mod)
    b <- lr(oper)
    res[[s]] <- data.frame(variant = v, stratum = paste0("regime=", s),
                           statistic = c(paste0("p", probs * 100), "share_na", "share_floored", "share_filled",
                                         "median_log_ratio_pblhmod", "mad_log_ratio_pblhmod", "n_log_pblhmod",
                                         "median_log_ratio_oper", "mad_log_ratio_oper", "n_log_oper"),
                           value = c(q, mean(is.na(h[sel])), mean(flag_floor[sel]), mean(flag_fill[sel]), a, b))
  }
  out <- do.call(rbind, res)
  out$stratum[out$stratum == "regime=all"] <- "all"
  out
}

for (site_id in site_ids) {
  t0 <- Sys.time()
  cli::cli_h1("Legacy regression: {site_id}")
  meta <- sites[sites$site_id == site_id, ]
  d <- read_icos_site(site_id, cfg, sites, sites_tv, exclusions)
  mod_path <- file.path(derived_dir, "pblhmod", paste0(site_id, "_pblhmod_default.csv.gz"))
  if (!file.exists(mod_path)) cli::cli_abort("Run {.code scripts/02_run_pblh_model.R} first ({.path {mod_path}} missing).")
  mod <- readr::read_csv(mod_path, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)

  # === Part A: reproduction ==============================================================
  cli::cli_h2("Part A: reproduction of the operational PBLH (operational legacy version)")
  variants <- list(
    rebuilt = list(split = run_split_local, ambient = "rebuilt"),
    AT_sonic = list(split = run_split_local, ambient = "AT_sonic"),
    AP_altitude = list(split = run_split_local, ambient = "AP_altitude"),
    rho_cp_legacy = list(split = run_split_local, ambient = "rho_cp_legacy"),
    one_run = list(split = NULL, ambient = "rebuilt")
  )
  cs_path <- resolve_path(cfg$cleanset_file, site = site_id, must_exist = FALSE)
  cs_inp <- NULL
  if (file.exists(cs_path)) {
    cs <- readr::read_csv(cs_path, col_select = c("TIMESTAMP_START", "TIMESTAMP_END", "MO_LENGTH", "USTAR", "H",
                                                  "AT", "AP", "RHO", "CP"),
                          col_types = readr::cols(TIMESTAMP_START = "c", TIMESTAMP_END = "c", .default = "d"),
                          na = c("-9999", "-9999.0", "NA", ""), progress = FALSE)
    cs_inp <- legacy_inputs_cleanset(cs, meta$utc_offset_h, d)
    variants$cleanset <- list(inputs = "cleanset")
  } else {
    cli::cli_alert_info("No cleanset file ({.path {cs_path}}): diagnostic run skipped.")
  }
  stats_all <- list()
  fills_all <- list()
  inputs <- list()
  results <- list()
  for (v in names(variants)) {
    vv <- variants[[v]]
    inp <- if (identical(vv$inputs, "cleanset")) cs_inp else
      legacy_inputs(d, meta, split_local = vv$split, ambient = vv$ambient, end_local = run_end_local)
    inputs[[v]] <- inp
    tm <- system.time(res <- run_original(inp, meta))
    results[[v]] <- res
    m <- match(as.numeric(inp$timestamp_start_utc), as.numeric(d$timestamp_start_utc))
    oper <- d$PBLH[m]
    year <- as.integer(substr(d$TIMESTAMP_START[m], 1, 4))
    regime <- mod$regime[match(as.numeric(inp$timestamp_start_utc), as.numeric(mod$timestamp_start_utc))]
    stats_all[[v]] <- rbind(repro_stats(res$pblh, oper, year, regime, v), fill_stats(res$flag_fill, oper, inp, v))
    fills_all[[v]] <- cbind(variant = v, res$fills)
    readr::write_csv(dplyr::tibble(timestamp_start_utc = inp$timestamp_start_utc, run_id = inp$run_id,
                                   pblh_legacy_m = res$pblh, flag_fill = res$flag_fill, pblh_oper_m = oper),
                     file.path(derived_dir, "legacy", paste0(site_id, "_legacy_reproduction_", v, ".csv.gz")))
    a <- stats_all[[v]][stats_all[[v]]$stratum == "all", ]
    cli::cli_inform("{v}: {nrow(inp)} rows in {length(unique(inp$run_id))} run{?s}, {round(tm[['elapsed']])} s; within 0.01 m {round(100 * a$value[a$statistic == 'share_within_0.01m'], 2)} %, median |rel. diff| {signif(100 * a$value[a$statistic == 'median_abs_rel_diff'], 3)} %, p95 {signif(100 * a$value[a$statistic == 'p95_abs_rel_diff'], 3)} %.")
  }
  stats_df <- do.call(rbind, stats_all)
  fills_df <- do.call(rbind, fills_all)
  fills_long <- data.frame(variant = rep(fills_df$variant, 2), stratum = paste0("run=", rep(fills_df$run_id, 2)),
                           statistic = rep(c("fill_day_m", "fill_night_m"), each = nrow(fills_df)),
                           value = c(fills_df$fill_day, fills_df$fill_night))
  readr::write_csv(rbind(stats_df, fills_long), file.path(tables_dir, paste0("03_reproduction_", site_id, ".csv")))
  cli::cli_h3("Processing runs, fill values and fill counts")
  print_tbl(fills_df)
  keep <- c("n", "share_within_0.01m", "median_abs_rel_diff", "p95_abs_rel_diff", "n_fill_day", "n_fill_night",
            "n_oper_fill_day", "n_oper_fill_night", "n_oper_fill_blocks", "blocks_reproduced_50pct",
            "blocks_reproduced_100pct")
  sa <- stats_df[stats_df$stratum == "all" & stats_df$statistic %in% keep, ]
  wa <- stats::reshape(sa, idvar = "statistic", timevar = "variant", direction = "wide", drop = "stratum")
  names(wa) <- sub("^value\\.", "", names(wa))
  print_tbl(wa[match(keep, wa$statistic), ])
  if (!is.null(cs_inp)) {
    cmp <- compare_inputs(cs_inp, inputs$rebuilt)
    readr::write_csv(cmp, file.path(tables_dir, paste0("03_cleanset_vs_l2_", site_id, ".csv")))
  }
  a <- stats_df[stats_df$variant == "rebuilt" & stats_df$stratum == "all", ]
  med <- a$value[a$statistic == "median_abs_rel_diff"]
  p95 <- a$value[a$statistic == "p95_abs_rel_diff"]
  fr <- fills_df[fills_df$variant == "rebuilt", ]
  fill_rep <- c(fr$fill_day[1], fr$fill_night[1], fr$fill_day[2], fr$fill_night[2])
  fill_dev <- abs(fill_rep / oper_fills - 1)
  accepted <- med < 0.005 && p95 < 0.02 && all(fill_dev < 0.005)
  cli::cli_inform("Acceptance (median < 0.5 %, p95 < 2 %, fills within 0.5 %): median {signif(100 * med, 3)} %, p95 {signif(100 * p95, 3)} %, max fill deviation {signif(100 * max(fill_dev), 3)} %: {if (accepted) 'ACCEPTED' else 'NOT accepted'}.")
  if (!accepted) {
    cli::cli_alert_danger("Part A not accepted: stop before parts B and C (spec 03 section 2).")
    next
  }

  # === Part B: port equality on the full series ==========================================
  cli::cli_h2("Part B: port vs original (rebuilt L2 inputs, two runs)")
  inp <- inputs$rebuilt
  tm <- system.time(port <- legacy_pblh(inp, meta))
  orig <- results$rebuilt$pblh
  same_na <- identical(is.na(port$pblh_m), is.na(orig))
  max_diff <- max(abs(port$pblh_m - orig), na.rm = TRUE)
  cli::cli_inform("Port: {round(tm[['elapsed']])} s; same NA pattern: {same_na}; max |port - original| = {signif(max_diff, 3)} m.")
  if (!same_na || max_diff >= 1e-6) cli::cli_abort("Part B failed: the port differs from the original (>= 1e-6 m).")
  cli::cli_alert_success("Part B: port equals the original (< 1e-6 m, same NA pattern).")

  # === Part C: attribution =================================================================
  cli::cli_h2("Part C: attribution, operational legacy baseline -> PBLHmod v2")
  forcing <- prepare_pblh_forcing(d, meta, pblh_options())
  mi <- match(as.numeric(inp$timestamp_start_utc), as.numeric(mod$timestamp_start_utc))
  ref_mod <- mod$pblh_m[mi]
  regime <- mod$regime[mi]
  md <- match(as.numeric(inp$timestamp_start_utc), as.numeric(d$timestamp_start_utc))
  oper <- d$PBLH[md]
  pm <- legacy_switches_pblhmod()
  base <- legacy_switches()
  set_sw <- function(sw, keys) utils::modifyList(sw, pm[keys])
  one_keys <- list(fill = "fill", floor_m = "floor_m", smoothing = "smoothing",
                   gap_interp_before_bg = "gap_interp_before_bg", mol_locf = "mol_locf", forcing = "forcing",
                   k_omega = c("k", "omega"), gamma = "gamma", neutral_formula = "neutral_formula",
                   stable_formula = "stable_formula", bg_integration = "bg_integration")
  path_keys <- list(fill = "fill", floor_m = "floor_m", smoothing = "smoothing",
                    gap_interp_before_bg = "gap_interp_before_bg", mol_locf = "mol_locf", forcing = "forcing",
                    k_omega = c("k", "omega"), gamma = "gamma", neutral_formula = "neutral_formula",
                    stable_formula = "stable_formula",
                    convective_trigger_stability_classes = c("convective_trigger", "stability_classes"),
                    bg_integration = "bg_integration")
  runs <- list(legacy_baseline = list(sw = base, inp = inp))
  for (k in names(one_keys)) runs[[paste0("legacy_", k)]] <- list(sw = set_sw(base, one_keys[[k]]), inp = inp)
  # at_meteo (one at a time): AT = METEO TA instead of the EddyPro ambient AT (sonic temperature offset, local
  # ETC QC report item 4); rows without METEO TA keep the EddyPro AT.
  at_met <- d$met_TA_K[md]
  inp_met <- inp
  inp_met$AT <- ifelse(is.na(at_met), inp$AT, at_met)
  n_at_fallback <- sum(is.na(at_met) & !is.na(inp$AT))
  runs$legacy_at_meteo <- list(sw = base, inp = inp_met)
  cum <- base
  path_names <- character(0)
  for (j in seq_along(path_keys)) {
    cum <- set_sw(cum, path_keys[[j]])
    nm <- sprintf("cum_%02d_%s", j, names(path_keys)[j])
    runs[[nm]] <- list(sw = cum, inp = inp)
    path_names <- c(path_names, nm)
  }
  one_v <- c(paste0("legacy_", names(one_keys)), "legacy_at_meteo")
  vstats <- list()
  vout <- list()
  for (v in names(runs)) {
    tm <- system.time(o <- legacy_pblh(runs[[v]]$inp, meta, runs[[v]]$sw, forcing = forcing, icos = d))
    vout[[v]] <- o
    vstats[[v]] <- variant_stats(v, o$pblh_m, o$flag_floor, o$flag_fill, ref_mod, oper, regime)
    readr::write_csv(dplyr::bind_cols(dplyr::tibble(timestamp_start_utc = inp$timestamp_start_utc, run_id = inp$run_id), o),
                     file.path(derived_dir, "legacy", paste0(site_id, "_legacy_", v, ".csv.gz")))
    cli::cli_inform("{v}: {round(tm[['elapsed']])} s.")
  }
  ref_stats <- variant_stats("pblhmod_default", ref_mod, rep(FALSE, length(ref_mod)), rep(FALSE, length(ref_mod)),
                             ref_mod, oper, regime)
  ops_stats <- variant_stats("operational", oper, !is.na(oper) & abs(oper - 50) < 1e-9,
                             !is.na(oper) & sprintf("%.2f", oper) %in% sprintf("%.2f", oper_fills), ref_mod, oper, regime)

  # Common rows (spec 03 amendment 2, review 2026-10-22): rows valid (pblh_m not NA) in every variant of the
  # analysis and in PBLHmod v2, so that medians compare values, not row sets.
  strata <- c("all", "stable", "unstable_diag", "convective")
  in_stratum <- function(g) if (g == "all") rep(TRUE, length(regime)) else regime %in% g
  common_mask <- function(vs) Reduce(`&`, c(lapply(vs, function(v) !is.na(vout[[v]]$pblh_m)), list(!is.na(ref_mod))))
  analyses <- list(one_at_a_time = c("legacy_baseline", one_v), cumulative = c("legacy_baseline", path_names))
  common_rows <- lapply(analyses, common_mask)
  common_stats <- list()
  for (an in names(analyses)) {
    cm <- common_rows[[an]]
    for (v in c(analyses[[an]], "pblhmod_default")) {
      h <- if (v == "pblhmod_default") ref_mod else vout[[v]]$pblh_m
      for (g in strata) {
        s <- cm & in_stratum(g)
        ok <- s & h > 0 & ref_mod > 0
        common_stats[[length(common_stats) + 1]] <- data.frame(
          variant = v, stratum = paste0("common_", an, "|", if (g == "all") "all" else paste0("regime=", g)),
          statistic = c("n", "p50", "median_log_ratio_pblhmod"),
          value = c(sum(s), stats::median(h[s]), stats::median(log(h[ok] / ref_mod[ok]))))
      }
    }
  }
  n_valid <- do.call(rbind, lapply(c(names(runs), "pblhmod_default"), function(v) {
    h <- if (v == "pblhmod_default") ref_mod else vout[[v]]$pblh_m
    do.call(rbind, lapply(strata, function(g) {
      data.frame(variant = v, stratum = if (g == "all") "all" else paste0("regime=", g), statistic = "n_valid",
                 value = sum(!is.na(h) & in_stratum(g)))
    }))
  }))
  vdf <- do.call(rbind, c(vstats, list(ref_stats, ops_stats, n_valid), common_stats))
  readr::write_csv(vdf, file.path(tables_dir, paste0("03_variants_", site_id, ".csv")))

  # Residual of the last cumulative step vs pblhmod_default
  last <- vout[[path_names[length(path_names)]]]$pblh_m
  res_same_na <- identical(is.na(last), is.na(ref_mod))
  res_max <- max(abs(last - ref_mod), na.rm = TRUE)
  cli::cli_h3("Residual of the fully switched legacy code vs pblhmod_default")
  cli::cli_inform("Same NA pattern: {res_same_na}; max |diff| = {signif(res_max, 3)} m.")

  # Summary tables (wide: one column per stratum)
  wide_tab <- function(df, stat, prefix = "") {
    st <- paste0(prefix, c("all", "regime=stable", "regime=unstable_diag", "regime=convective"))
    x <- df[df$statistic == stat & df$stratum %in% st, ]
    w <- stats::reshape(x[c("variant", "stratum", "value")], idvar = "variant", timevar = "stratum", direction = "wide")
    names(w) <- sub("^value\\.", "", names(w))
    names(w) <- sub(paste0("^", gsub("|", "\\|", prefix, fixed = TRUE)), "", names(w))
    names(w) <- sub("^regime=", "", names(w))
    w[c("variant", strata)]
  }
  delta <- function(tab, rows, base_row) {
    b <- tab[tab$variant == base_row, strata]
    data.frame(variant = rows, round(tab[match(rows, tab$variant), strata] - b[rep(1, length(rows)), ], 1))
  }
  steps <- function(tab, vs) {
    t2 <- tab[match(vs, tab$variant), strata]
    data.frame(step = vs[-1], round(t2[-1, ] - t2[-nrow(t2), ], 1))
  }
  nv <- wide_tab(vdf, "n_valid")
  cli::cli_h3("Valid rows (pblh_m not NA) per variant and regime")
  print_tbl(nv)
  med_all <- wide_tab(vdf, "p50")
  pv <- c("legacy_baseline", path_names)
  for (an in names(analyses)) {
    pre <- paste0("common_", an, "|")
    med_c <- wide_tab(vdf, "p50", pre)
    n_c <- wide_tab(vdf, "n", pre)
    cli::cli_h3("{an}: COMMON ROWS (n = {n_c$all[1]}; stable {n_c$stable[1]}, unstable_diag {n_c$unstable_diag[1]}, convective {n_c$convective[1]}): median PBLH [m] and change")
    print_tbl(med_c)
    if (an == "one_at_a_time") print_tbl(delta(med_c, one_v, "legacy_baseline")) else print_tbl(steps(med_c, pv))
    cli::cli_h3("{an}: ALL VALID ROWS of each variant (row sets differ): change of the median [m]")
    if (an == "one_at_a_time") print_tbl(delta(med_all, one_v, "legacy_baseline")) else print_tbl(steps(med_all, pv))
  }
  cli::cli_inform("at_meteo: METEO TA missing on {n_at_fallback} rows (EddyPro AT kept there).")

  # === Figures ===========================================================================
  cli::cli_h2("Figures")
  rr <- results$rebuilt
  figs <- list()
  figs$reproduction <- plot_legacy_reproduction(inp$timestamp_start_utc, oper, rr$pblh)
  reg_long <- function(tab, value = "median") {
    out <- do.call(rbind, lapply(strata, function(g) data.frame(step = tab$variant, regime = g, v = tab[[g]])))
    names(out)[3] <- value
    out
  }
  short <- function(x) sub("^legacy_", "", sub("^cum_[0-9]+_", "", x))
  n_lab <- function(vs) paste0(short(vs), "\n(n = ", format(nv$all[match(vs, nv$variant)], big.mark = ","), ")")
  # waterfall: common rows (main) and all valid rows (secondary file, n per step on the axis)
  med_cp <- wide_tab(vdf, "p50", "common_cumulative|")
  n_cp <- wide_tab(vdf, "n", "common_cumulative|")$all[1]
  wf <- reg_long(med_cp[match(pv, med_cp$variant), ])
  wf$step <- factor(short(wf$step), levels = short(pv))
  figs$waterfall <- plot_legacy_waterfall(wf, reg_long(med_cp[med_cp$variant == "pblhmod_default", ])[c("regime", "median")],
                                          subtitle = sprintf("Common rows (valid in every step and in PBLHmod v2): n = %s",
                                                             format(n_cp, big.mark = ",")))
  wfa <- reg_long(med_all[match(pv, med_all$variant), ])
  wfa$step <- factor(n_lab(wfa$step), levels = n_lab(pv))
  figs$waterfall_allrows <- plot_legacy_waterfall(wfa, reg_long(med_all[med_all$variant == "pblhmod_default", ])[c("regime", "median")],
                                                  subtitle = "All valid rows of each step (row sets differ; n per step)")
  # one at a time: common rows (main) and all valid rows (secondary)
  ov <- c("legacy_baseline", one_v)
  lr_c <- wide_tab(vdf, "median_log_ratio_pblhmod", "common_one_at_a_time|")
  n_co <- wide_tab(vdf, "n", "common_one_at_a_time|")$all[1]
  lo <- reg_long(lr_c[match(ov, lr_c$variant), ], "median_log_ratio")
  names(lo)[1] <- "variant"
  lo$variant <- factor(short(lo$variant), levels = rev(short(ov)))
  figs$one_at_a_time <- plot_legacy_one_at_a_time(lo, subtitle = sprintf("Common rows: n = %s", format(n_co, big.mark = ",")))
  lr_a <- wide_tab(vdf, "median_log_ratio_pblhmod")
  la <- reg_long(lr_a[match(ov, lr_a$variant), ], "median_log_ratio")
  names(la)[1] <- "variant"
  la$variant <- factor(n_lab(la$variant), levels = rev(n_lab(ov)))
  figs$one_at_a_time_allrows <- plot_legacy_one_at_a_time(la, subtitle = "All valid rows of each variant (n per variant)")
  dw <- dplyr::tibble(timestamp_start_utc = inp$timestamp_start_utc, operational = oper,
                      legacy_baseline = vout$legacy_baseline$pblh_m, legacy_no_postproc = vout$cum_03_smoothing$pblh_m,
                      pblhmod_v2 = ref_mod)
  dw$season <- season_of(dw$timestamp_start_utc)
  dw$hour_utc <- as.integer(format(dw$timestamp_start_utc, "%H", tz = "UTC"))
  figs$diurnal <- plot_legacy_diurnal(diurnal_summary(dw, c("operational", "legacy_baseline", "legacy_no_postproc", "pblhmod_v2")))
  sizes <- list(reproduction = c(24, 10), waterfall = c(18, 11), waterfall_allrows = c(18, 12), one_at_a_time = c(16, 11),
                one_at_a_time_allrows = c(16, 12), diurnal = c(17, 13))
  for (nm in names(figs)) {
    p <- save_figure(figs[[nm]], paste0("03_", nm, "_", site_id), sizes[[nm]][1], sizes[[nm]][2])
    cli::cli_alert_success("{.path {p[['png']]}}")
  }
  cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 2))}.")
}
