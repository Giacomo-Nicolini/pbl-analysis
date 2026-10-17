# Spec 04: ERA5 readers, standard PBLH files and a first, exploratory
# comparison of PBLHmod vs ERA5 (and vs the operational ICOS PBLH).
# Inputs: ICOS L2 (spec 01), ERA5 (hourly single levels, pressure levels,
# static fields), PBLHmod outputs of scripts/02_run_pblh_model.R.
# Outputs:
#   <derived_dir>/era5/era5_cell_metadata.csv (one row per site)
#   <standard_dir>/<site>/<site>__<source>.csv.gz + .yml (5 sources)
#   <output_dir>/tables/04_metrics_<site>.csv, 04_diurnal_<site>.csv,
#     04_H_metrics_<site>.csv, 04_gamma_850_700_<site>.csv
#   <output_dir>/figures/04_*_<site>.png / .pdf (7 figures)
#
# Usage: Rscript scripts/04_era5_compare.R [SITE ...]
#   Without arguments the sites in cfg$test_sites are used.
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
args <- commandArgs(trailingOnly = TRUE)
site_ids <- if (length(args) > 0) args else unlist(cfg$test_sites)
if (length(site_ids) == 0) {
  cli::cli_abort("No sites: pass station codes as arguments or set {.field test_sites} in {.file config.yml}.")
}

sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
derived_dir <- resolve_path(cfg$derived_dir, must_exist = FALSE)
standard_dir <- resolve_path(cfg$standard_dir, must_exist = FALSE)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# PBLHmod runs compared here (same options as scripts/02_run_pblh_model.R).
runs <- list(
  default = pblh_options(),
  gamma0025 = pblh_options(gamma = 0.0025, label = "gamma0025")
)

print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)
fmt_t <- function(t) format(t, "%Y-%m-%d %H:%M", tz = "UTC")
period_row <- function(source, std) {
  v <- std[!is.na(std$pblh_m), ]
  data.frame(source = source, first_utc = fmt_t(min(v$timestamp_start_utc)),
             last_utc = fmt_t(max(v$timestamp_end_utc)), n_rows = nrow(std), n_valid = nrow(v))
}
as_hourly <- function(std) { # instantaneous standard table -> timestamp_utc
  dplyr::rename(std, timestamp_utc = timestamp_start_utc)
}

for (site_id in site_ids) {
  t0 <- Sys.time()
  cli::cli_h1("ERA5 comparison: {site_id}")
  meta <- sites[sites$site_id == site_id, ]
  if (nrow(meta) != 1) cli::cli_abort("Site {.val {site_id}} not found in {.file {cfg$sites_file}}.")

  # --- 1. Inputs ---------------------------------------------------------------
  cli::cli_h2("Inputs")
  d <- read_icos_site(site_id, cfg, sites, sites_tv, exclusions)
  sl <- find_era5_files(site_id, cfg, "single_level", sites = sites)
  blh <- read_era5_single_level(sl, "boundary_layer_height", cells = "all")
  cbh <- read_era5_single_level(sl, "cloud_base_height", cells = "center")
  hsf <- read_era5_single_level(sl, "surface_sensible_heat_flux", cells = "center")
  pl <- read_era5_pressure_levels(find_era5_files(site_id, cfg, "pressure_level", sites = sites), cells = "center")
  static <- read_era5_static(site_id, cfg, sites)

  # --- 2. Cell metadata ----------------------------------------------------------
  cm <- era5_cell_metadata(site_id, cfg, sites, static)
  cm_dir <- file.path(derived_dir, "era5")
  dir.create(cm_dir, recursive = TRUE, showWarnings = FALSE)
  cm_path <- file.path(cm_dir, "era5_cell_metadata.csv")
  if (file.exists(cm_path)) {
    old <- readr::read_csv(cm_path, show_col_types = FALSE, progress = FALSE)
    cm_all <- dplyr::bind_rows(old[old$site_id != site_id, ], cm)
  } else {
    cm_all <- cm
  }
  readr::write_csv(cm_all[order(cm_all$site_id), ], cm_path)

  # --- 3. PBLHmod outputs ---------------------------------------------------------
  mods <- list()
  for (lab in names(runs)) {
    p <- file.path(derived_dir, "pblhmod", paste0(site_id, "_pblhmod_", lab, ".csv.gz"))
    if (!file.exists(p)) {
      cli::cli_abort(c("PBLHmod output not found: {.path {p}}.", "i" = "Run {.code Rscript scripts/02_run_pblh_model.R {site_id}} first."))
    }
    mods[[lab]] <- readr::read_csv(p, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
  }

  # --- 4. Standard files ----------------------------------------------------------
  cli::cli_h2("Standard files")
  center <- blh[blh$cell_row == 0 & blh$cell_col == 0, ]
  bilin <- era5_bilinear(blh, meta$lat, meta$lon, "blh_m")
  spread <- era5_spread(blh, "blh_m")
  hdr <- sl$header[[which(sl$variable == "boundary_layer_height" & sl$cell_row == 0 & sl$cell_col == 0)]]
  std <- list(
    pblhmod_default = standardize_pblhmod(mods$default, "default", runs$default),
    pblhmod_gamma0025 = standardize_pblhmod(mods$gamma0025, "gamma0025", runs$gamma0025),
    era5_blh_center = standardize_era5_blh(center, bilin, spread, "center", exclusions, meta$utc_offset_h, hdr),
    era5_blh_bilinear = standardize_era5_blh(center, bilin, spread, "bilinear", exclusions, meta$utc_offset_h, hdr),
    icos_oper = standardize_icos_oper(d)
  )
  site_std_dir <- file.path(standard_dir, site_id)
  for (s in names(std)) {
    paths <- write_standard_pblh(std[[s]]$data, std[[s]]$meta, site_std_dir)
    cli::cli_alert_success("{s}: {nrow(std[[s]]$data)} rows -> {.path {paths[['data']]}}")
  }

  # --- 5. Matching and comparisons ------------------------------------------------
  era_c <- as_hourly(std$era5_blh_center$data)
  era_b <- as_hourly(std$era5_blh_bilinear$data)
  hours <- era_c$timestamp_utc
  h_def <- match_halfhours_to_hours(std$pblhmod_default$data, hours)
  h_g <- match_halfhours_to_hours(std$pblhmod_gamma0025$data, hours)
  h_op <- match_halfhours_to_hours(std$icos_oper$data, hours)
  h_op_qc0 <- h_op[!is.na(h_op$qc_flag) & h_op$qc_flag == 0, ]
  op_qc0_hh <- std$icos_oper$data[!is.na(std$icos_oper$data$qc_flag) & std$icos_oper$data$qc_flag == 0, ]
  def_hh <- dplyr::rename(std$pblhmod_default$data, timestamp_utc = timestamp_start_utc)

  pairs <- dplyr::bind_rows(
    pair_table(h_def, era_c, "1_pblhmod_default_vs_era5_center"),
    pair_table(h_def, era_b, "2_pblhmod_default_vs_era5_bilinear"),
    pair_table(h_g, era_c, "3_pblhmod_gamma0025_vs_era5_center"),
    pair_table(h_op_qc0, era_c, "4a_icos_oper_qc0_vs_era5_center"),
    pair_table(h_op, era_c, "4b_icos_oper_all_vs_era5_center"),
    # pair 5: both half-hourly, compared on the half-hour (start time)
    pair_table(def_hh, dplyr::rename(op_qc0_hh, timestamp_utc = timestamp_start_utc),
               "5_pblhmod_default_vs_icos_oper_qc0_halfhourly")
  )
  # PBLHmod regime as stratum for pairs 4a/4b (operational has no regime)
  reg_h <- stats::setNames(h_def$x_regime, as.character(as.numeric(h_def$timestamp_utc)))
  is4 <- startsWith(pairs$pair, "4")
  pairs$regime[is4] <- reg_h[as.character(as.numeric(pairs$timestamp_utc[is4]))]
  metrics <- pblh_metrics_long(pairs)
  readr::write_csv(metrics, file.path(tables_dir, paste0("04_metrics_", site_id, ".csv")))

  # Diurnal cycles over the overlap of PBLHmod (default), ERA5 (centre) and operational PBLH
  per <- do.call(rbind, lapply(names(std), function(s) period_row(s, std[[s]]$data)))
  ov_start <- max(as.POSIXct(per$first_utc[per$source %in% c("pblhmod_default", "era5_blh_center", "icos_oper")], tz = "UTC"))
  ov_end <- min(as.POSIXct(per$last_utc[per$source %in% c("pblhmod_default", "era5_blh_center", "icos_oper")], tz = "UTC"))
  wide <- dplyr::tibble(timestamp_utc = hours, era5_blh_center = era_c$pblh_m[match(hours, era_c$timestamp_utc)])
  wide$pblhmod_default <- h_def$pblh_m[match(hours, h_def$timestamp_utc)]
  wide$icos_oper_all <- h_op$pblh_m[match(hours, h_op$timestamp_utc)]
  wide$icos_oper_qc0 <- ifelse(h_op$qc_flag[match(hours, h_op$timestamp_utc)] %in% 0L, wide$icos_oper_all, NA)
  wide <- wide[wide$timestamp_utc >= ov_start & wide$timestamp_utc <= ov_end, ]
  wide$season <- season_of(wide$timestamp_utc)
  wide$hour_utc <- as.integer(format(wide$timestamp_utc, "%H", tz = "UTC"))
  ds <- diurnal_summary(wide, c("pblhmod_default", "era5_blh_center", "icos_oper_qc0", "icos_oper_all"))
  readr::write_csv(ds, file.path(tables_dir, paste0("04_diurnal_", site_id, ".csv")))

  # --- 6. Extras -------------------------------------------------------------------
  # EC H (as used by PBLHmod) vs ERA5 H of the tower cell
  ech <- ec_h_from_pblhmod(mods$default[!mods$default$excl_flag, ])
  hm <- match_halfhours_to_accum(ech, hsf$timestamp_end_utc, "H_ec_Wm2")
  hd <- dplyr::inner_join(dplyr::rename(hm, H_ec_Wm2 = value),
                          dplyr::select(hsf, timestamp_utc = timestamp_end_utc, H_era5_Wm2), by = "timestamp_utc")
  hd$season <- season_of(hd$timestamp_utc - 1800)
  hd_ok <- hd[!is.na(hd$H_ec_Wm2) & !is.na(hd$H_era5_Wm2), ]
  abs_metrics <- c("n", "bias", "mae", "rmse", "r_pearson", "r_spearman")
  hmet <- dplyr::bind_rows(
    dplyr::bind_cols(dplyr::tibble(stratum = "all"), pblh_metrics(hd_ok$H_ec_Wm2, hd_ok$H_era5_Wm2)),
    dplyr::mutate(pblh_metrics_by(hd_ok, "H_ec_Wm2", "H_era5_Wm2", "season"),
                  stratum = paste0("season=", season), .keep = "unused")
  )
  hmet <- hmet[c("stratum", abs_metrics)]
  hmet_long <- dplyr::tibble(
    pair = "EC_H_vs_ERA5_H", stratum = rep(hmet$stratum, times = length(abs_metrics)),
    metric = rep(abs_metrics, each = nrow(hmet)),
    value = unlist(lapply(abs_metrics, function(k) as.numeric(hmet[[k]])), use.names = FALSE)
  )
  readr::write_csv(hmet_long, file.path(tables_dir, paste0("04_H_metrics_", site_id, ".csv")))

  # Cloud-base classes
  p1 <- pairs[pairs$pair == "1_pblhmod_default_vs_era5_center", ]
  p1$cloud <- cloud_class(cbh$cbh_m[match(as.numeric(p1$timestamp_utc), as.numeric(cbh$timestamp_utc))])

  # Free-atmosphere theta gradient
  g <- gamma_850_700(pl)
  gm <- gamma_monthly(g)
  readr::write_csv(gm, file.path(tables_dir, paste0("04_gamma_850_700_", site_id, ".csv")))

  # --- 7. Figures ------------------------------------------------------------------
  cli::cli_h2("Figures")
  figs <- list(
    list(plot_diurnal_by_season(ds), "04_diurnal_by_season", 17, 13),
    list(plot_scatter_pair(p1), "04_scatter_mod_era5", 24, 8),
    list(plot_logratio_by_hour(p1), "04_logratio_by_hour", 17, 13),
    list(plot_era5_spread(std$era5_blh_center$data, std$era5_blh_bilinear$data), "04_era5_spread", 20, 9),
    list(plot_h_ec_vs_era5(hd), "04_H_ec_vs_era5", 22, 10),
    list(plot_logratio_by_cloud(p1), "04_logratio_by_cloud", 17, 8),
    list(plot_gamma_monthly(gm), "04_gamma_850_700", 14, 9)
  )
  for (fg in figs) {
    pth <- save_figure(fg[[1]], paste0(fg[[2]], "_", site_id), fg[[3]], fg[[4]])
    cli::cli_alert_success("{.path {pth[['png']]}}")
  }

  # --- 8. Summary ------------------------------------------------------------------
  cli::cli_h2("Summary: {site_id}")
  cli::cli_h3("Periods (valid PBLH) and overlap")
  print_tbl(per)
  cli::cli_inform("Overlap of PBLHmod default, ERA5 centre and operational PBLH: {fmt_t(ov_start)} to {fmt_t(ov_end)} UTC.")

  cli::cli_h3("Rows dropped by exclusions")
  print_tbl(data.frame(source = names(std), n_excluded = vapply(std, function(s) as.integer(s$n_excluded), 1L)))

  cli::cli_h3("Matched hours per pair")
  era_ok <- !is.na(era_c$pblh_m)
  lost <- function(h) sum(h$n_samples[match(era_c$timestamp_utc[era_ok], h$timestamp_utc)] == 1, na.rm = TRUE)
  mt <- dplyr::summarise(pairs, .by = pair, n_matched = dplyr::n())
  mt$lost_require_both <- c(lost(h_def), lost(h_def), lost(h_g), NA, lost(h_op), NA)[match(mt$pair, unique(pairs$pair))]
  print_tbl(mt)
  cli::cli_inform("lost_require_both: ERA5 hours with only one valid half-hour (NA for pairs 4a and 5: qc subset / half-hourly).")

  cli::cli_h3("Operational PBLH: floor and fill values")
  cli::cli_inform("{std$icos_oper$n_floor} floored rows (50 m).")
  cli::cli_bullets(stats::setNames(format_fill_values(std$icos_oper$fill_values), rep("*", nrow(std$icos_oper$fill_values))))

  cli::cli_h3("Pair 1 (PBLHmod default vs ERA5 centre): metrics by regime and season")
  keep <- c("n", "bias", "mae", "median_log_ratio", "mad_log_ratio", "frac_within_factor2", "r_spearman")
  m1 <- metrics[metrics$pair == "1_pblhmod_default_vs_era5_center" & metrics$metric %in% keep &
                  (metrics$stratum == "all" | startsWith(metrics$stratum, "regime=") | startsWith(metrics$stratum, "season=")), ]
  m1w <- stats::reshape(as.data.frame(m1[c("stratum", "metric", "value")]), idvar = "stratum", timevar = "metric",
                        direction = "wide")
  names(m1w) <- sub("^value\\.", "", names(m1w))
  print_tbl(m1w[c("stratum", keep)])

  cli::cli_h3("All pairs: overall metrics")
  ma <- metrics[metrics$stratum == "all" & metrics$metric %in% keep, ]
  maw <- stats::reshape(as.data.frame(ma[c("pair", "metric", "value")]), idvar = "pair", timevar = "metric", direction = "wide")
  names(maw) <- sub("^value\\.", "", names(maw))
  print_tbl(maw[c("pair", keep)])

  cli::cli_h3("EC H vs ERA5 H (descriptive)")
  print_tbl(hmet)

  cli::cli_h3("log(PBLHmod / ERA5) by cloud class and regime (descriptive)")
  pc <- p1[p1$mod > 0 & p1$ref > 0 & !is.na(p1$regime), ]
  cs <- dplyr::summarise(pc, .by = c(regime, cloud), n = dplyr::n(), median_log_ratio = stats::median(log(mod / ref)))
  print_tbl(cs[order(cs$regime, cs$cloud), ])

  cli::cli_h3("gamma_850_700 monthly median [K m-1]")
  print_tbl(gm)
  cli::cli_inform("Monthly median range: 00 UTC {signif(min(gm$median[gm$hour_utc == 0]), 3)} to {signif(max(gm$median[gm$hour_utc == 0]), 3)}; 12 UTC {signif(min(gm$median[gm$hour_utc == 12]), 3)} to {signif(max(gm$median[gm$hour_utc == 12]), 3)} K m-1 (PBLHmod: 0.01 default, 0.0025 legacy).")

  cli::cli_h3("ERA5 cell metadata")
  print_tbl(cm)

  cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 2))}.")
}
