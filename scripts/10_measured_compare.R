# Spec 10 A5: PBLHmod and ERA5 against a measured MLH with the generalised spec-07 matching
# (measured_time_base()), plus the range-compression slope (Theil-Sen of ln(source) on ln(measured), day-block
# bootstrap 95 % interval) and the median log ratio by tercile of the measured MLH.
#
# Aerosol-based reference (provider SIRTA-ReOBS in config/measured_sites.csv; FR-Gri):
#   standard sources sirta_mlh_chm15k (primary) and sirta_mlh_cl31 (script 10a), hourly means [t, t + 1 h).
#   PBLHmod: mean of the two half-hours inside the hour (both required), convective if both are convective;
#   ERA5: mean of the values at t and t + 1 h (both required), an approximation of the hourly mean.
#   Primary stratum: PBLHmod convective and apparent solar time at the interval midpoint in opts$aerosol_window
#   (10-16 h; Schween2014, Kotthaus2018); sensitivity: all convective hours. The distance class (23 km:
#   daytime_convective) keeps night and transition hours out of every metric; they appear in the diurnal figure as
#   context only.
#   Pairs: 1 PBLHmod vs CHM15k; 2 / 2b ERA5 tower cell / bilinear vs CHM15k; 3 PBLHmod vs CL31; 4 / 4b ERA5 vs CL31;
#   5 CL31 vs CHM15k (inter-method spread, "definitional uncertainty"); 6 PBLHmod vs ERA5 tower cell at the ERA5
#   hours (record, as spec 04). Figures: ERA5 tower cell only; bilinear in the tables only.
# Turbulence-based reference (provider cloudnet; FI-Hyy): only the compression table, on the exact values (qc 0/2)
#   of the spec-07 pairs 1, 2, 2b (spec-07 tables unchanged; G 2026-10-29). Exact values over-represent shallow
#   layers (deep layers are often only lower bounds, qc 3), so the slope is biased; noted in the table.
#
# Output (<output_dir>): tables/10_compression_<site>.csv; for aerosol references also tables/10_metrics_<site>.csv
# (long), tables/10_diurnal_<site>.csv, figures/10_diurnal_by_season_<site>.png, 10_scatter_<site>.png,
# 10_logratio_by_measured_tercile_<site>.png.
#
# Usage: Rscript scripts/10_measured_compare.R SITE [--n-boot=N]   (default N = 1000; small N for a quick check)
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
nb_arg <- grep("^--n-boot=", args, value = TRUE)
site_id <- setdiff(args, nb_arg)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val FR-Gri}.")
opts <- measured_compare_options()
if (length(nb_arg) == 1) opts$n_boot <- as.integer(sub("^--n-boot=", "", nb_arg))
cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
meta_site <- sites[sites$site_id == site_id, ]
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$icos_site == site_id, ]
std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)
rd <- function(src) {
  p <- file.path(std_dir, paste0(site_id, "__", src, ".csv.gz"))
  if (!file.exists(p)) cli::cli_abort("Standard file {.path {p}} not found (scripts 04, 07c, 10a).")
  read_standard_pblh(p)
}
t0 <- Sys.time()

if (any(ms$provider == "SIRTA-ReOBS")) {
  ms <- ms[ms$provider == "SIRTA-ReOBS", ]
  mode <- "aerosol"
} else if (any(ms$provider == "cloudnet")) {
  ms <- ms[ms$provider == "cloudnet", ]
  mode <- "turbulence"
} else {
  cli::cli_abort("{site_id}: no SIRTA-ReOBS or Cloudnet pairing in {.file measured_sites.csv}.")
}
if (nrow(ms) != 1) cli::cli_abort("{site_id}: more than one {ms$provider[1]} pairing.")
dclass <- distance_class(ms$distance_km)
if (dclass == "not_used") cli::cli_abort("{site_id}: measurement {ms$distance_km} km away: not used (distance rule).")
cli::cli_h1("Measured MLH comparison (spec 10): {site_id}, {ms$provider} {ms$provider_site}, {ms$distance_km} km ({dclass})")
cli::cli_inform("Bootstrap: {opts$n_boot} day-block replicates, seed {opts$seed}; aerosol window {opts$aerosol_window[1]}-{opts$aerosol_window[2]} h apparent solar time.")

mod <- rd("pblhmod_default")
era <- list(era5_blh_center = rd("era5_blh_center"), era5_blh_bilinear = rd("era5_blh_bilinear"))

# Window flag and season/hour for pair rows; `mid_s` = offset from the time stamp to the interval midpoint [s].
add_time <- function(pt, mid_s) {
  mid <- pt$timestamp_utc + mid_s
  pt$ast_mid <- apparent_solar_hour(mid, meta_site$lon)
  pt$in_window <- pt$ast_mid >= opts$aerosol_window[1] & pt$ast_mid < opts$aerosol_window[2]
  pt
}
exact_ref <- function(v) v[v$kind %in% "exact", c("timestamp_utc", "pblh_m")]
compress_rows <- function(pt, pair, set, note) {
  cs <- compression_slope(pt$mod, pt$ref, pt$timestamp_utc, opts)
  lt <- logratio_by_tercile(pt$mod, pt$ref)
  tw <- if (nrow(lt) == 3) {
    dplyr::tibble(ref_q33 = lt$ref_hi[1], ref_q67 = lt$ref_hi[2], n_T1 = lt$n[1], n_T2 = lt$n[2], n_T3 = lt$n[3],
                  mlr_T1 = lt$median_log_ratio[1], mlr_T2 = lt$median_log_ratio[2],
                  mlr_T3 = lt$median_log_ratio[3], f2_T1 = lt$frac_within_factor2[1],
                  f2_T2 = lt$frac_within_factor2[2], f2_T3 = lt$frac_within_factor2[3])
  } else {
    dplyr::tibble(ref_q33 = NA_real_)
  }
  cli::cli_inform("  {pair} [{set}]: n {cs$n}, slope {round(cs$slope, 3)} [{round(cs$slope_lo, 3)}, {round(cs$slope_hi, 3)}], reverse {round(cs$slope_reverse, 3)} ({format(round(difftime(Sys.time(), t0, units = 'mins'), 1))})")
  out <- dplyr::bind_cols(dplyr::tibble(site_id = site_id, pair = pair, stratum_set = set), cs, tw,
                          dplyr::tibble(values_note = note))
  attr(out, "terciles") <- dplyr::mutate(lt, pair = pair, stratum_set = set, .before = 1)
  out
}

if (mode == "turbulence") {
  # --- FI-Hyy: compression on the spec-07 pairs (exact values) ---
  meas <- rd("cloudnet_eps_mlh")
  t_rng <- range(meas$timestamp_start_utc)
  in_per <- function(d) d[d$timestamp_start_utc >= t_rng[1] & d$timestamp_start_utc <= t_rng[2], ]
  tb <- measured_time_base(meas, in_per(mod), lapply(era, in_per))
  prs <- list(
    `1_pblhmod_default_vs_measured` = pair_table(tb$mod, exact_ref(tb$meas), "1_pblhmod_default_vs_measured"),
    `2_era5_center_vs_measured` = pair_table(tb$era$era5_blh_center, exact_ref(tb$era_meas), "2_era5_center_vs_measured"),
    `2b_era5_bilinear_vs_measured` = pair_table(tb$era$era5_blh_bilinear, exact_ref(tb$era_meas), "2b_era5_bilinear_vs_measured")
  )
  note <- paste("exact values only (qc 0/2): shallow layers over-represented, deep layers often only lower bounds",
                "(qc 3, excluded); slope biased (spec 07)")
  cli::cli_h2("Compression (exact values)")
  comp <- list()
  for (p in names(prs)) {
    pt <- prs[[p]]
    pt <- pt[keep_regime_by_distance(pt$regime, dclass), ]
    comp[[length(comp) + 1]] <- compress_rows(pt[pt$regime %in% "convective", ], p, "convective", note)
    comp[[length(comp) + 1]] <- compress_rows(pt, p, "all_regimes", note)
  }
  comp <- dplyr::bind_rows(comp)
  comp$reference <- "cloudnet_eps_mlh"
  readr::write_csv(comp, file.path(tables_dir, paste0("10_compression_", site_id, ".csv")))
  print_tbl(dplyr::mutate(comp[c("pair", "stratum_set", "n", "n_days", "slope", "slope_lo", "slope_hi",
                                 "slope_reverse", "mlr_T1", "mlr_T2", "mlr_T3")],
                          dplyr::across(dplyr::where(is.double), ~ round(.x, 3))))
  cli::cli_alert_success("{.path {file.path(tables_dir, paste0('10_compression_', site_id, '.csv'))}}")
  cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
  quit(save = "no")
}

# --- Aerosol reference (FR-Gri, SIRTA ReOBS) ---
chm <- rd("sirta_mlh_chm15k")
cl <- rd("sirta_mlh_cl31")
v_rng <- function(d) range(d$timestamp_start_utc[!is.na(d$pblh_m)])
p_start <- max(v_rng(mod)[1], v_rng(chm)[1], v_rng(cl)[1])
p_end <- min(v_rng(mod)[2], v_rng(chm)[2], v_rng(cl)[2])
cli::cli_inform("Comparison period (valid PBLHmod, CHM15k and CL31): {format(p_start, '%Y-%m-%d %H:%M')} to {format(p_end, '%Y-%m-%d %H:%M')} UTC.")
in_per <- function(d) d[d$timestamp_start_utc >= p_start & d$timestamp_start_utc <= p_end, ]
chm <- in_per(chm)
cl <- in_per(cl)
tb_ch <- measured_time_base(chm, mod, era)
tb_cl <- measured_time_base(cl, mod, era)
cli::cli_inform(c(
  "CHM15k hours: {nrow(chm)}, valid {sum(!is.na(chm$pblh_m))} (qc 0 {sum(chm$qc_flag %in% 0L)}, qc 1 {sum(chm$qc_flag %in% 1L)}); PBLHmod valid on {sum(!is.na(tb_ch$mod$pblh_m))} of them.",
  "CL31 hours: {nrow(cl)}, valid {sum(!is.na(cl$pblh_m))} (qc 0 {sum(cl$qc_flag %in% 0L)}, qc 1 {sum(cl$qc_flag %in% 1L)})."
))

# CL31 on the CHM15k intervals, with the PBLHmod regime of the interval (pair 5).
cl_on_ch <- dplyr::tibble(timestamp_utc = tb_ch$meas$timestamp_utc,
                          pblh_m = cl$pblh_m[match(as.numeric(tb_ch$meas$timestamp_utc), as.numeric(cl$timestamp_start_utc))],
                          x_regime = tb_ch$mod$x_regime[match(as.numeric(tb_ch$meas$timestamp_utc), as.numeric(tb_ch$mod$timestamp_utc))])
era_c_in <- dplyr::rename(in_per(era$era5_blh_center), timestamp_utc = timestamp_start_utc)
mod_h <- match_halfhours_to_hours(mod, era_c_in$timestamp_utc)
pairs_def <- list(
  `1_pblhmod_default_vs_chm15k` = list(x = tb_ch$mod, ref = exact_ref(tb_ch$meas), mid = 1800, cls = "parameterisation vs aerosol"),
  `2_era5_center_vs_chm15k` = list(x = tb_ch$era$era5_blh_center, ref = exact_ref(tb_ch$meas), mid = 1800, cls = "cross-class: bulk Ri vs aerosol"),
  `2b_era5_bilinear_vs_chm15k` = list(x = tb_ch$era$era5_blh_bilinear, ref = exact_ref(tb_ch$meas), mid = 1800, cls = "cross-class: bulk Ri vs aerosol"),
  `3_pblhmod_default_vs_cl31` = list(x = tb_cl$mod, ref = exact_ref(tb_cl$meas), mid = 1800, cls = "parameterisation vs aerosol"),
  `4_era5_center_vs_cl31` = list(x = tb_cl$era$era5_blh_center, ref = exact_ref(tb_cl$meas), mid = 1800, cls = "cross-class: bulk Ri vs aerosol"),
  `4b_era5_bilinear_vs_cl31` = list(x = tb_cl$era$era5_blh_bilinear, ref = exact_ref(tb_cl$meas), mid = 1800, cls = "cross-class: bulk Ri vs aerosol"),
  `5_cl31_vs_chm15k` = list(x = cl_on_ch, ref = exact_ref(tb_ch$meas), mid = 1800, cls = "aerosol vs aerosol (CABAM vs STRATfinder): definitional uncertainty"),
  `6_pblhmod_default_vs_era5_center` = list(x = mod_h, ref = era_c_in, mid = 0, cls = "parameterisation vs model (record)")
)
pairs <- lapply(names(pairs_def), function(p) add_time(pair_table(pairs_def[[p]]$x, pairs_def[[p]]$ref, p), pairs_def[[p]]$mid))
names(pairs) <- names(pairs_def)
cls <- vapply(pairs_def, `[[`, "", "cls")

# --- Metrics: window convective (primary), all convective hours (sensitivity); pair 6 also all regimes ---
sets <- list(window_convective = function(pt) pt$regime %in% "convective" & pt$in_window,
             convective = function(pt) pt$regime %in% "convective")
metrics <- list()
for (p in names(pairs)) {
  for (s in names(sets)) {
    m <- pblh_metrics_long(pairs[[p]][sets[[s]](pairs[[p]]), ], strata = c("season", "hour_utc"))
    metrics[[length(metrics) + 1]] <- dplyr::mutate(m, stratum_set = s, .after = "pair")
  }
}
m6 <- pblh_metrics_long(pairs[["6_pblhmod_default_vs_era5_center"]])
metrics[[length(metrics) + 1]] <- dplyr::mutate(m6, stratum_set = "all_regimes", .after = "pair")
metrics <- dplyr::bind_rows(metrics)
metrics$comparison_class <- cls[metrics$pair]
metrics$distance_class <- ifelse(metrics$pair == "6_pblhmod_default_vs_era5_center", "tower (model vs model)", dclass)
metrics$hour_note <- "hour_utc = start hour of the interval [h, h + 1 h) (pairs 1-5) or ERA5 hour (pair 6)"
readr::write_csv(metrics, file.path(tables_dir, paste0("10_metrics_", site_id, ".csv")))

# --- Compression ---
cli::cli_h2("Compression slopes")
note <- "hourly means, not censored"
comp <- list()
for (p in names(pairs)) {
  for (s in names(sets)) comp[[length(comp) + 1]] <- compress_rows(pairs[[p]][sets[[s]](pairs[[p]]), ], p, s, note)
}
terc <- dplyr::bind_rows(lapply(comp, attr, "terciles"))
comp <- dplyr::bind_rows(comp)
comp$reference <- ifelse(grepl("chm15k$", comp$pair), "sirta_mlh_chm15k",
                         ifelse(grepl("cl31$", comp$pair), "sirta_mlh_cl31", "era5_blh_center"))
comp$comparison_class <- cls[comp$pair]
readr::write_csv(comp, file.path(tables_dir, paste0("10_compression_", site_id, ".csv")))

# --- Diurnal cycle (context; on the hours with a valid CHM15k value) ---
use <- tb_ch$meas$timestamp_utc[!is.na(tb_ch$meas$pblh_m)]
w <- dplyr::tibble(timestamp_utc = use)
w$sirta_mlh_chm15k <- chm$pblh_m[match(use, chm$timestamp_start_utc)]
w$sirta_mlh_cl31 <- cl$pblh_m[match(use, cl$timestamp_start_utc)]
w$pblhmod_default <- tb_ch$mod$pblh_m[match(use, tb_ch$mod$timestamp_utc)]
w$era5_blh_center <- tb_ch$era$era5_blh_center$pblh_m[match(use, tb_ch$era$era5_blh_center$timestamp_utc)]
w$season <- season_of(w$timestamp_utc)
w$hour_utc <- as.integer(format(w$timestamp_utc, "%H", tz = "UTC"))
ds <- diurnal_summary(w, c("sirta_mlh_chm15k", "sirta_mlh_cl31", "pblhmod_default", "era5_blh_center"))
ds$x_utc <- ds$hour_utc + 0.5 # interval midpoint
readr::write_csv(ds, file.path(tables_dir, paste0("10_diurnal_", site_id, ".csv")))
# Window and night per season in UTC: window = AST limits minus the median AST - UTC offset; sunrise/sunset =
# median geometric times over the days of the season in the comparison period.
days <- seq(as.Date(p_start), as.Date(p_end), by = "day")
st <- sun_times(days, meta_site$lat, meta_site$lon)
h_of <- function(t) as.numeric(t - as.POSIXct(format(t, "%Y-%m-%d"), tz = "UTC"), units = "hours")
st$season <- season_of(as.POSIXct(st$date, tz = "UTC"))
tt <- as.POSIXct(days, tz = "UTC") + 12 * 3600
off <- dplyr::tibble(season = season_of(tt), d = apparent_solar_hour(tt, meta_site$lon) - 12)
bands <- dplyr::summarise(dplyr::group_by(st, season), sunrise = stats::median(h_of(sunrise_utc)),
                          sunset = stats::median(h_of(sunset_utc)), .groups = "drop")
bands <- dplyr::left_join(bands, dplyr::summarise(dplyr::group_by(off, season), d = stats::median(d), .groups = "drop"),
                          by = "season")
bands$win_start <- opts$aerosol_window[1] - bands$d
bands$win_end <- opts$aerosol_window[2] - bands$d

# --- Figures (ERA5 tower cell only) ---
cli::cli_h2("Figures")
fig_pairs <- c(`1_pblhmod_default_vs_chm15k` = "1 PBLHmod vs CHM15k", `2_era5_center_vs_chm15k` = "2 ERA5 vs CHM15k",
               `3_pblhmod_default_vs_cl31` = "3 PBLHmod vs CL31", `4_era5_center_vs_cl31` = "4 ERA5 vs CL31",
               `5_cl31_vs_chm15k` = "5 CL31 vs CHM15k")
sc <- dplyr::bind_rows(lapply(names(fig_pairs), function(p) {
  x <- pairs[[p]][sets$window_convective(pairs[[p]]), ]
  x$panel <- fig_pairs[[p]]
  x
}))
cs <- comp[comp$pair %in% names(fig_pairs) & comp$stratum_set == "window_convective", ]
cs$panel <- unname(fig_pairs[cs$pair])
lt <- terc[terc$pair %in% names(fig_pairs), ]
lt$pair_label <- unname(fig_pairs[lt$pair])
figs <- list(
  list(plot_aerosol_diurnal(ds[c("season", "x_utc", "variable", "median")], bands), "10_diurnal_by_season", 17, 13),
  list(plot_compression_scatter(sc, cs), "10_scatter", 24, 17),
  list(plot_logratio_tercile(lt), "10_logratio_by_measured_tercile", 20, 13)
)
for (fg in figs) {
  pth <- save_figure(fg[[1]], paste0(fg[[2]], "_", site_id), fg[[3]], fg[[4]])
  cli::cli_alert_success("{.path {pth[['png']]}}")
}

# --- Report ---
get_m <- function(p, set, strata, mets = c("n", "bias", "mae", "median_log_ratio", "mad_log_ratio",
                                           "frac_within_factor2", "r_spearman")) {
  x <- metrics[metrics$pair == p & metrics$stratum_set == set & metrics$stratum %in% strata & metrics$metric %in% mets, ]
  w <- stats::reshape(as.data.frame(x[c("stratum", "metric", "value")]), idvar = "stratum", timevar = "metric",
                      direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w <- w[match(strata, w$stratum), c("stratum", intersect(mets, names(w)))]
  w[!is.na(w$stratum), ]
}
for (s in names(sets)) {
  cli::cli_h2("Pairs, {s}: all seasons")
  print_tbl(dplyr::bind_rows(lapply(names(pairs), function(p) dplyr::bind_cols(dplyr::tibble(pair = p), get_m(p, s, "all")))))
}
cli::cli_h2("Primary (window convective) by season")
for (p in names(pairs)) {
  cli::cli_h3("{p} [{cls[[p]]}]")
  print_tbl(get_m(p, "window_convective", c("all", "season=DJF", "season=MAM", "season=JJA", "season=SON")))
}
cli::cli_h3("6_pblhmod_default_vs_era5_center, all regimes (record)")
print_tbl(get_m("6_pblhmod_default_vs_era5_center", "all_regimes",
                c("all", "regime=convective", "regime=unstable_diag", "regime=mixed", "regime=stable")))
cli::cli_h2("Compression: Theil-Sen slope of ln(source) on ln(reference), day-block bootstrap 95 % interval")
print_tbl(dplyr::mutate(comp[c("pair", "stratum_set", "n", "n_days", "slope", "slope_lo", "slope_hi", "slope_reverse",
                               "ref_q33", "ref_q67", "mlr_T1", "mlr_T2", "mlr_T3")],
                        dplyr::across(dplyr::where(is.double), ~ round(.x, 3))))
cli::cli_h2("Diurnal medians [m], JJA and DJF (all hours with a CHM15k value; night = context only)")
for (s in c("JJA", "DJF")) {
  cat("--", s, "\n")
  x <- ds[ds$season == s, ]
  t <- stats::reshape(as.data.frame(x[c("hour_utc", "variable", "median")]), idvar = "hour_utc", timevar = "variable",
                      direction = "wide")
  names(t) <- sub("^median\\.", "", names(t))
  print_tbl(dplyr::mutate(t, dplyr::across(-hour_utc, round)))
}
print_tbl(dplyr::mutate(bands, dplyr::across(dplyr::where(is.double), ~ round(.x, 2))))
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
