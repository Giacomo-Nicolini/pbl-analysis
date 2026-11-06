# Spec 10 A5: PBLHmod and ERA5 against a measured MLH with the generalised spec-07 matching
# (measured_time_base()), plus the range-compression statistics (G 2026-10-29): symmetric (RMA) slope of
# ln(source) on ln(measured) and Bland-Altman trend of ln(source / measured) on the mean of the logs, both with a
# day-block bootstrap 95 % interval (main evidence); forward and reverse Theil-Sen slopes and the median log ratio
# by tercile of the measured MLH as diagnostics (regression dilution / regression to the mean). Range compression
# is a hypothesis under test.
#
# Aerosol-based reference (primary source with method_class aerosol in config/measured_sites.csv):
#   FR-Gri (SIRTA-ReOBS): sirta_mlh_chm15k (primary) and sirta_mlh_cl31 (script 10a), hourly means [t, t + 1 h);
#   PBLHmod: mean of the two half-hours inside the hour (both required), convective if both are convective;
#   ERA5: mean of the values at t and t + 1 h (both required), an approximation of the hourly mean.
#   MeteoSwiss sites (spec 13; script 11b): one half-hourly source, meteoswiss_stratfinder_mlh or meteoswiss_cabam_mlh;
#   PBLHmod on the same half-hours; ERA5 at its hours against the mean of the two adjacent half-hours
#   (measured_time_base()). Pairs 1, 2, 2b and 6 only (no second instrument).
#   Primary stratum: PBLHmod convective and apparent solar time at the interval midpoint in opts$aerosol_window
#   (10-16 h; Schween2014, Kotthaus2018); sensitivity: all convective hours. The distance class (23 km:
#   daytime_convective) keeps night and transition hours out of every metric; in the figures the hours outside the
#   primary set are drawn faded as context (spec 11 E, plot_utils.R).
#   Pairs: 1 PBLHmod vs CHM15k; 2 / 2b ERA5 tower cell / bilinear vs CHM15k; 3 PBLHmod vs CL31; 4 / 4b ERA5 vs CL31;
#   5 CL31 vs CHM15k (inter-method spread, "definitional uncertainty"); 6 PBLHmod vs ERA5 tower cell at the ERA5
#   hours (record, as spec 04). Figures: ERA5 tower cell only; bilinear in the tables only.
# Turbulence-based reference (provider cloudnet; FI-Hyy): compression statistics on the exact values (qc 0/2) of
#   the spec-07 pairs 1, 2, 2b (spec-07 tables unchanged; G 2026-10-29); exact values over-represent shallow layers
#   (deep layers are often only lower bounds, qc 3), so the statistics are biased (noted in the table). Plus the
#   scatter-only null for the censored tails (censored_tail_null()).
#
# Metrics (pblh_metrics()) include the fractional bias and normalised RMSD of Chang2004 as in Sinclair2022
# (Eqs. 6-7), with the measured (or, in pair 6, the ERA5) value as the observation: fb > 0 when the observation is
# deeper.
#
# Output (<output_dir>): tables/10_compression_<site>.csv, 10_bland_altman_bins_<site>.csv,
# figures/10_bland_altman_<site>.png; for aerosol references also tables/10_metrics_<site>.csv (long),
# 10_diurnal_<site>.csv, figures/10_diurnal_by_season_<site>.png, 10_scatter_<site>.png; for turbulence references
# also tables/10_censored_null_<site>.csv.
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
fb_note <- paste("fb = (mean_obs - mean_source) / (0.5 (mean_obs + mean_source)), obs = reference column (pair 6:",
                 "ERA5): fb > 0 when the observation is deeper (Sinclair2022 Eq. 6, Chang2004);",
                 "nrmsd = sqrt(mean((obs - source)^2) / (mean_obs mean_source)) (Sinclair2022 Eq. 7)")
tercile_note <- "terciles of the measured MLH: affected by regression to the mean (diagnostic only)"

# The site's primary measured source decides the mode (spec 13): aerosol (SIRTA-ReOBS, MeteoSwiss) or turbulence
# (Cloudnet).
ms <- ms[!is.na(ms$primary_source_id), ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: {nrow(ms)} primary measured source{?s} in {.file measured_sites.csv} (need 1).")
mode <- ms$method_class
if (mode == "turbulence" && ms$provider != "cloudnet") cli::cli_abort("{site_id}: turbulence references only from Cloudnet.")
if (!mode %in% c("aerosol", "turbulence")) cli::cli_abort("{site_id}: method_class {.val {mode}} not supported.")
dclass <- distance_class(ms$distance_km)
if (dclass == "not_used") cli::cli_abort("{site_id}: measurement {ms$distance_km} km away: not used (distance rule).")
# Tables carry the pairing distance and the borderline flag (spec 13, G 2026-11-06)
write_tab <- function(x, path) readr::write_csv(tag_measured_pairing(x, site_id, ms), path)
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
# Compression statistics of one pair and stratum set; Bland-Altman band and bins kept as attributes.
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
  cli::cli_inform("  {pair} [{set}]: n {cs$n}, symmetric {round(cs$slope_sym, 3)} [{round(cs$slope_sym_lo, 3)}, {round(cs$slope_sym_hi, 3)}], BA trend {round(cs$ba_slope, 3)} [{round(cs$ba_slope_lo, 3)}, {round(cs$ba_slope_hi, 3)}]")
  out <- dplyr::bind_cols(dplyr::tibble(site_id = site_id, pair = pair, stratum_set = set), cs, tw,
                          dplyr::tibble(values_note = note, tercile_note = tercile_note))
  attr(out, "ba_band") <- dplyr::mutate(attr(cs, "ba_band"), pair = pair, stratum_set = set, .before = 1)
  attr(out, "ba_bins") <- dplyr::mutate(bland_altman_bins(pt$mod, pt$ref, opts$n_bins), pair = pair,
                                        stratum_set = set, .before = 1)
  out
}
bind_comp <- function(comp) {
  list(comp = dplyr::bind_rows(comp), band = dplyr::bind_rows(lapply(comp, attr, "ba_band")),
       bins = dplyr::bind_rows(lapply(comp, attr, "ba_bins")))
}
comp_cols <- c("pair", "stratum_set", "n", "n_days", "r_log", "slope_sym", "slope_sym_lo", "slope_sym_hi", "ba_mean",
               "ba_slope", "ba_slope_lo", "ba_slope_hi", "ts_slope", "ts_slope_reverse")
# Bland-Altman figure of the pairs in `labels` (stratum set `set`).
ba_figure <- function(pairs, cb, labels, set, subtitle) {
  lab <- function(d) {
    if ("stratum_set" %in% names(d)) d <- d[d$stratum_set == set, ]
    d <- d[d$pair %in% names(labels), ]
    d$panel <- factor(unname(labels[d$pair]), levels = unname(labels))
    d
  }
  pt <- dplyr::bind_rows(lapply(names(labels), function(p) dplyr::mutate(pairs[[p]], pair = p)))
  plot_bland_altman(lab(pt), lab(cb$bins), lab(cb$band), lab(cb$comp), subtitle)
}

if (mode == "turbulence") {
  # --- FI-Hyy: compression on the spec-07 pairs (exact values) and the censored-tail null ---
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
                "(qc 3, excluded); statistics biased (spec 07)")
  cli::cli_h2("Compression (exact values)")
  comp <- list()
  for (p in names(prs)) {
    pt <- prs[[p]]
    pt <- pt[keep_regime_by_distance(pt$regime, dclass), ]
    prs[[p]] <- pt
    comp[[length(comp) + 1]] <- compress_rows(pt[pt$regime %in% "convective", ], p, "convective", note)
    comp[[length(comp) + 1]] <- compress_rows(pt, p, "all_regimes", note)
  }
  cb <- bind_comp(comp)
  cb$comp$reference <- "cloudnet_eps_mlh"
  write_tab(cb$comp, file.path(tables_dir, paste0("10_compression_", site_id, ".csv")))
  write_tab(cb$bins, file.path(tables_dir, paste0("10_bland_altman_bins_", site_id, ".csv")))

  # --- Scatter-only null for the censored tails (all values: exact, below z_min, lower bound) ---
  cli::cli_h2("Censored tails: observed shares vs scatter-only null")
  null_df <- function(m, v) {
    i <- match(as.numeric(m$timestamp_utc), as.numeric(v$timestamp_utc))
    d <- dplyr::tibble(timestamp_utc = m$timestamp_utc, mod = m$pblh_m, kind = v$kind[i], ref = v$pblh_m[i],
                       z_min = v$z_min[i], regime = m$x_regime, season = season_of(m$timestamp_utc))
    d[!is.na(d$mod) & !is.na(d$kind) & keep_regime_by_distance(d$regime, dclass), ]
  }
  nd <- list(`1_pblhmod_default_vs_measured` = null_df(tb$mod, tb$meas),
             `2_era5_center_vs_measured` = null_df(tb$era$era5_blh_center, tb$era_meas),
             `2b_era5_bilinear_vs_measured` = null_df(tb$era$era5_blh_bilinear, tb$era_meas))
  null_sets <- list(all_regimes = function(d) rep(TRUE, nrow(d)), convective = function(d) d$regime %in% "convective",
                    stable = function(d) d$regime %in% "stable")
  nul <- dplyr::bind_rows(lapply(names(nd), function(p) {
    dplyr::bind_rows(lapply(names(null_sets), function(s) {
      d <- nd[[p]][null_sets[[s]](nd[[p]]), ]
      dplyr::bind_cols(dplyr::tibble(site_id = site_id, pair = p, regime_set = s), censored_tail_null(d, opts))
    }))
  }))
  nul$tail_meaning <- ifelse(nul$tail == "below", "share of source <= z_min where the measurement is below z_min",
                             "share of source < bound where the measurement is a lower bound (qc 3)")
  nul$sim_note <- ifelse(nul$sim_bound == "lower",
                         "null share is a lower bound (true depth <= z_min): observed below it = source too deep",
                         "null share is an upper bound (true depth >= bound): observed above it = source too shallow")
  write_tab(nul, file.path(tables_dir, paste0("10_censored_null_", site_id, ".csv")))
  print_tbl(dplyr::mutate(nul[c("pair", "regime_set", "tail", "n", "n_days", "obs_share", "obs_lo", "obs_hi",
                                "sim_share", "sim_lo", "sim_hi", "sim_bound", "n_fallback_regime", "n_fallback_all")],
                          dplyr::across(dplyr::where(is.double), ~ round(.x, 3))))

  labels <- c(`1_pblhmod_default_vs_measured` = "1 PBLHmod vs Cloudnet MLH",
              `2_era5_center_vs_measured` = "2 ERA5 vs Cloudnet MLH")
  g <- ba_figure(lapply(prs, function(pt) pt[pt$regime %in% "convective", ]), cb, labels, "convective",
                 "PBLHmod convective; exact measured values only (biased, spec 07)")
  pth <- save_figure(g, paste0("10_bland_altman_", site_id), 20, 11)
  cli::cli_alert_success("{.path {pth[['png']]}}")
  cli::cli_h2("Compression statistics")
  print_tbl(dplyr::mutate(cb$comp[c(comp_cols, "mlr_T1", "mlr_T2", "mlr_T3")],
                          dplyr::across(dplyr::where(is.double), ~ round(.x, 3))))
  cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
  quit(save = "no")
}


# --- Aerosol reference: the primary source, plus the second SIRTA instrument at FR-Gri ---
# Short tags (pair names, figure labels): FR-Gri pair names unchanged from spec 10.
aero_tag <- c(sirta_mlh_chm15k = "chm15k", sirta_mlh_cl31 = "cl31", meteoswiss_stratfinder_mlh = "stratfinder",
              meteoswiss_cabam_mlh = "cabam")
aero_lab <- c(sirta_mlh_chm15k = "CHM15k", sirta_mlh_cl31 = "CL31", meteoswiss_stratfinder_mlh = "STRATfinder",
              meteoswiss_cabam_mlh = "CABAM")
ref1 <- ms$primary_source_id
ref2 <- if (ms$provider == "SIRTA-ReOBS") "sirta_mlh_cl31" else NULL
if (!all(c(ref1, ref2) %in% names(aero_tag))) cli::cli_abort("{site_id}: aerosol source {.val {ref1}} not supported.")
refs <- c(ref1, ref2)
meas <- lapply(stats::setNames(refs, refs), rd)
v_rng <- function(d) range(d$timestamp_start_utc[!is.na(d$pblh_m)])
p_start <- do.call(max, lapply(c(list(mod), meas), function(d) v_rng(d)[1]))
p_end <- do.call(min, lapply(c(list(mod), meas), function(d) v_rng(d)[2]))
if (p_start >= p_end) cli::cli_abort("{site_id}: no common period of PBLHmod and {paste(refs, collapse = ', ')}.")
cli::cli_inform("Comparison period (valid PBLHmod and {paste(aero_lab[refs], collapse = ', ')}): {format(p_start, '%Y-%m-%d %H:%M')} to {format(p_end, '%Y-%m-%d %H:%M')} UTC.")
in_per <- function(d) d[d$timestamp_start_utc >= p_start & d$timestamp_start_utc <= p_end, ]
meas <- lapply(meas, in_per)
tb <- lapply(meas, measured_time_base, mod = mod, era = era)
for (r in refs) {
  m <- meas[[r]]
  ok <- !is.na(m$pblh_m)
  pm <- tb[[r]]$mod$pblh_m[match(as.numeric(m$timestamp_start_utc[ok]), as.numeric(tb[[r]]$mod$timestamp_utc))]
  cli::cli_inform("{aero_lab[[r]]} ({c(halfhour = 'half-hourly', hour = 'hourly')[[tb[[r]]$base]]} intervals): {nrow(m)}, valid {sum(ok)} (qc 0 {sum(m$qc_flag %in% 0L)}, qc 1 {sum(m$qc_flag %in% 1L)}); PBLHmod valid on {sum(!is.na(pm))} of the valid ones.")
}

era_c_in <- dplyr::rename(in_per(era$era5_blh_center), timestamp_utc = timestamp_start_utc)
mod_h <- match_halfhours_to_hours(mod, era_c_in$timestamp_utc)
# Pairs 1 / 2 / 2b (primary), 3 / 4 / 4b (second instrument), 5 second vs primary, 6 PBLHmod vs ERA5 (record).
src_pairs <- function(r, k) {
  t <- tb[[r]]
  n <- c(paste0(k[1], "_pblhmod_default_vs_", aero_tag[[r]]), paste0(k[2], "_era5_center_vs_", aero_tag[[r]]),
         paste0(k[2], "b_era5_bilinear_vs_", aero_tag[[r]]))
  stats::setNames(list(
    list(x = t$mod, ref = exact_ref(t$meas), mid = t$mid_s, cls = "parameterisation vs aerosol", reference = r),
    list(x = t$era$era5_blh_center, ref = exact_ref(t$era_meas), mid = t$era_mid_s, cls = "cross-class: bulk Ri vs aerosol",
         reference = r),
    list(x = t$era$era5_blh_bilinear, ref = exact_ref(t$era_meas), mid = t$era_mid_s,
         cls = "cross-class: bulk Ri vs aerosol", reference = r)
  ), n)
}
pairs_def <- src_pairs(ref1, c("1", "2"))
if (!is.null(ref2)) {
  pairs_def <- c(pairs_def, src_pairs(ref2, c("3", "4")))
  # second instrument on the primary intervals, with the PBLHmod regime of the interval (pair 5)
  t1 <- tb[[ref1]]
  x5 <- dplyr::tibble(timestamp_utc = t1$meas$timestamp_utc,
                      pblh_m = meas[[ref2]]$pblh_m[match(as.numeric(t1$meas$timestamp_utc), as.numeric(meas[[ref2]]$timestamp_start_utc))],
                      x_regime = t1$mod$x_regime[match(as.numeric(t1$meas$timestamp_utc), as.numeric(t1$mod$timestamp_utc))])
  pairs_def[[paste0("5_", aero_tag[[ref2]], "_vs_", aero_tag[[ref1]])]] <-
    list(x = x5, ref = exact_ref(t1$meas), mid = t1$mid_s, reference = ref1,
         cls = "aerosol vs aerosol (CABAM vs STRATfinder): definitional uncertainty") # SIRTA CL31 vs CHM15k
}
pairs_def[["6_pblhmod_default_vs_era5_center"]] <- list(x = mod_h, ref = era_c_in, mid = 0, reference = "era5_blh_center",
                                                         cls = "parameterisation vs model (record)")
pairs <- lapply(names(pairs_def), function(p) add_time(pair_table(pairs_def[[p]]$x, pairs_def[[p]]$ref, p), pairs_def[[p]]$mid))
names(pairs) <- names(pairs_def)
cls <- vapply(pairs_def, `[[`, "", "cls")
ref_of <- vapply(pairs_def, `[[`, "", "reference")
p6 <- "6_pblhmod_default_vs_era5_center"


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
metrics$hour_note <- if (tb[[ref1]]$base == "hour") {
  "hour_utc = start hour of the interval [h, h + 1 h) (pairs 1-5) or ERA5 hour (pair 6)"
} else {
  "hour_utc = hour of the half-hour start (pair 1, half-hours [t, t + 30 min)) or ERA5 hour (pairs 2, 2b, 6)"
}
metrics$fb_convention <- fb_note
write_tab(metrics, file.path(tables_dir, paste0("10_metrics_", site_id, ".csv")))

# --- CL31 low-value cluster at FR-Gri (G 2026-10-29: report, keep the values) ---
p5n <- grep("^5_", names(pairs), value = TRUE)
if (length(p5n) == 1) {
  p5 <- pairs[[p5n]]
  p5w <- p5[sets$window_convective(p5), ]
  n_low <- sum(p5w$mod < 150 & p5w$ref > 300)
  cli::cli_inform("{aero_lab[[ref2]]} < 150 m while {aero_lab[[ref1]]} > 300 m: {n_low} of {nrow(p5w)} window convective intervals ({round(100 * n_low / nrow(p5w), 1)} %); all convective: {sum(p5$regime %in% 'convective' & p5$mod < 150 & p5$ref > 300)} of {sum(p5$regime %in% 'convective')}.")
  low <- p5w[p5w$mod < 150 & p5w$ref > 300, ]
  if (nrow(low) > 0) {
    print_tbl(dplyr::count(dplyr::mutate(low, year = format(timestamp_utc, "%Y")), year, season))
  }
}

# --- Compression ---
cli::cli_h2("Compression statistics")
note <- if (tb[[ref1]]$base == "hour") "hourly means, not censored" else
  "half-hourly values (ERA5 pairs: means of the two half-hours around the ERA5 hour), not censored"
comp <- list()
for (p in names(pairs)) {
  for (s in names(sets)) comp[[length(comp) + 1]] <- compress_rows(pairs[[p]][sets[[s]](pairs[[p]]), ], p, s, note)
}
cb <- bind_comp(comp)
cb$comp$reference <- unname(ref_of[cb$comp$pair])
cb$comp$comparison_class <- cls[cb$comp$pair]
comp <- cb$comp
write_tab(comp, file.path(tables_dir, paste0("10_compression_", site_id, ".csv")))
write_tab(cb$bins, file.path(tables_dir, paste0("10_bland_altman_bins_", site_id, ".csv")))

# --- Diurnal cycle (context; on the intervals with a valid primary value, in hour bins) ---
# ERA5: hourly interval mean (mean of t and t + 1 h) of the hour containing the interval.
t1 <- tb[[ref1]]
use <- t1$meas$timestamp_utc[!is.na(t1$meas$pblh_m)]
w <- dplyr::tibble(timestamp_utc = use)
for (r in refs) w[[r]] <- meas[[r]]$pblh_m[match(as.numeric(use), as.numeric(meas[[r]]$timestamp_start_utc))]
w$pblhmod_default <- t1$mod$pblh_m[match(as.numeric(use), as.numeric(t1$mod$timestamp_utc))]
hr <- as.POSIXct(floor(as.numeric(use) / 3600) * 3600, origin = "1970-01-01", tz = "UTC")
ei <- match_instants_to_intervals(era$era5_blh_center, hr) # full series: t + 1 h may fall after the period end
w$era5_blh_center <- ei$pblh_m[match(as.numeric(hr), as.numeric(ei$timestamp_utc))]
w$season <- season_of(w$timestamp_utc)
w$hour_utc <- as.integer(format(w$timestamp_utc, "%H", tz = "UTC"))
ds <- diurnal_summary(w, c(refs, "pblhmod_default", "era5_blh_center"))
ds$x_utc <- ds$hour_utc + 0.5 # midpoint of the hour bin
write_tab(ds, file.path(tables_dir, paste0("10_diurnal_", site_id, ".csv")))
# Window and night per season in UTC: window = AST limits minus the median AST - UTC offset; sunrise/sunset =
# median geometric times over the days of the season in the comparison period.
days <- seq(as.Date(p_start), as.Date(p_end), by = "day")
st <- sun_times(days, meta_site$lat, meta_site$lon)
h_of <- function(t) as.numeric(t - as.POSIXct(format(t, "%Y-%m-%d"), tz = "UTC"), units = "hours")
st$season <- season_of(as.POSIXct(st$date, tz = "UTC"))
tt <- as.POSIXct(days, tz = "UTC") + 12 * 3600
off <- dplyr::tibble(season = season_of(tt), d = apparent_solar_hour(tt, meta_site$lon) - 12)
bands <- dplyr::summarise(dplyr::group_by(st, season), sunrise = stats::median(h_of(sunrise_utc), na.rm = TRUE),
                          sunset = stats::median(h_of(sunset_utc), na.rm = TRUE), .groups = "drop")
bands <- dplyr::left_join(bands, dplyr::summarise(dplyr::group_by(off, season), d = stats::median(d), .groups = "drop"),
                          by = "season")
bands$win_start <- opts$aerosol_window[1] - bands$d
bands$win_end <- opts$aerosol_window[2] - bands$d

# --- Figures (ERA5 tower cell only) ---
cli::cli_h2("Figures")
fig_pairs <- c(stats::setNames(paste(c("1 PBLHmod vs", "2 ERA5 vs"), aero_lab[[ref1]]),
                               paste0(c("1_pblhmod_default_vs_", "2_era5_center_vs_"), aero_tag[[ref1]])))
if (!is.null(ref2)) {
  fig_pairs <- c(fig_pairs, stats::setNames(c(paste(c("3 PBLHmod vs", "4 ERA5 vs"), aero_lab[[ref2]]),
                                              paste("5", aero_lab[[ref2]], "vs", aero_lab[[ref1]])),
                                            c(paste0(c("3_pblhmod_default_vs_", "4_era5_center_vs_"), aero_tag[[ref2]]), p5n)))
}
# all rows of each pair: the primary set (window, convective) is the comparable set, the rest drawn faded (spec 11 E)
sc <- dplyr::bind_rows(lapply(names(fig_pairs), function(p) {
  x <- pairs[[p]]
  x$comparable <- sets$window_convective(x)
  x$panel <- fig_pairs[[p]]
  x
}))
cs <- comp[comp$pair %in% names(fig_pairs) & comp$stratum_set == "window_convective", ]
cs$panel <- unname(fig_pairs[cs$pair])
win_pairs <- lapply(pairs, function(pt) dplyr::mutate(pt, comparable = sets$window_convective(pt)))
# figure sizes [cm]: scatter width, height; Bland-Altman width, height (FR-Gri: five pairs; one source: two)
nr <- if (length(fig_pairs) > 2) c(24, 17, 24, 15) else c(14, 19, 18, 12)
figs <- list(
  list(plot_aerosol_diurnal(ds[c("season", "x_utc", "variable", "median")], bands), "10_diurnal_by_season", 17, 13),
  list(plot_compression_scatter(sc, cs), "10_scatter", nr[1], nr[2]),
  list(ba_figure(win_pairs, cb, fig_pairs, "window_convective",
                 sprintf("Window %g-%g h apparent solar time, PBLHmod convective; %s",
                         opts$aerosol_window[1], opts$aerosol_window[2],
                         if (tb[[ref1]]$base == "hour") "hourly means" else "half-hourly values")), "10_bland_altman", nr[3], nr[4])
)
for (fg in figs) {
  pth <- save_figure(fg[[1]], paste0(fg[[2]], "_", site_id), fg[[3]], fg[[4]])
  cli::cli_alert_success("{.path {pth[['png']]}}")
}

# --- Report ---
get_m <- function(p, set, strata, mets = c("n", "bias", "mae", "median_log_ratio", "mad_log_ratio",
                                           "frac_within_factor2", "r_spearman", "fb", "nrmsd")) {
  x <- metrics[metrics$pair == p & metrics$stratum_set == set & metrics$stratum %in% strata & metrics$metric %in% mets, ]
  w <- stats::reshape(as.data.frame(x[c("stratum", "metric", "value")]), idvar = "stratum", timevar = "metric",
                      direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w <- w[match(strata, w$stratum), c("stratum", intersect(mets, names(w)))]
  w[!is.na(w$stratum), ]
}
cli::cli_inform("FB convention: {fb_note}.")
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
cli::cli_h2("Compression: symmetric slope and Bland-Altman trend (day-block bootstrap 95 %); Theil-Sen diagnostics")
print_tbl(dplyr::mutate(comp[comp_cols], dplyr::across(dplyr::where(is.double), ~ round(.x, 3))))
cli::cli_h2("Diurnal medians [m], JJA and DJF (hour bins of the intervals with a {aero_lab[[ref1]]} value; night = context only)")
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
