# Spec 07 section 5: PBLHmod, ERA5 and the operational ICOS PBLH against the measured turbulence-based MLH
# (Cloudnet Doppler-lidar epsilon, source `cloudnet_eps_mlh`), for one station.
# Inputs: standard files in <standard_dir>/<site>/ written by scripts 04 (pblhmod_default, era5_blh_center,
# era5_blh_bilinear, icos_oper) and 07c (cloudnet_eps_mlh).
# Measured MLH (G 2026-10-23): qc 0/2 exact (metrics); qc 1 censored below z_min (share of model values below
# z_min reported); qc 3 signal-limited = lower bound (share of model values below the bound, median model/bound).
# Pairs (model vs measured): 1 pblhmod_default (half-hours); 2 / 2b era5_blh_center / era5_blh_bilinear (ERA5
# hours, measured = mean of the two adjacent half-hours, both required; cross-class: bulk Richardson vs turbulence);
# 3a / 3b icos_oper qc 0 / all (half-hours); 4 pblhmod_default vs era5_blh_center over the measured period (record).
# Strata as spec 04 (PBLHmod regime, season, hour UTC); primary statistics = convective regime; stable and
# transition hours separately. The distance class of the pairing (config/measured_sites.csv) limits the regimes.
# Output: <output_dir>/tables/07_metrics_<site>.csv (long), 07_lower_bound_<site>.csv, 07_censoring_<site>.csv,
#         07_diurnal_<site>.csv; figures 07_diurnal_by_season_, 07_scatter_, 07_censoring_, 07_logratio_by_hour_
#         <site>.png.
#
# Usage: Rscript scripts/07_measured_compare.R SITE
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
site_id <- commandArgs(trailingOnly = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val FI-Hyy}.")
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$provider == "cloudnet" & ms$icos_site == site_id, ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) Cloudnet pairing in {.file measured_sites.csv}.")
dclass <- distance_class(ms$distance_km)
if (dclass == "not_used") cli::cli_abort("{site_id}: measurement {ms$distance_km} km away: not used (distance rule).")
std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)
rd <- function(src) {
  p <- file.path(std_dir, paste0(site_id, "__", src, ".csv.gz"))
  if (!file.exists(p)) cli::cli_abort("Standard file {.path {p}} not found (run scripts 04 and 07c).")
  read_standard_pblh(p)
}

cli::cli_h1("Measured MLH comparison: {site_id} ({ms$provider_site}, {ms$distance_km} km, {dclass})")
meas <- rd("cloudnet_eps_mlh")
mod <- rd("pblhmod_default")
era_c <- rd("era5_blh_center")
era_b <- rd("era5_blh_bilinear")
oper <- rd("icos_oper")
t_rng <- range(meas$timestamp_start_utc)
cli::cli_inform("Measured period: {format(t_rng[1], '%Y-%m-%d %H:%M')} to {format(t_rng[2], '%Y-%m-%d %H:%M')} UTC; {nrow(meas)} half-hours with profiles.")
in_per <- function(d) d[d$timestamp_start_utc >= t_rng[1] & d$timestamp_start_utc <= t_rng[2], ]
mod <- in_per(mod)
oper <- in_per(oper)
era_c <- in_per(era_c)
era_b <- in_per(era_b)

# --- Measured, half-hourly and hourly ---
mh <- dplyr::bind_cols(
  dplyr::tibble(timestamp_utc = meas$timestamp_start_utc, pblh_m = meas$pblh_m, qc_flag = meas$qc_flag,
                z_min = meas$x_censor_below_m),
  measured_interval(meas$pblh_m, meas$qc_flag, meas$x_censor_below_m)[c("left", "right", "kind")]
)
hours <- era_c$timestamp_start_utc
mhr <- match_measured_to_hours(meas, hours)
mhr <- dplyr::bind_cols(mhr, measured_interval(mhr$pblh_m, mhr$qc_flag, mhr$x_censor_below_m)[c("left", "right", "kind")])
mhr$z_min <- mhr$x_censor_below_m
cli::cli_inform(c(
  "Half-hours: exact {sum(mh$kind %in% 'exact')}, below z_min {sum(mh$kind %in% 'below')}, lower bound {sum(mh$kind %in% 'lower_bound')}, not usable {sum(is.na(mh$kind))}.",
  "ERA5 hours: exact {sum(mhr$kind %in% 'exact')}, below z_min {sum(mhr$kind %in% 'below')}, lower bound {sum(mhr$kind %in% 'lower_bound')}, not usable {sum(is.na(mhr$kind))}."
))

# --- Model / reanalysis / operational series on the measurement time base ---
mod_hh <- dplyr::rename(mod, timestamp_utc = timestamp_start_utc)
reg_hh <- stats::setNames(mod_hh$x_regime, as.character(as.numeric(mod_hh$timestamp_utc)))
oper_hh <- dplyr::rename(oper, timestamp_utc = timestamp_start_utc)
oper_hh$x_regime <- reg_hh[as.character(as.numeric(oper_hh$timestamp_utc))]
oper_qc0 <- oper_hh[oper_hh$qc_flag %in% 0L, ]
mod_h <- match_halfhours_to_hours(mod, hours)
reg_h <- stats::setNames(mod_h$x_regime, as.character(as.numeric(mod_h$timestamp_utc)))
as_h <- function(e) {
  e <- dplyr::rename(e, timestamp_utc = timestamp_start_utc)
  e$x_regime <- reg_h[as.character(as.numeric(e$timestamp_utc))]
  e
}
era_ch <- as_h(era_c)
era_bh <- as_h(era_b)
sources <- list(
  `1_pblhmod_default_vs_measured` = list(x = mod_hh, m = mh, src = "PBLHmod", cls = "parameterisation vs turbulence"),
  `2_era5_center_vs_measured` = list(x = era_ch, m = mhr, src = "ERA5 (tower cell)", cls = "cross-class: bulk Ri vs turbulence"),
  `2b_era5_bilinear_vs_measured` = list(x = era_bh, m = mhr, src = "ERA5 (bilinear)", cls = "cross-class: bulk Ri vs turbulence"),
  `3a_icos_oper_qc0_vs_measured` = list(x = oper_qc0, m = mh, src = "ICOS oper. (qc 0)", cls = "parameterisation vs turbulence"),
  `3b_icos_oper_all_vs_measured` = list(x = oper_hh, m = mh, src = "ICOS oper. (all)", cls = "parameterisation vs turbulence")
)
# Distance rule: restrict the regimes compared (spec 07, project plan).
keep_regime <- function(r) {
  if (dclass == "all_regimes") rep(TRUE, length(r)) else r %in% "convective"
}
if (dclass == "daytime_climatology") cli::cli_alert_warning("Distance class {dclass}: pair metrics are indicative only.")

# --- Pairs on exact values (qc 0/2) ---
pairs <- dplyr::bind_rows(lapply(names(sources), function(p) {
  s <- sources[[p]]
  ref <- s$m[s$m$kind %in% "exact", c("timestamp_utc", "pblh_m")]
  pt <- pair_table(s$x, ref, p)
  pt[keep_regime(pt$regime), ]
}))
pairs$source <- vapply(sources, `[[`, "", "src")[pairs$pair]
p4 <- pair_table(mod_h, era_ch, "4_pblhmod_default_vs_era5_center")
metrics <- dplyr::bind_rows(pblh_metrics_long(pairs), pblh_metrics_long(p4))
cls <- c(vapply(sources, `[[`, "", "cls"), `4_pblhmod_default_vs_era5_center` = "parameterisation vs model")
metrics$comparison_class <- cls[metrics$pair]
metrics$distance_class <- dclass
readr::write_csv(metrics, file.path(tables_dir, paste0("07_metrics_", site_id, ".csv")))

# --- Lower bounds (qc 3) and censored below z_min (qc 1) ---
lb <- list()
cb <- list()
for (p in names(sources)) {
  s <- sources[[p]]
  x <- dplyr::left_join(s$m, dplyr::transmute(s$x, timestamp_utc, mod = pblh_m, regime = x_regime), by = "timestamp_utc")
  x$season <- season_of(x$timestamp_utc)
  x$hour_utc <- as.integer(format(x$timestamp_utc, "%H", tz = "UTC"))
  x <- x[keep_regime(x$regime), ]
  xl <- x[x$kind %in% "lower_bound", ]
  xl$bound <- xl$pblh_m
  lb[[p]] <- dplyr::bind_rows(
    dplyr::bind_cols(dplyr::tibble(pair = p, regime = "all", season = "all"), lower_bound_summary(xl)),
    dplyr::bind_cols(dplyr::tibble(pair = p, season = "all"), lower_bound_summary(xl, "regime")),
    dplyr::bind_cols(dplyr::tibble(pair = p, regime = "all"), dplyr::mutate(lower_bound_summary(xl, "season"), season = as.character(season))),
    dplyr::bind_cols(dplyr::tibble(pair = p), dplyr::mutate(lower_bound_summary(xl, c("regime", "season")), season = as.character(season)))
  )
  xc <- x[x$kind %in% "below", ]
  cb[[p]] <- dplyr::bind_rows(
    dplyr::bind_cols(dplyr::tibble(pair = p, season = "all", hour_utc = NA_integer_), censored_below_summary(xc)),
    dplyr::bind_cols(dplyr::tibble(pair = p), dplyr::mutate(censored_below_summary(xc, c("season", "hour_utc")),
                                                           season = as.character(season)))
  )
}
lb <- dplyr::bind_rows(lb)
cb <- dplyr::bind_rows(cb)
readr::write_csv(lb, file.path(tables_dir, paste0("07_lower_bound_", site_id, ".csv")))
cens_meas <- dplyr::summarise(
  dplyr::group_by(dplyr::mutate(mh, season = as.character(season_of(timestamp_utc)),
                                hour_utc = as.integer(format(timestamp_utc, "%H", tz = "UTC"))), season, hour_utc),
  n_usable = sum(!is.na(kind)), n_below = sum(kind %in% "below"), n_lower_bound = sum(kind %in% "lower_bound"),
  pct_below = round(100 * n_below / max(1, n_usable), 1), pct_lower_bound = round(100 * n_lower_bound / max(1, n_usable), 1),
  .groups = "drop"
)
w <- cb[!is.na(cb$hour_utc), c("pair", "season", "hour_utc", "frac_mod_below_zmin")]
w$frac_mod_below_zmin <- round(w$frac_mod_below_zmin, 3)
w <- stats::reshape(as.data.frame(w), idvar = c("season", "hour_utc"), timevar = "pair", direction = "wide")
cens_tab <- dplyr::left_join(cens_meas, w, by = c("season", "hour_utc"))
readr::write_csv(cens_tab, file.path(tables_dir, paste0("07_censoring_", site_id, ".csv")))

# --- Diurnal cycles (on the measured usable half-hours / hours) ---
mh$season <- season_of(mh$timestamp_utc)
mh$hour_utc <- as.integer(format(mh$timestamp_utc, "%H", tz = "UTC"))
dm <- measured_diurnal(mh)
use_hh <- mh$timestamp_utc[!is.na(mh$kind)]
use_h <- mhr$timestamp_utc[!is.na(mhr$kind)]
w_hh <- dplyr::tibble(timestamp_utc = use_hh)
w_hh$pblhmod_default <- mod_hh$pblh_m[match(use_hh, mod_hh$timestamp_utc)]
w_hh$icos_oper_qc0 <- oper_qc0$pblh_m[match(use_hh, oper_qc0$timestamp_utc)]
w_h <- dplyr::tibble(timestamp_utc = use_h)
w_h$era5_blh_center <- era_ch$pblh_m[match(use_h, era_ch$timestamp_utc)]
w_h$era5_blh_bilinear <- era_bh$pblh_m[match(use_h, era_bh$timestamp_utc)]
for (nm in c("w_hh", "w_h")) {
  w <- get(nm)
  w$season <- season_of(w$timestamp_utc)
  w$hour_utc <- as.integer(format(w$timestamp_utc, "%H", tz = "UTC"))
  assign(nm, w)
}
ds <- dplyr::bind_rows(diurnal_summary(w_hh, c("pblhmod_default", "icos_oper_qc0")),
                       diurnal_summary(w_h, c("era5_blh_center", "era5_blh_bilinear")))
readr::write_csv(dplyr::bind_rows(
  dplyr::transmute(dm, season, hour_utc, variable = "cloudnet_eps_mlh_censored", median = median_censored, n),
  dplyr::transmute(dm, season, hour_utc, variable = "cloudnet_eps_mlh_qc02", median = median_exact, n = n_exact),
  ds[c("season", "hour_utc", "variable", "median", "n")]
), file.path(tables_dir, paste0("07_diurnal_", site_id, ".csv")))

# --- Figures ---
cli::cli_h2("Figures")
pairs_fig <- pairs[pairs$pair %in% c("1_pblhmod_default_vs_measured", "2_era5_center_vs_measured",
                                     "3a_icos_oper_qc0_vs_measured"), ]
figs <- list(
  list(plot_measured_diurnal(dm, ds[ds$variable != "era5_blh_bilinear", ]), "07_diurnal_by_season", 17, 13),
  list(plot_measured_scatter(pairs_fig), "07_scatter", 24, 18),
  list(plot_measured_qc(meas), "07_censoring", 17, 13),
  list(plot_measured_logratio(pairs_fig), "07_logratio_by_hour", 20, 13)
)
for (fg in figs) {
  pth <- save_figure(fg[[1]], paste0(fg[[2]], "_", site_id), fg[[3]], fg[[4]])
  cli::cli_alert_success("{.path {pth[['png']]}}")
}

# --- Report ---
get_m <- function(p, strata, mets = c("n", "bias", "mae", "median_log_ratio", "mad_log_ratio", "frac_within_factor2",
                                      "r_spearman")) {
  x <- metrics[metrics$pair == p & metrics$stratum %in% strata & metrics$metric %in% mets, ]
  w <- stats::reshape(as.data.frame(x[c("stratum", "metric", "value")]), idvar = "stratum", timevar = "metric",
                      direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w <- w[match(strata, w$stratum), c("stratum", intersect(mets, names(w)))]
  w[!is.na(w$stratum), ]
}
cli::cli_h2("Primary: convective regime (qc 0/2)")
prim <- dplyr::bind_rows(lapply(c(names(sources), "4_pblhmod_default_vs_era5_center"), function(p) {
  dplyr::bind_cols(dplyr::tibble(pair = p), get_m(p, "regime=convective"))
}))
print_tbl(prim)
for (p in c(names(sources), "4_pblhmod_default_vs_era5_center")) {
  cli::cli_h3("{p} [{cls[[p]]}]")
  print_tbl(get_m(p, c("all", "regime=convective", "regime=unstable_diag", "regime=mixed", "regime=stable",
                       "season=DJF", "season=MAM", "season=JJA", "season=SON")))
}
cli::cli_h2("Lower bounds (measured qc 3): model below the bound, median model / bound")
print_tbl(lb[lb$season == "all" | lb$regime == "all", ])
cli::cli_h2("Censored below z_min (measured qc 1): model value <= z_min")
print_tbl(cb[cb$season == "all", ])
print_tbl(dplyr::summarise(dplyr::group_by(cb[!is.na(cb$hour_utc), ], pair, season),
                           frac_mod_below_zmin = round(sum(frac_mod_below_zmin * n) / sum(n), 3), n = sum(n),
                           .groups = "drop"))
cli::cli_h2("Diurnal medians [m], JJA and DJF (measured censored / qc 0-2, PBLHmod, ERA5 centre, operational qc 0)")
for (s in c("JJA", "DJF")) {
  cat("--", s, "\n")
  t <- dm[dm$season == s, c("hour_utc", "median_censored", "median_exact", "n_lower_bound")]
  for (v in c("pblhmod_default", "era5_blh_center", "icos_oper_qc0")) {
    x <- ds[ds$season == s & ds$variable == v, ]
    t[[v]] <- x$median[match(t$hour_utc, x$hour_utc)]
  }
  print_tbl(dplyr::mutate(t, dplyr::across(-hour_utc, round)))
}
