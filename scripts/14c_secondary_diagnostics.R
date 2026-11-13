# Spec 14 Part C: diagnostics of the secondary measured sources (internal; not for the paper).
#   DE-Tha (tud_skyvue_mlh, layers 1-3) and CZ-Lnz (vaisala_blview_mlh, all versions):
#   1. Time-convention test: on clear days (kt >= 0.65, SW_IN on >= 90 % of the daylight half-hours; spec 10), the
#      lag L (-3 ... +3 h) that maximises the Pearson correlation between the reference (PBLHmod, ERA5) in the
#      morning hours (sunrise to solar noon) and the measured value L hours later, per season and per series (DE-Tha
#      layer, BL-View version). L = 0 supports UTC; +1 / +2 means stamps in CET / CEST. Reported, nothing corrected.
#      Plus the median clear-day diurnal cycle by UTC hour.
#   2. DE-Tha layers: per UTC hour, median of layers 1-3 and of their confidence vs PBLHmod and ERA5 (all days), and
#      window metrics per layer (aerosol window 10-16 h apparent solar time, PBLHmod convective), so that G can
#      decide a layer rule (or none).
#   3. CZ-Lnz quality: window metrics of PBLHmod and ERA5 against the BL-View value by quality set and version.
#      fix-15 (G 2026-11-12): main set = quality >= quality_min of measured_sites.csv (CZ-Lnz 2; the instrument's own
#      reliability index, not the agreement with the references); sensitivity rows all and q3; classes q1, q2 kept
#      (quality_role). Lag test and diurnal cycle on the main set, plus blview_allq and blview_q3.
# Inputs: <standard_dir>/<site>/<site>__{tud_skyvue_mlh|vaisala_blview_mlh, pblhmod_default, era5_blh_center}
#   (scripts 14a / 14b, 04), <derived_dir>/pblhmod/<site>_pblhmod_default.csv.gz (day_id), the ICOS L2 (SW_IN,
#   SW_IN_POT) for the clear days.
# Output (internal): <output_dir>/tables/internal/14_lag_<site>.csv (series x reference x season x lag),
#   14_lag_best_<site>.csv, 14_diurnal_clear_<site>.csv, 14_window_metrics_<site>.csv, DE-Tha also
#   14_skyvue_layers_hourly_<site>.csv; <output_dir>/figures/internal/14_lag_<site>.png, 14_diurnal_clear_<site>.png,
#   DE-Tha 14_skyvue_layers_hourly_<site>.png.
#
# fix-14c (2026-11-12): lag table without data masking (one row per site x series x reference x season in
#   14_lag_best), lags -3 ... +6 h, references named (PBLHmod, ERA5), window metrics per SkyVUE layer
#   (secondary_window_metrics()), group aesthetic of the lag figure; --blview-322-tz=Europe/Prague re-reads the 3.2.2
#   stamps as local time with summer time (hypothesis of G 2026-11-12; diagnostics only, outputs with suffix _tz322;
#   the standard file is not changed).
# fix-14a item 11: 14b now writes the BL-View standard file time-corrected by default (3.2.2 Europe/Prague -> UTC,
#   all versions -1 h; x_time_shift_h); the lag test is then expected at 0 h in every season with enough data (a
#   check, not a target). --blview-322-tz is refused on a corrected file (only for a file with stored stamps).
#
# Usage: Rscript scripts/14c_secondary_diagnostics.R SITE [--from=YYYY-MM-DD] [--to=YYYY-MM-DD]
#          [--blview-322-tz=Europe/Prague]
#   e.g. DE-Tha --from=2026-07-01 --to=2026-07-31 (test month), then the full record without --from / --to.
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a) == 1) sub(paste0("^--", name, "="), "", a) else default
}
pos <- args[!startsWith(args, "--")]
bad <- setdiff(sub("=.*$", "", args[startsWith(args, "--")]), c("--from", "--to", "--blview-322-tz"))
if (length(bad) > 0) cli::cli_abort("Unknown option{?s} {.val {bad}}.")
if (length(pos) != 1) cli::cli_abort("Pass one ICOS station code (DE-Tha or CZ-Lnz).")
site_id <- pos[1]
cfg <- get_cfg()
t0 <- Sys.time()
opts <- measured_compare_options()
kt_clear <- 0.65 # clear day (spec 10 item 3b, G 2026-10-30)
kt_cov <- 0.9
sites <- read_sites(resolve_path(cfg$sites_file))
meta <- sites[sites$site_id == site_id, ]
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
# the row by provider: the SkyVUE MLH is inventory since G 2026-11-12, its diagnostics stay available
ms <- ms[ms$icos_site == site_id & ms$provider %in% c("TU Dresden", "CzechGlobe"), ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) SkyVUE / BL-View row in {.file measured_sites.csv}.")
src <- c(`TU Dresden` = "tud_skyvue_mlh", CzechGlobe = "vaisala_blview_mlh")[[ms$provider]]
cli::cli_inform("{site_id}: {src}, status {ms$status} in {.file measured_sites.csv}.")
tz322 <- opt("blview-322-tz")
if (!is.null(tz322) && src != "vaisala_blview_mlh") cli::cli_abort("--blview-322-tz applies to BL-View (CZ-Lnz) only.")
if (!is.null(tz322) && !tz322 %in% OlsonNames()) cli::cli_abort("Unknown time zone {.val {tz322}}.")
sfx <- if (is.null(tz322)) "" else "_tz322" # outputs of the conversion run kept apart
std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
rd <- function(s) read_standard_pblh(file.path(std_dir, paste0(site_id, "__", s, ".csv.gz")))
tab_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables", "internal")
fig_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "figures", "internal")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
out_tab <- function(x, name) readr::write_csv(x, file.path(tab_dir, paste0("14_", name, "_", site_id, sfx, ".csv")))
print_tbl <- function(x) print(as.data.frame(dplyr::mutate(x, dplyr::across(dplyr::where(is.double), ~ signif(.x, 3)))),
                               row.names = FALSE)

meas <- rd(src)
t_from <- if (is.null(opt("from"))) min(meas$timestamp_start_utc) else as.POSIXct(opt("from"), tz = "UTC")
t_to <- if (is.null(opt("to"))) max(meas$timestamp_end_utc) else as.POSIXct(opt("to"), tz = "UTC") + 86400
in_per <- function(d, col = "timestamp_start_utc") d[d[[col]] >= t_from & d[[col]] < t_to, ]
meas <- in_per(meas)
mod <- in_per(rd("pblhmod_default"))
era <- list(era5_blh_center = in_per(rd("era5_blh_center")))
if (nrow(meas) == 0) cli::cli_abort("{site_id}: no {src} row between {format(t_from)} and {format(t_to)}.")
if (src == "vaisala_blview_mlh") {
  # fix-14a item 11: the standard file is time-corrected by 14b by default (x_time_shift_h); report it, and refuse a
  # second 3.2.2 conversion
  corrected <- "x_time_shift_h" %in% names(meas) && any(!is.na(meas$x_time_shift_h))
  if (corrected) {
    sh <- unique(meas[c("x_version", "x_time_shift_h")])
    cli::cli_alert_info("BL-View standard file time-corrected by 14b (fix-14a item 11); shifts by version (h): {paste(sh$x_version, sh$x_time_shift_h, collapse = '; ')}. Expected lag: 0 h.")
  } else {
    cli::cli_alert_warning("BL-View standard file with stored stamps (written before fix-14a item 11, or with --blview-322-tz=none --shift-h=0).")
  }
  if (corrected && !is.null(tz322)) cli::cli_abort("--blview-322-tz: the standard file is already time-corrected (x_time_shift_h); rerun without it.")
}
if (!is.null(tz322)) {
  # fix-14c (G 2026-11-12): hypothesis = 3.2.2 stamps are local time of tz322 (with summer time); diagnostics only
  cv <- blview_convert_version(meas, "3.2.2", tz322)
  meas <- cv$data
  cli::cli_alert_info("BL-View 3.2.2 stamps re-read as {tz322} local time: {cv$n_dropped} nonexistent hour{?s} dropped, {cv$n_duplicated} duplicated hour{?s} (repeated autumn hour) dropped; outputs with suffix {.val {sfx}}.")
}
cli::cli_h1("{site_id} {src}: {format(t_from, '%Y-%m-%d')} to {format(t_to - 1, '%Y-%m-%d')} ({nrow(meas)} rows)")
cli::cli_inform("PBLHmod rows {nrow(mod)} (valid {sum(!is.na(mod$pblh_m))}); ERA5 hours {nrow(era$era5_blh_center)}.")

# series of the measured source: SkyVUE layers (half-hourly), BL-View all versions and each version (hourly)
# fix-15 (G 2026-11-12): BL-View main set = quality >= quality_min of measured_sites.csv (CZ-Lnz 2). Lag test and
# diurnal cycle: blview_all and the versions on the main set, plus the sensitivity series blview_allq (every class)
# and blview_q3 (class 3 only). Window metrics: unfiltered series, quality sets inside secondary_window_metrics().
q_min <- ms$quality_min
series_wm <- if (src == "tud_skyvue_mlh") {
  stats::setNames(lapply(1:3, function(k) dplyr::mutate(meas, pblh_m = .data[[paste0("x_layer", k, "_m")]])),
                  paste0("skyvue_layer", 1:3))
} else {
  v <- sort(unique(meas$x_version))
  c(list(blview_all = meas), stats::setNames(lapply(v, function(x) meas[meas$x_version %in% x, ]), paste0("blview_", v)))
}
series <- series_wm
if (src == "vaisala_blview_mlh" && !is.na(q_min)) {
  keep_q <- function(d, qs) d[d$x_quality %in% qs, ]
  series <- c(lapply(series_wm, keep_q, qs = seq.int(q_min, 3L)),
              list(blview_allq = meas, blview_q3 = keep_q(meas, 3L)))
  cli::cli_inform("BL-View main set: quality >= {q_min} ({sum(meas$x_quality %in% seq.int(q_min, 3L) & !is.na(meas$pblh_m))} of {sum(!is.na(meas$pblh_m))} valid hours); sensitivity series blview_allq (every class) and blview_q3.")
}

# --- clear days ---
pm_raw <- readr::read_csv(file.path(resolve_path(cfg$derived_dir), "pblhmod", paste0(site_id, "_pblhmod_default.csv.gz")),
                          col_select = c("timestamp_start_utc", "day_id"), show_col_types = FALSE, progress = FALSE)
icos <- suppressMessages(suppressWarnings(read_icos_site(site_id, cfg, sites,
                                                         read_sites_timevarying(resolve_path(cfg$sites_timevarying_file)),
                                                         read_exclusions(resolve_path(cfg$exclusions_file)))))
icos <- in_per(icos)
day_of <- as.Date(pm_raw$day_id)[match(as.numeric(icos$timestamp_start_utc), as.numeric(as_utc_time(pm_raw$timestamp_start_utc)))]
kt <- daily_clearness(day_of, icos$met_SW_IN, icos$fn_SW_IN_POT)
clear <- kt$day_id[!is.na(kt$kt) & kt$coverage >= kt_cov & kt$kt >= kt_clear]
cli::cli_inform("Clear days (kt >= {kt_clear}, coverage >= {kt_cov}): {length(clear)} of {nrow(kt)}.")
mb <- morning_bins(clear, meta$lat, meta$lon)

# --- 1. time-convention test ---
ref_bins <- list(pblhmod_default = hour_bins(mod$timestamp_start_utc, mod$timestamp_end_utc, mod$pblh_m),
                 era5_blh_center = instant_hour_bins(era$era5_blh_center$timestamp_start_utc, era$era5_blh_center$pblh_m))
meas_bins <- lapply(series, function(s) hour_bins(s$timestamp_start_utc, s$timestamp_end_utc, s$pblh_m))
lags <- -3:6 # fix-14c: widened (the 3.2.2 peak sat at the +3 h edge)
lag <- lag_table(meas_bins, ref_bins, mb, site_id, lags = lags)
best <- best_lag(lag, by = c("site_id", "series", "reference"))
# fix-15: quality set of each series (BL-View with quality_min: main / sensitivity; otherwise all)
qset <- function(x) if (src == "vaisala_blview_mlh" && !is.na(q_min)) {
  ifelse(x %in% c("blview_allq", "blview_q3"), paste0("sensitivity (", sub("blview_", "", x), ")"), paste0("main (q >= ", q_min, ")"))
} else rep("all", length(x))
lag$quality_set <- qset(lag$series)
best$quality_set <- qset(best$series)
out_tab(lag, "lag")
out_tab(best, "lag_best")
cli::cli_h2("Best lag of the morning growth (L > 0: measured stamps later than UTC, e.g. +1 = CET); lags {min(lags)} to +{max(lags)} h")
if (nrow(best) > 0) {
  bp <- best
  bp$reference <- reference_label(bp$reference)
  bp$season <- factor(bp$season, levels = c("all", "DJF", "MAM", "JJA", "SON"))
  bp <- bp[order(bp$series, bp$reference, bp$season), ]
  bp$season <- as.character(bp$season)
  bp$at_edge <- bp$best_lag_h %in% range(lags)
  print_tbl(bp[c("series", "quality_set", "reference", "season", "best_lag_h", "r_best", "r_lag0", "n_best", "at_edge")])
} else {
  cli::cli_alert_warning("No correlation (too few clear-day morning pairs).")
}

clear_bins <- function(b) b[as.Date(b$bin_utc + 1800) %in% clear | as.Date(b$bin_utc) %in% clear, ]
diur <- dplyr::bind_rows(lapply(c(names(series), names(ref_bins)), function(s) {
  b <- clear_bins(if (s %in% names(series)) meas_bins[[s]] else ref_bins[[s]])
  if (nrow(b) == 0) return(NULL)
  b$season <- as.character(season_of(b$bin_utc))
  b$hour_utc <- as.integer(format(b$bin_utc, "%H", tz = "UTC"))
  dplyr::summarise(dplyr::group_by(b, season, hour_utc), median_m = stats::median(value), n = dplyr::n(),
                   .groups = "drop") |> dplyr::mutate(series = s, .before = 1)
}))
out_tab(diur, "diurnal_clear")
sty <- function(ids) {
  grp <- ifelse(ids %in% style_ids("secondary"), "secondary", "source")
  list(col = vapply(seq_along(ids), function(i) unname(style_of(ids[i], "colour", grp[i])), ""),
       lab = vapply(seq_along(ids), function(i) unname(style_of(ids[i], "label", grp[i])), ""),
       lty = vapply(seq_along(ids), function(i) unname(style_of(ids[i], "linetype", grp[i])), ""))
}
if (nrow(diur) > 0) {
  ids <- unique(diur$series)
  st <- sty(ids)
  g <- ggplot2::ggplot(diur, ggplot2::aes(hour_utc + 0.5, median_m, colour = series, linetype = series)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~season) +
    ggplot2::scale_colour_manual(values = stats::setNames(st$col, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(st$lty, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::labs(x = "UTC hour (bin midpoint)", y = "Median height [m]",
                  subtitle = paste0(site_id, ", clear days (kt >= ", kt_clear, "); internal")) +
    theme_pblh()
  p <- save_figure(g, paste0("14_diurnal_clear_", site_id, sfx), 19, 14, dir = fig_dir)
  cli::cli_alert_success("{.path {p[['png']]}}")
}
lagp <- lag[!is.na(lag$r), ]
if (nrow(lagp) > 0) {
  ids <- unique(lagp$series)
  st <- sty(ids)
  lagp$reference <- reference_label(lagp$reference)
  lagp$season <- factor(lagp$season, levels = c("all", "DJF", "MAM", "JJA", "SON"))
  g <- ggplot2::ggplot(lagp, ggplot2::aes(lag_h, r, colour = series, linetype = series, group = series)) +
    ggplot2::geom_vline(xintercept = 0, colour = style_of("reference_line", "colour", "element"), linewidth = 0.3) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::facet_grid(reference ~ season) +
    ggplot2::scale_colour_manual(values = stats::setNames(st$col, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(st$lty, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::scale_x_continuous(breaks = lags) +
    ggplot2::labs(x = "Lag L [h] (measured value L hours after the reference hour)", y = "Pearson r, morning growth",
                  subtitle = paste0(site_id, ", clear-day mornings (sunrise to solar noon); L = 0: UTC, +1: CET, +2: CEST")) +
    theme_pblh()
  p <- save_figure(g, paste0("14_lag_", site_id, sfx), 19, 12, dir = fig_dir)
  cli::cli_alert_success("{.path {p[['png']]}}")
}

# --- window metrics (aerosol window, PBLHmod convective) per series and, for BL-View, per quality class ---
# SkyVUE: one series per layer (pblh_m = the layer); BL-View: per quality class (secondary_window_metrics())
wm <- dplyr::bind_rows(lapply(names(series_wm), function(s) {
  secondary_window_metrics(series_wm[[s]], s, mod, era, meta$lon, opts, site_id, quality_min = q_min)
}))
if (nrow(wm) == 0) cli::cli_alert_warning("No window metrics: no measured value in the period.")
if (nrow(wm) > 0) {
  wm$note <- paste("ln(source / measured); window", paste(opts$aerosol_window, collapse = "-"),
                   "h apparent solar time, PBLHmod convective; BL-View quality sets (quality_role: main = quality >=",
                   q_min, ", sensitivity = all / q3, class = single classes); SkyVUE: one series per layer")
}
out_tab(wm, "window_metrics")
cli::cli_h2("Window metrics (source vs measured; SkyVUE per layer, BL-View per quality class)")
if (nrow(wm) > 0) {
  wp <- wm
  wp$source <- reference_label(wp$source)
  print_tbl(wp[c("series", "source", "quality", "quality_role", "n", "n_days", "median_log_ratio", "frac_within_factor2")])
}

# --- 2. DE-Tha layers by UTC hour ---
if (src == "tud_skyvue_mlh") {
  hr <- function(t) as.integer(format(t, "%H", tz = "UTC"))
  lay <- dplyr::bind_rows(lapply(1:3, function(k) {
    x <- meas[!is.na(meas[[paste0("x_layer", k, "_m")]]), ]
    dplyr::summarise(dplyr::group_by(dplyr::tibble(hour_utc = hr(x$timestamp_start_utc), h = x[[paste0("x_layer", k, "_m")]],
                                                   cf = x[[paste0("x_layer", k, "_conf")]], c3 = x[[paste0("x_layer", k, "_share_conf3")]]),
                                     hour_utc),
                     median_m = stats::median(h), median_conf = stats::median(cf, na.rm = TRUE),
                     mean_share_conf3 = mean(c3, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(series = paste0("skyvue_layer", k), .before = 1)
  }))
  ref_h <- dplyr::bind_rows(lapply(names(ref_bins), function(r) {
    b <- ref_bins[[r]]
    dplyr::summarise(dplyr::group_by(dplyr::tibble(hour_utc = hr(b$bin_utc), v = b$value), hour_utc),
                     median_m = stats::median(v), n = dplyr::n(), .groups = "drop") |> dplyr::mutate(series = r, .before = 1)
  }))
  lay_all <- dplyr::bind_rows(lay, ref_h)
  out_tab(lay_all, "skyvue_layers_hourly")
  ids <- unique(lay_all$series)
  st <- sty(ids)
  g <- ggplot2::ggplot(lay_all, ggplot2::aes(hour_utc + 0.5, median_m, colour = series, linetype = series)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_manual(values = stats::setNames(st$col, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(st$lty, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::labs(x = "UTC hour (bin midpoint)", y = "Median height [m]",
                  subtitle = paste0(site_id, ": SkyVUE layers 1-3 (half-hourly medians), PBLHmod, ERA5; all days; internal")) +
    theme_pblh()
  p <- save_figure(g, paste0("14_skyvue_layers_hourly_", site_id, sfx), 19, 12, dir = fig_dir)
  cli::cli_alert_success("{.path {p[['png']]}}")
  cli::cli_h2("Layer medians by UTC hour (confidence 1-3: 3 = highest)")
  print_tbl(lay[lay$hour_utc %in% c(0, 6, 9, 12, 15, 18), ])
}
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
