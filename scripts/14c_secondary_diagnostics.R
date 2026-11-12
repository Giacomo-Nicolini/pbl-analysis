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
#   3. CZ-Lnz quality: window metrics of PBLHmod and ERA5 against the BL-View value by quality class (1, 2, 3, all)
#      and version.
# Inputs: <standard_dir>/<site>/<site>__{tud_skyvue_mlh|vaisala_blview_mlh, pblhmod_default, era5_blh_center}
#   (scripts 14a / 14b, 04), <derived_dir>/pblhmod/<site>_pblhmod_default.csv.gz (day_id), the ICOS L2 (SW_IN,
#   SW_IN_POT) for the clear days.
# Output (internal): <output_dir>/tables/internal/14_lag_<site>.csv (series x reference x season x lag),
#   14_lag_best_<site>.csv, 14_diurnal_clear_<site>.csv, 14_window_metrics_<site>.csv, DE-Tha also
#   14_skyvue_layers_hourly_<site>.csv; <output_dir>/figures/internal/14_lag_<site>.png, 14_diurnal_clear_<site>.png,
#   DE-Tha 14_skyvue_layers_hourly_<site>.png.
#
# Usage: Rscript scripts/14c_secondary_diagnostics.R SITE [--from=YYYY-MM-DD] [--to=YYYY-MM-DD]
#   e.g. DE-Tha --from=2026-07-01 --to=2026-07-31 (test month), then the full record without --from / --to.
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a) == 1) sub(paste0("^--", name, "="), "", a) else default
}
pos <- args[!startsWith(args, "--")]
bad <- setdiff(sub("=.*$", "", args[startsWith(args, "--")]), c("--from", "--to"))
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
ms <- ms[ms$icos_site == site_id & ms$status %in% "secondary", ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) secondary source in {.file measured_sites.csv}.")
src <- c(`TU Dresden` = "tud_skyvue_mlh", CzechGlobe = "vaisala_blview_mlh")[[ms$provider]]
std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
rd <- function(s) read_standard_pblh(file.path(std_dir, paste0(site_id, "__", s, ".csv.gz")))
tab_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables", "internal")
fig_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "figures", "internal")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
out_tab <- function(x, name) readr::write_csv(x, file.path(tab_dir, paste0("14_", name, "_", site_id, ".csv")))
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
cli::cli_h1("{site_id} {src}: {format(t_from, '%Y-%m-%d')} to {format(t_to - 1, '%Y-%m-%d')} ({nrow(meas)} rows)")
cli::cli_inform("PBLHmod rows {nrow(mod)} (valid {sum(!is.na(mod$pblh_m))}); ERA5 hours {nrow(era$era5_blh_center)}.")

# series of the measured source: SkyVUE layers (half-hourly), BL-View all versions and each version (hourly)
series <- if (src == "tud_skyvue_mlh") {
  stats::setNames(lapply(1:3, function(k) dplyr::mutate(meas, pblh_m = .data[[paste0("x_layer", k, "_m")]])),
                  paste0("skyvue_layer", 1:3))
} else {
  v <- sort(unique(meas$x_version))
  c(list(blview_all = meas), stats::setNames(lapply(v, function(x) meas[meas$x_version %in% x, ]), paste0("blview_", v)))
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
lag <- dplyr::bind_rows(lapply(names(series), function(s) dplyr::bind_rows(lapply(names(ref_bins), function(r) {
  dplyr::mutate(lag_correlation(meas_bins[[s]], ref_bins[[r]], mb), site_id = site_id, series = s, reference = r,
                .before = 1)
}))))
out_tab(lag, "lag")
best <- best_lag(lag, by = c("site_id", "series", "reference"))
out_tab(best, "lag_best")
cli::cli_h2("Best lag of the morning growth (L > 0: measured stamps later than UTC, e.g. +1 = CET)")
if (nrow(best) > 0) print_tbl(best) else cli::cli_alert_warning("No correlation (too few clear-day morning pairs).")

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
  p <- save_figure(g, paste0("14_diurnal_clear_", site_id), 19, 14, dir = fig_dir)
  cli::cli_alert_success("{.path {p[['png']]}}")
}
lagp <- lag[!is.na(lag$r), ]
if (nrow(lagp) > 0) {
  ids <- unique(lagp$series)
  st <- sty(ids)
  g <- ggplot2::ggplot(lagp, ggplot2::aes(lag_h, r, colour = series, linetype = series)) +
    ggplot2::geom_vline(xintercept = 0, colour = style_of("reference_line", "colour", "element"), linewidth = 0.3) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::facet_grid(reference ~ season) +
    ggplot2::scale_colour_manual(values = stats::setNames(st$col, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(st$lty, ids), labels = stats::setNames(st$lab, ids), name = NULL) +
    ggplot2::scale_x_continuous(breaks = -3:3) +
    ggplot2::labs(x = "Lag L [h] (measured value L hours after the reference hour)", y = "Pearson r, morning growth",
                  subtitle = paste0(site_id, ", clear-day mornings (sunrise to solar noon); L = 0: UTC, +1: CET, +2: CEST")) +
    theme_pblh()
  p <- save_figure(g, paste0("14_lag_", site_id), 19, 12, dir = fig_dir)
  cli::cli_alert_success("{.path {p[['png']]}}")
}

# --- window metrics (aerosol window, PBLHmod convective) per series and, for BL-View, per quality class ---
exact_ref <- function(v) v[v$kind %in% "exact", c("timestamp_utc", "pblh_m")]
win_metrics <- function(s, label) {
  s <- s[!is.na(s$timestamp_start_utc), ]
  if (sum(!is.na(s$pblh_m)) == 0) return(NULL)
  tb <- measured_time_base(s, mod, era, step_s = stats::median(as.numeric(s$timestamp_end_utc - s$timestamp_start_utc,
                                                                          units = "secs")))
  prs <- list(pblhmod_default = list(x = tb$mod, ref = exact_ref(tb$meas), mid = tb$mid_s),
              era5_blh_center = list(x = tb$era$era5_blh_center, ref = exact_ref(tb$era_meas), mid = tb$era_mid_s))
  dplyr::bind_rows(lapply(names(prs), function(r) {
    pt <- pair_table(prs[[r]]$x, prs[[r]]$ref, r)
    a <- apparent_solar_hour(pt$timestamp_utc + prs[[r]]$mid, meta$lon)
    pt <- pt[pt$regime %in% "convective" & a >= opts$aerosol_window[1] & a < opts$aerosol_window[2], ]
    q <- if ("x_quality" %in% names(s)) s$x_quality[match(as.numeric(pt$timestamp_utc), as.numeric(s$timestamp_start_utc))]
    cls <- if (is.null(q)) list(all = rep(TRUE, nrow(pt))) else
      c(list(all = rep(TRUE, nrow(pt))), stats::setNames(lapply(1:3, function(k) q %in% k), paste0("q", 1:3)))
    dplyr::bind_rows(lapply(names(cls), function(k) {
      x <- pt[cls[[k]], ]
      dplyr::bind_cols(dplyr::tibble(site_id = site_id, series = label, source = r, quality = k,
                                     n_days = length(unique(as.Date(x$timestamp_utc)))),
                       pblh_metrics(x$mod, x$ref))
    }))
  }))
}
wm <- dplyr::bind_rows(lapply(names(series), function(s) win_metrics(series[[s]], s)))
wm$note <- paste("ln(source / measured); window", paste(opts$aerosol_window, collapse = "-"),
                 "h apparent solar time, PBLHmod convective; quality all / q1-q3 (BL-View only)")
out_tab(wm, "window_metrics")
cli::cli_h2("Window metrics (source vs measured)")
if (nrow(wm) > 0) print_tbl(wm[c("series", "source", "quality", "n", "n_days", "median_log_ratio", "frac_within_factor2")])

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
  p <- save_figure(g, paste0("14_skyvue_layers_hourly_", site_id), 19, 12, dir = fig_dir)
  cli::cli_alert_success("{.path {p[['png']]}}")
  cli::cli_h2("Layer medians by UTC hour (confidence 1-3: 3 = highest)")
  print_tbl(lay[lay$hour_utc %in% c(0, 6, 9, 12, 15, 18), ])
}
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
