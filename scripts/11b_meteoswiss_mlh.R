# Spec 11 C / spec 13: MeteoSwiss MLH as a standard source for the paired ICOS station:
#   STRATfinder files (CHM15k, ~60 s; CH-Dav Davos DAV, DE-RuS Juelich JUE, DE-HoH Magdeburg MAG):
#     `meteoswiss_stratfinder_mlh`, half-hourly median of the quality-0 1-min values;
#   CABAM files (CL31, 15-min means; FI-Lom Kenttarova KEN, FR-Fon Melun MEL, FR-Lus Poitiers POI):
#     `meteoswiss_cabam_mlh`, half-hourly mean of the quality-0 15-min steps (qc 0 both, 1 one).
# Input: the single MeteoSwiss_STRATfinder_*.nc or MeteoSwiss_CABAM_*.nc file under <measured_dir> (recursive).
# Pairing: provider `MeteoSwiss` in config/measured_sites.csv (station code and distance checked against the file).
# Output: <standard_dir>/<site>/<site>__<source_id>.csv.gz + .yml.
# Reports: file attributes, distance and altitude offset to the tower, quality levels and qc bits, coverage per year,
# qc counts, instrument events (STRATfinder: window-transmission days of the file vs the profiles with bit 4), and the
# night and midday medians by season.
#
# Usage: Rscript scripts/11b_meteoswiss_mlh.R SITE
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
site_id <- commandArgs(trailingOnly = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val CH-Dav}.")
sites <- read_sites(resolve_path(cfg$sites_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
meas_excl <- read_measured_exclusions(resolve_path(cfg$measured_exclusions_file)) # spec 13
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$provider == "MeteoSwiss" & ms$icos_site == site_id, ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) MeteoSwiss pairing in {.file measured_sites.csv}.")
meta_site <- sites[sites$site_id == site_id, ]
meas_dir <- resolve_path(cfg$measured_dir, site = site_id)
std_dir <- file.path(resolve_path(cfg$standard_dir, must_exist = FALSE), site_id)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

f_nc <- list.files(meas_dir, pattern = "^MeteoSwiss_(STRATfinder|CABAM)_.*\\.nc$", recursive = TRUE, full.names = TRUE)
if (length(f_nc) != 1) {
  cli::cli_abort("{.path {meas_dir}}: expected one MeteoSwiss_STRATfinder_*.nc or MeteoSwiss_CABAM_*.nc file, found {length(f_nc)}.")
}

r <- read_meteoswiss_mlh(f_nc)
mm <- attr(r, "meteoswiss_meta")
cabam <- mm$method == "CABAM"
cli::cli_h1("MeteoSwiss {mm$method} MLH for {site_id} ({ms$provider_site}, {ms$distance_km} km)")
if (!identical(mm$station_code, ms$provider_site)) {
  cli::cli_abort("File station {.val {mm$station_code}} differs from {.file measured_sites.csv} ({.val {ms$provider_site}}).")
}
t_first <- if (cabam) min(r$interval_start_utc) else min(r$time_utc)
t_last <- if (cabam) max(r$interval_end_utc) else max(r$time_utc)
cli::cli_inform(c(
  "File {.file {mm$file}}: {mm$title}; product {mm$product_version}, created {mm$date_created}.",
  "Licence: {meteoswiss_licence(mm$licence)}",
  "Data owner: {if (is.na(mm$data_owner)) mm$institution else mm$data_owner}.",
  "Station {mm$station_name} ({mm$station_code}{if (!is.na(mm$wigos_station_id)) paste0(', WIGOS ', mm$wigos_station_id) else ''}) {round(mm$station_lat, 5)} N, {round(mm$station_lon, 5)} E, {mm$station_alt} m; sensor {mm$sensor_height} m{if (!cabam) paste0('; daytime search ceiling ', mm$search_upper_limit, ' m') else ''}.",
  "{if (cabam) '15-min steps' else 'Profiles'}: {nrow(r)}, {format(t_first, '%Y-%m-%d %H:%M:%S')} to {format(t_last, '%Y-%m-%d %H:%M:%S')} UTC."
))
d_km <- haversine_km(meta_site$lat, meta_site$lon, mm$station_lat, mm$station_lon)
if (abs(d_km - ms$distance_km) > 0.05) {
  cli::cli_abort("Distance {round(d_km, 2)} km from the file differs from {.file measured_sites.csv} ({ms$distance_km} km).")
}
cli::cli_inform("Distance to the tower {round(d_km, 3)} km ({distance_class(ms$distance_km)}); tower altitude {meta_site$altitude_m} m, ceilometer {mm$station_alt} m: offset {meta_site$altitude_m - mm$station_alt} m (heights used as is, above ground at the ceilometer; G 2026-10-30).")
if (cabam) {
  cli::cli_inform("Intervals: all 15 min, on the quarter hour, time at the centre (checked by read_meteoswiss_cabam()); tilt {paste(range(r$tilt_deg, na.rm = TRUE), collapse = '-')} deg.")
} else {
  dt <- diff(as.numeric(r$time_utc))
  cli::cli_inform("Time step: median {stats::median(dt)} s, 1-99 % {stats::quantile(dt, 0.01)}-{stats::quantile(dt, 0.99)} s; gaps > 1 h: {sum(dt > 3600)} (longest {round(max(dt) / 86400, 1)} d).")
}
print_tbl(dplyr::count(dplyr::tibble(quality_mlh = r$quality_mlh), quality_mlh) |> dplyr::mutate(pct = round(100 * n / sum(n), 1)))
bits <- meteoswiss_qc_bits(mm$method)
print_tbl(dplyr::tibble(bit = unname(bits), meaning = names(bits),
                        n = vapply(bits, function(b) sum(bitwAnd(r$qc_flag_mlh, b) > 0, na.rm = TRUE), 1L),
                        n_quality0 = vapply(bits, function(b) sum(bitwAnd(r$qc_flag_mlh[r$quality_mlh %in% 0L], b) > 0, na.rm = TRUE), 1L)))

# Instrument events (file attribute). STRATfinder E-PROFILE: no maintenance log; bit 4 = window transmission below
# 60 % (daily median), the "intervention" days of the file (= qc_events_stratfinder.csv): check that they agree.
ev <- meteoswiss_events(mm$instrument_events)
if (nrow(ev) > 0) print_tbl(dplyr::mutate(ev, description = substr(description, 1, 70)))
if (!cabam) {
  day <- as.Date(r$time_utc)
  b4 <- tapply(bitwAnd(r$qc_flag_mlh, 4L) > 0, day, mean)
  b4_days <- as.Date(names(b4)[b4 > 0])
  iv <- ev[ev$category %in% "intervention", ]
  ev_days <- if (nrow(iv) > 0) as.Date(unlist(lapply(seq_len(nrow(iv)), function(i) {
    seq(as.Date(iv$start_utc[i]), as.Date(iv$end_utc[i]), by = "day")
  })), origin = "1970-01-01") else as.Date(character())
  ev_days <- ev_days[ev_days %in% as.Date(names(b4))]
  cli::cli_inform("Bit 4 (maintenance / degraded instrument): {length(b4_days)} day{?s} with a flagged profile; intervention days in the file: {length(ev_days)}; flagged but not listed: {length(setdiff(b4_days, ev_days))}; listed but not flagged: {length(setdiff(ev_days, b4_days))}.")
  b8192 <- sum(bitwAnd(r$qc_flag_mlh, 8192L) > 0)
  if (b8192 > 0) {
    tt <- r$time_utc[bitwAnd(r$qc_flag_mlh, 8192L) > 0]
    cli::cli_inform("Bit 8192 (no temperature overlap model, informational, kept): {b8192} profiles, {format(min(tt), '%Y-%m-%d')} to {format(max(tt), '%Y-%m-%d')}; share per half-hour in x_share_no_temp_overlap.")
  }
}

std <- standardize_meteoswiss_mlh(r, site_id, exclusions, meta_site$utc_offset_h, meas_excl)
p <- write_standard_pblh(std$data, std$meta, std_dir)
d <- std$data
cli::cli_h2("{std$meta$source_id}: {nrow(d)} half-hours, {sum(!is.na(d$pblh_m))} with a value; excluded {std$n_excluded} (site periods), {std$n_excluded_measured} (measured-source periods)")
cli::cli_inform("qc 0: {sum(d$qc_flag %in% 0L)}, qc 1: {sum(d$qc_flag %in% 1L)}; quality-0 values but too few for a value: {sum(d$n_samples > 0 & is.na(d$qc_flag))}.")
y <- format(d$timestamp_start_utc, "%Y", tz = "UTC")
print_tbl(dplyr::summarise(dplyr::group_by(dplyr::tibble(year = y, ok = !is.na(d$pblh_m)), year),
                           n_halfhours_valid = sum(ok), share = round(mean(ok), 3), .groups = "drop"))
cli::cli_alert_success("{.path {p[['data']]}}")

cli::cli_h2("Night and midday (half-hourly values; UTC hours)")
h <- as.integer(format(d$timestamp_start_utc, "%H", tz = "UTC"))
se <- season_of(d$timestamp_start_utc)
q <- function(x) paste(round(stats::quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE)), collapse = " / ")
print_tbl(dplyr::bind_rows(lapply(levels(se), function(s) {
  dplyr::tibble(season = s, night_22_04_UTC_p10_p50_p90 = q(d$pblh_m[se == s & h %in% c(22, 23, 0:3)]),
                midday_11_14_UTC_p10_p50_p90 = q(d$pblh_m[se == s & h %in% 11:13]),
                share_lowest_night = round(mean(d$x_share_lowest[se == s & h %in% c(22, 23, 0:3)], na.rm = TRUE), 3))
})))
