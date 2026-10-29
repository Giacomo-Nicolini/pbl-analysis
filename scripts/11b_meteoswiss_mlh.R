# Spec 11 C: MeteoSwiss STRATfinder MLH (CHM15k ceilometer) as the standard source `meteoswiss_stratfinder_mlh` for
# the paired ICOS station (CH-Dav, Davos DAV, 0.98 km): half-hourly median of the quality-0 1-min values.
# Input: the single MeteoSwiss_STRATfinder_*.nc file under <measured_dir> (recursive). Pairing: provider `MeteoSwiss`
# in config/measured_sites.csv.
# Output: <standard_dir>/<site>/<site>__meteoswiss_stratfinder_mlh.csv.gz + .yml.
# Reports: file attributes, distance and altitude offset to the tower, quality levels, coverage per year, qc counts,
# and the night floor and midday medians by season (file facts of spec 11 C).
#
# Usage: Rscript scripts/11b_meteoswiss_mlh.R SITE
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
site_id <- commandArgs(trailingOnly = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val CH-Dav}.")
sites <- read_sites(resolve_path(cfg$sites_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$provider == "MeteoSwiss" & ms$icos_site == site_id, ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) MeteoSwiss pairing in {.file measured_sites.csv}.")
meta_site <- sites[sites$site_id == site_id, ]
meas_dir <- resolve_path(cfg$measured_dir, site = site_id)
std_dir <- file.path(resolve_path(cfg$standard_dir, must_exist = FALSE), site_id)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

f_nc <- list.files(meas_dir, pattern = "^MeteoSwiss_STRATfinder_.*\\.nc$", recursive = TRUE, full.names = TRUE)
if (length(f_nc) != 1) cli::cli_abort("{.path {meas_dir}}: expected one MeteoSwiss_STRATfinder_*.nc file, found {length(f_nc)}.")

cli::cli_h1("MeteoSwiss STRATfinder MLH for {site_id} ({ms$provider_site}, {ms$distance_km} km)")
r <- read_meteoswiss_stratfinder(f_nc)
mm <- attr(r, "meteoswiss_meta")
cli::cli_inform(c(
  "File {.file {mm$file}}: {mm$title}; product {mm$product_version}, created {mm$date_created}.",
  "Licence: {mm$licence}. Station {mm$station_name} ({mm$station_code}) {round(mm$station_lat, 5)} N, {round(mm$station_lon, 5)} E, {mm$station_alt} m; sensor {mm$sensor_height} m; daytime search ceiling {mm$search_upper_limit} m.",
  "Profiles: {nrow(r)}, {format(min(r$time_utc), '%Y-%m-%d %H:%M:%S')} to {format(max(r$time_utc), '%Y-%m-%d %H:%M:%S')} UTC."
))
d_km <- haversine_km(meta_site$lat, meta_site$lon, mm$station_lat, mm$station_lon)
if (abs(d_km - ms$distance_km) > 0.05) {
  cli::cli_abort("Distance {round(d_km, 2)} km from the file differs from {.file measured_sites.csv} ({ms$distance_km} km).")
}
cli::cli_inform("Distance to the tower {round(d_km, 3)} km; tower altitude {meta_site$altitude_m} m, ceilometer {mm$station_alt} m: offset {meta_site$altitude_m - mm$station_alt} m (heights used as is, above ground at the ceilometer; G 2026-10-30).")
dt <- diff(as.numeric(r$time_utc))
cli::cli_inform("Time step: median {stats::median(dt)} s, 1-99 % {stats::quantile(dt, 0.01)}-{stats::quantile(dt, 0.99)} s; gaps > 1 h: {sum(dt > 3600)} (longest {round(max(dt) / 86400, 1)} d).")
print_tbl(dplyr::count(dplyr::tibble(quality_mlh = r$quality_mlh), quality_mlh) |> dplyr::mutate(pct = round(100 * n / sum(n), 1)))

std <- standardize_meteoswiss_mlh(r, site_id, exclusions, meta_site$utc_offset_h)
p <- write_standard_pblh(std$data, std$meta, std_dir)
d <- std$data
cli::cli_h2("{std$meta$source_id}: {nrow(d)} half-hours, {sum(!is.na(d$pblh_m))} with a value; excluded {std$n_excluded}")
cli::cli_inform("qc 0: {sum(d$qc_flag %in% 0L)}, qc 1: {sum(d$qc_flag %in% 1L)}; quality-0 values but < 50 % of the minutes: {sum(d$n_samples > 0 & is.na(d$qc_flag))}.")
y <- format(d$timestamp_start_utc, "%Y", tz = "UTC")
print_tbl(dplyr::summarise(dplyr::group_by(dplyr::tibble(year = y, ok = !is.na(d$pblh_m)), year),
                           n_halfhours_valid = sum(ok), share = round(mean(ok), 3), .groups = "drop"))
cli::cli_alert_success("{.path {p[['data']]}}")

cli::cli_h2("Night floor and midday (half-hourly values, spec 11 C)")
h <- as.integer(format(d$timestamp_start_utc, "%H", tz = "UTC"))
se <- season_of(d$timestamp_start_utc)
q <- function(x) paste(round(stats::quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE)), collapse = " / ")
print_tbl(dplyr::bind_rows(lapply(levels(se), function(s) {
  dplyr::tibble(season = s, night_22_04_UTC_p10_p50_p90 = q(d$pblh_m[se == s & h %in% c(22, 23, 0:3)]),
                midday_11_14_UTC_p10_p50_p90 = q(d$pblh_m[se == s & h %in% 11:13]),
                share_lowest_night = round(mean(d$x_share_lowest[se == s & h %in% c(22, 23, 0:3)], na.rm = TRUE), 3))
})))
