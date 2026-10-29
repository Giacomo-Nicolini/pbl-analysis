# Spec 10 A3: SIRTA ReOBS hourly aerosol-based mixing-layer heights (Palaiseau) as two standard sources for the
# paired ICOS station: `sirta_mlh_chm15k` (Lufft CHM15k, STRATfinder; primary aerosol reference) and
# `sirta_mlh_cl31` (Vaisala CL31, CABAM v2019.3).
# Input: the single zip under <measured_dir> (recursive) that contains a `reobs_*.nc` file; the NetCDF is extracted
# to a temporary directory at read time (never permanently). Pairing: provider `SIRTA-ReOBS` in
# config/measured_sites.csv.
# Output: <standard_dir>/<site>/<site>__sirta_mlh_chm15k.csv.gz + .yml, <site>__sirta_mlh_cl31.csv.gz + .yml.
# Reports: file attributes (product version, uuid, licence, station), coverage per year, qc counts, and the timing
# check of spec 10 A2 (JJA hour of the median minimum and maximum, CHM15k vs CL31 Spearman since 2020).
#
# Usage: Rscript scripts/10a_sirta_reobs.R SITE
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
site_id <- commandArgs(trailingOnly = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val FR-Gri}.")
sites <- read_sites(resolve_path(cfg$sites_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$provider == "SIRTA-ReOBS" & ms$icos_site == site_id, ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) SIRTA-ReOBS pairing in {.file measured_sites.csv}.")
meta_site <- sites[sites$site_id == site_id, ]
meas_dir <- resolve_path(cfg$measured_dir, site = site_id)
std_dir <- file.path(resolve_path(cfg$standard_dir, must_exist = FALSE), site_id)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

zips <- list.files(meas_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
zips <- zips[vapply(zips, function(z) any(grepl("(^|/)reobs_[^/]*\\.nc$", utils::unzip(z, list = TRUE)$Name)),
                    logical(1))]
if (length(zips) != 1) cli::cli_abort("{.path {meas_dir}}: expected one zip with a reobs_*.nc file, found {length(zips)}.")

cli::cli_h1("SIRTA ReOBS MLH for {site_id} ({ms$provider_site}, {ms$distance_km} km)")
r <- read_reobs(zips)
rm <- attr(r, "reobs_meta")
cli::cli_inform(c(
  "File {.file {rm$file}}: {rm$title}; product {rm$product_version}, ReOBS {rm$reobs_version}, uuid {rm$uuid}.",
  "Licence: {rm$licence}. Station {rm$station_name} {rm$station_lat} N, {rm$station_lon} E, {rm$station_alt} m.",
  "Hours: {nrow(r)}, {format(min(r$timestamp_start_utc), '%Y-%m-%d %H:%M')} to {format(max(r$timestamp_start_utc), '%Y-%m-%d %H:%M')} UTC (start of hour)."
))
d_km <- haversine_km(meta_site$lat, meta_site$lon, rm$station_lat, rm$station_lon)
if (abs(d_km - ms$distance_km) > 0.05) {
  cli::cli_abort("Distance {round(d_km, 2)} km from the file differs from {.file measured_sites.csv} ({ms$distance_km} km).")
}

for (inst in c("chm15k", "cl31")) {
  std <- standardize_reobs_mlh(r, inst, site_id, exclusions, meta_site$utc_offset_h)
  p <- write_standard_pblh(std$data, std$meta, std_dir)
  d <- std$data
  cli::cli_h2("{std$meta$source_id}: {nrow(d)} hours, {sum(!is.na(d$pblh_m))} with a value; excluded {std$n_excluded}")
  cli::cli_inform("qc 0: {sum(d$qc_flag %in% 0L)}, qc 1: {sum(d$qc_flag %in% 1L)}; value in ReOBS but nb < 50 %: {sum(is.na(d$qc_flag) & !is.na(r[[paste0('mlh_', inst)]][match(d$timestamp_start_utc, r$timestamp_start_utc)]))}; NaN value with nb > 0: {sum(is.na(r[[paste0('mlh_', inst)]]) & !is.na(r[[paste0('mlh_', inst, '_nb')]]))}.")
  y <- format(d$timestamp_start_utc, "%Y", tz = "UTC")
  print_tbl(dplyr::summarise(dplyr::group_by(dplyr::tibble(year = y, ok = !is.na(d$pblh_m)), year),
                             n_hours_valid = sum(ok), .groups = "drop"))
  cli::cli_alert_success("{.path {p[['data']]}}")
}

cli::cli_h2("Timing check (spec 10 A2)")
h <- as.integer(format(r$timestamp_start_utc, "%H", tz = "UTC"))
mo <- as.integer(format(r$timestamp_start_utc, "%m", tz = "UTC"))
for (s in list(JJA = 6:8, DJF = c(12, 1, 2))) {
  i <- mo %in% s
  med <- vapply(0:23, function(k) c(stats::median(r$mlh_chm15k[i & h == k], na.rm = TRUE),
                                    stats::median(r$mlh_cl31[i & h == k], na.rm = TRUE)), numeric(2))
  print_tbl(data.frame(hour_utc_start = 0:23, chm15k = round(med[1, ]), cl31 = round(med[2, ])))
}
s20 <- r$timestamp_start_utc >= as.POSIXct("2020-01-01", tz = "UTC")
cli::cli_inform("CHM15k vs CL31 Spearman (2020 onward): {round(stats::cor(r$mlh_chm15k[s20], r$mlh_cl31[s20], method = 'spearman', use = 'complete.obs'), 3)}.")
