# Spec 07 section 1: inventory of open ACTRIS/Cloudnet data at the Cloudnet sites paired with ICOS stations
# (config/measured_sites.csv, provider "cloudnet"). No data download.
# For each site and product (doppler-lidar, epsilon-lidar, doppler-lidar-wind, lidar, mwr, classification) the
# Cloudnet API /api/files is queried one calendar year at a time, one request at a time with a short pause
# (legacy files included, flagged). Also checks distance_km in measured_sites.csv against the API site coordinates.
# Output: <output_dir>/tables/07_cloudnet_inventory.csv (site x product x instrument x month) and
#         <output_dir>/tables/07_cloudnet_files.csv (one row per file, API metadata).
#
# Usage: Rscript scripts/07a_cloudnet_inventory.R [FROM_DATE [TO_DATE]]   (default 2020-01-01 to today)
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
args <- commandArgs(trailingOnly = TRUE)
date_from <- as.Date(if (length(args) >= 1) args[1] else "2020-01-01")
date_to <- as.Date(if (length(args) >= 2) args[2] else Sys.Date())
products <- c("doppler-lidar", "epsilon-lidar", "doppler-lidar-wind", "lidar", "mwr", "classification")

sites <- read_sites(resolve_path(cfg$sites_file))
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$provider == "cloudnet", ]
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

cli::cli_h1("Cloudnet inventory, {date_from} to {date_to}")

# --- Site coordinates and distance check ---
api_sites <- cloudnet_api_get("sites")
chk <- merge(ms[, c("icos_site", "provider_site", "provider_lat", "provider_lon", "distance_km")],
             api_sites[, c("id", "latitude", "longitude", "altitude")], by.x = "provider_site", by.y = "id",
             all.x = TRUE)
chk <- merge(chk, sites[, c("site_id", "lat", "lon")], by.x = "icos_site", by.y = "site_id", all.x = TRUE)
chk$distance_api_km <- round(haversine_km(chk$lat, chk$lon, chk$latitude, chk$longitude), 2)
chk$distance_class <- distance_class(chk$distance_km)
cli::cli_h2("Pairing and distance rule")
print_tbl(chk[, c("icos_site", "provider_site", "altitude", "distance_km", "distance_api_km", "distance_class")])
bad <- is.na(chk$distance_api_km) | abs(chk$distance_api_km - chk$distance_km) > 0.05
if (any(bad)) {
  cli::cli_alert_danger("distance_km in measured_sites.csv differs from the API coordinates for {chk$icos_site[bad]}.")
} else {
  cli::cli_alert_success("distance_km consistent with the API site coordinates (within 0.05 km).")
}

# --- File lists ---
files <- list()
for (s in ms$provider_site) {
  for (p in products) {
    x <- cloudnet_files(s, p, date_from, date_to, retries = 2)
    cli::cli_inform("{s} / {p}: {nrow(x)} file{?s}.")
    files[[paste(s, p)]] <- x
  }
}
files <- dplyr::bind_rows(files)
inv <- cloudnet_inventory_summary(files)
readr::write_csv(inv, file.path(tables_dir, "07_cloudnet_inventory.csv"))
readr::write_csv(files, file.path(tables_dir, "07_cloudnet_files.csv"))
cli::cli_alert_success("Written {.path {file.path(tables_dir, '07_cloudnet_inventory.csv')}} ({nrow(inv)} rows) and {.path {file.path(tables_dir, '07_cloudnet_files.csv')}} ({nrow(files)} files).")

# --- Report ---
cli::cli_h2("Totals per site and product")
tot <- dplyr::summarise(dplyr::group_by(files, .data$site, .data$product),
                        n_files = dplyr::n(), size_GB = sum(.data$size_bytes) / 1e9,
                        first = min(.data$measurement_date), last = max(.data$measurement_date),
                        n_legacy = sum(.data$legacy %in% TRUE), n_volatile = sum(.data$volatile %in% TRUE),
                        .groups = "drop")
print_tbl(tot)

cli::cli_h2("Instruments per site and product")
ins <- dplyr::summarise(dplyr::group_by(files, .data$site, .data$product, .data$instrument_name, .data$instrument_pid),
                        n_files = dplyr::n(), first = min(.data$measurement_date),
                        last = max(.data$measurement_date), .groups = "drop")
print_tbl(ins)

for (p in c("epsilon-lidar", "doppler-lidar")) {
  cli::cli_h2("Coverage periods (consecutive months with files): {p}")
  per <- cloudnet_coverage_periods(inv[inv$product == p, ])
  if (nrow(per) == 0) cli::cli_inform("No files.") else print_tbl(per)
}
