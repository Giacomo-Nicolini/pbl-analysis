# Spec 07 section 2: download the daily Cloudnet files of one product for the Cloudnet site paired with an ICOS
# station (config/measured_sites.csv), over the station's reference period (FLUXES range, CLAUDE.md).
# Only days with files in the API are requested, i.e. the months with coverage in the inventory (script 07a).
# Without --confirm: size report only (number of files, total size from the API metadata), nothing downloaded.
# With --confirm: downloads the files not yet present with the API SHA-256 (failures logged, retried once) and
# writes the provenance table <cloudnet_dir>/cloudnet_files_provenance_<product>.csv (uuid, checksum, version,
# PID, instrument PID ...). Aborts if the total size exceeds 10 GB (spec 07: stop and ask G).
# Output: <cloudnet_dir>/<filename>.nc (gitignored, under data/).
#
# Usage: Rscript scripts/07b_cloudnet_download.R SITE [--confirm] [--product=epsilon-lidar]
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

max_total_gb <- 10 # spec 07 section 2

cfg <- get_cfg()
args <- commandArgs(trailingOnly = TRUE)
confirm <- "--confirm" %in% args
product <- sub("^--product=", "", grep("^--product=", args, value = TRUE))
if (length(product) == 0) product <- "epsilon-lidar"
site_id <- setdiff(grep("^--", args, value = TRUE, invert = TRUE), "")
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val FI-Hyy}.")

sites <- read_sites(resolve_path(cfg$sites_file))
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$provider == "cloudnet" & ms$icos_site == site_id, ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) Cloudnet pairing in {.file measured_sites.csv}.")
meta <- sites[sites$site_id == site_id, ]
dest_dir <- resolve_path(cfg$cloudnet_dir, site = site_id, must_exist = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

ref <- icos_reference_period(site_id, cfg, meta$utc_offset_h)
d_from <- as.Date(ref[1], tz = "UTC")
d_to <- as.Date(ref[2], tz = "UTC")
cli::cli_h1("Cloudnet {product}: {ms$provider_site} for {site_id}")
cli::cli_inform("Reference period (FLUXES range): {format(ref[1], '%Y-%m-%d %H:%M')} to {format(ref[2], '%Y-%m-%d %H:%M')} UTC; files requested for {d_from} to {d_to}.")

files <- cloudnet_files(ms$provider_site, product, d_from, d_to, retries = 2)
if (nrow(files) == 0) cli::cli_abort("No {product} files for {ms$provider_site} in {d_from} to {d_to}.")
plan <- cloudnet_download_plan(files, dest_dir)

# --- Size report (API metadata) ---
total_gb <- sum(plan$size_bytes) / 1e9
todo_gb <- sum(plan$size_bytes[!plan$present]) / 1e9
cli::cli_h2("Size report (API metadata)")
cli::cli_inform(c(
  "{nrow(plan)} daily file{?s}, {format(round(total_gb * 1000, 1), nsmall = 1)} MB in total ({round(total_gb, 3)} GB); first {min(plan$measurement_date)}, last {max(plan$measurement_date)}.",
  "Already present with the API checksum: {sum(plan$present)}; to download: {sum(!plan$present)} ({format(round(todo_gb * 1000, 1), nsmall = 1)} MB).",
  "Volatile (no PID yet, may change): {sum(plan$volatile %in% TRUE)}; legacy: {sum(plan$legacy %in% TRUE)}."
))
per_year <- dplyr::summarise(dplyr::group_by(plan, year = format(.data$measurement_date, "%Y")),
                             n_files = dplyr::n(), size_MB = sum(.data$size_bytes) / 1e6,
                             n_months = dplyr::n_distinct(format(.data$measurement_date, "%Y-%m")), .groups = "drop")
print_tbl(per_year)
ins <- dplyr::count(plan, .data$instrument_name, .data$instrument_pid, name = "n_files")
print_tbl(ins)
n_day_multi <- sum(table(plan$measurement_date) > 1)
if (n_day_multi > 0) cli::cli_alert_warning("{n_day_multi} day{?s} with more than one file (several instruments).")

if (total_gb > max_total_gb) {
  cli::cli_abort("Total size {round(total_gb, 2)} GB exceeds {max_total_gb} GB: stop and ask G (spec 07 section 2).")
}
if (!confirm) {
  cli::cli_alert_info("Size report only. Rerun with {.code --confirm} to download into {.path {dest_dir}}.")
  quit(save = "no", status = 0)
}

# --- Download ---
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
res <- cloudnet_download(plan, retries = 1)
prov <- res[, c("uuid", "filename", "site", "product", "measurement_date", "instrument_name", "instrument_pid",
                "size_bytes", "checksum_sha256", "version", "pid", "volatile", "legacy", "coverage",
                "error_level", "download_url", "status", "message", "downloaded_utc")]
prov_file <- file.path(dest_dir, paste0("cloudnet_files_provenance_", product, ".csv"))
if (file.exists(prov_file)) { # keep the download time of files already present
  old <- readr::read_csv(prov_file, col_types = readr::cols(.default = "c"), progress = FALSE)
  keep <- match(prov$uuid, old$uuid)
  was <- !is.na(keep) & prov$status == "present"
  prov$downloaded_utc[was] <- as.POSIXct(old$downloaded_utc[keep[was]], tz = "UTC")
}
readr::write_csv(prov, prov_file)
cli::cli_h2("Download summary")
print_tbl(dplyr::count(res, .data$status))
if (any(res$status == "failed")) {
  print_tbl(res[res$status == "failed", c("filename", "message")])
  cli::cli_alert_danger("{sum(res$status == 'failed')} file{?s} failed (see {.path {prov_file}}).")
} else {
  cli::cli_alert_success("All {nrow(res)} files present with the API checksum. Provenance: {.path {prov_file}}.")
}
