# Spec 07 (G / Cowork 2026-10-23): timing of the sensible heat flux, EC vs ERA5, per station.
# EC: FLUXNET H_F_MDS, each half-hour at its midpoint (start + 15 min); with the L2 FLUXES sensible-heat storage
# term SH when present (H_F_MDS + SH, on rows with SH). ERA5: tower-cell surface sensible heat flux (hourly
# accumulation, H = -sshf / 3600), at the midpoint of the accumulation hour (timestamp - 30 min).
# Per season: mean diurnal cycle, morning and afternoon crossings of H = 0 and of 20 % of the maximum, best-fit
# lag of ERA5 against each EC series (positive = ERA5 earlier). Diagnostic only: no change to the PBLHmod forcing.
# Output: <output_dir>/tables/07_h_timing_<site>.csv (times in minutes after 00:00 UTC and as HH:MM UTC).
#
# Usage: Rscript scripts/07d_h_timing.R SITE [SITE ...] [--from=2020-01-01] [--to=2025-12-31]
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
t_from <- as.POSIXct(opt("from", "2020-01-01"), tz = "UTC")
t_to <- as.POSIXct(opt("to", "2025-12-31"), tz = "UTC") + 86400
site_ids <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_ids) == 0) cli::cli_abort("Pass at least one station code.")

sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
hhmm <- function(m) ifelse(is.na(m), NA_character_, sprintf("%02d:%02d", (round(m) %/% 60) %% 24, round(m) %% 60))
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

for (site_id in site_ids) {
  cli::cli_h1("Sensible heat timing, EC vs ERA5: {site_id} ({format(t_from, '%Y-%m-%d')} to {format(t_to - 1, '%Y-%m-%d')})")
  d <- read_icos_site(site_id, cfg, sites, sites_tv, exclusions)
  cli::cli_inform("Exclusions: {sum(d$excl_flag)} half-hours dropped.")
  d <- d[!d$excl_flag & d$timestamp_start_utc >= t_from & d$timestamp_start_utc < t_to, ]
  t_ec <- d$timestamp_start_utc + 900 # half-hour midpoint
  series <- list(EC_H = dplyr::tibble(t_mid = t_ec, h = d$fn_H_F_MDS))
  has_sh <- "SH" %in% names(d) && any(!is.na(d$SH))
  if (has_sh) {
    ok <- !is.na(d$SH) & !is.na(d$fn_H_F_MDS)
    series$EC_H_rows_with_SH <- dplyr::tibble(t_mid = t_ec[ok], h = d$fn_H_F_MDS[ok])
    series$EC_H_plus_SH <- dplyr::tibble(t_mid = t_ec[ok], h = d$fn_H_F_MDS[ok] + d$SH[ok])
    cli::cli_inform("SH (L2 FLUXES storage term) valid on {sum(!is.na(d$SH))} of {nrow(d)} half-hours ({round(100 * mean(!is.na(d$SH)), 1)} %); median {round(stats::median(d$SH, na.rm = TRUE), 2)} W m-2, p5-p95 {paste(round(stats::quantile(d$SH, c(.05, .95), na.rm = TRUE), 1), collapse = ' to ')} W m-2.")
  } else {
    cli::cli_alert_info("No SH (sensible-heat storage) in L2 FLUXES: EC H only.")
  }
  sl <- find_era5_files(site_id, cfg, "single_level", sites = sites)
  hsf <- read_era5_single_level(sl, "surface_sensible_heat_flux", cells = "center")
  hsf <- hsf[hsf$timestamp_end_utc > t_from & hsf$timestamp_end_utc <= t_to, ]
  series$ERA5_H <- dplyr::tibble(t_mid = hsf$timestamp_end_utc - 1800, h = hsf$H_era5_Wm2) # accumulation midpoint

  tab <- h_timing_table(series)
  tab <- dplyr::bind_cols(dplyr::tibble(site_id = site_id), tab)
  tcols <- c("t_max_min", grep("^(morning|afternoon)_", names(tab), value = TRUE))
  for (k in tcols) tab[[paste0(sub("_min$", "", k), "_utc")]] <- hhmm(tab[[k]])
  readr::write_csv(tab, file.path(tables_dir, paste0("07_h_timing_", site_id, ".csv")))
  cli::cli_alert_success("Written {.path {file.path(tables_dir, paste0('07_h_timing_', site_id, '.csv'))}}.")
  show <- tab[, c("season", "series", "max_H", "t_max_utc", "morning_0pct_utc", "morning_20pct_utc",
                  "afternoon_20pct_utc", "afternoon_0pct_utc", "lag_era5_min")]
  show$max_H <- round(show$max_H)
  print_tbl(show)
  # EC - ERA5 differences of the crossing times [min] (positive = ERA5 earlier)
  dif <- dplyr::bind_rows(lapply(setdiff(unique(tab$series), "ERA5_H"), function(s) {
    a <- tab[tab$series == s, ]
    b <- tab[tab$series == "ERA5_H", ][match(a$season, tab$season[tab$series == "ERA5_H"]), ]
    dplyr::tibble(season = a$season, series = s,
                  morning_0 = round(a$morning_0pct - b$morning_0pct),
                  morning_20 = round(a$morning_20pct - b$morning_20pct),
                  afternoon_20 = round(a$afternoon_20pct - b$afternoon_20pct),
                  afternoon_0 = round(a$afternoon_0pct - b$afternoon_0pct), lag = a$lag_era5_min)
  }))
  cli::cli_inform("EC minus ERA5 crossing times [min] (positive = ERA5 earlier):")
  print_tbl(dif)
}
