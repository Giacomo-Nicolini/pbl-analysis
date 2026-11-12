# Spec 14 Part A: SkyVUE PRO built-in mixing-layer height (DE-Tha, TU Dresden) as the secondary standard source
# `tud_skyvue_mlh` (half-hourly, layers 1-3 kept as x_ columns, pblh_m NA until G decides a layer rule).
# Input: the daily pairs E1227_<ddmmyy>_Midnight.csv + .txt under <measured_dir> of the site (recursive; the empty
#   .nc files are ignored), or under --dir. Every day: CSV time stamps and layers paired with the TXT message 006 by
#   record order and checked (counts, layers and confidences, cloud base, profile decode vs CSV, CRC-16;
#   read_skyvue_day()). A day that fails a check is listed with the error and left out; the script then ends with
#   exit status 1 after writing everything else. fix-14c item 9: the 7-column CSV before the MLH option (until about
#   2025-06-24) is read too (cloud base 1 and vertical visibility checked, no MLH, profile check "not available"); a
#   day with only a CSV or only a TXT (e.g. 2025-06-25) is skipped with skip_reason in the day table, not an error.
#   fix-14a item 10: a day whose CSV and TXT record counts differ (June 2025, around the MLH-option installation) is
#   skipped too (skip_reason "record count mismatch: CSV n, TXT m", n_csv / n_txt in the day table; no alignment);
#   exit status 1 only for real check failures (CRC, cloud base, vertical visibility, layers, decode). The log gives
#   the mismatch days per month.
# Pairing: provider "TU Dresden" in config/measured_sites.csv (status secondary).
# Output: <standard_dir>/<site>/<site>__tud_skyvue_mlh.csv.gz + .yml; <output_dir>/tables/internal/
#   14_skyvue_days_<site>.csv (one row per day file pair: records, CRC failures, records with MLH, profile check,
#   header values, error). Internal use only (measured data are not released).
#
# Usage: Rscript scripts/14a_skyvue_mlh.R [SITE] [--from=YYYY-MM-DD] [--to=YYYY-MM-DD] [--dir=PATH]
#   Default site DE-Tha; --from / --to select days (inclusive; e.g. the July 2026 test month).
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a) == 1) sub(paste0("^--", name, "="), "", a) else default
}
pos <- args[!startsWith(args, "--")]
bad <- setdiff(sub("=.*$", "", args[startsWith(args, "--")]), c("--from", "--to", "--dir"))
if (length(bad) > 0) cli::cli_abort("Unknown option{?s} {.val {bad}}.")
site_id <- if (length(pos) >= 1) pos[1] else "DE-Tha"
cfg <- get_cfg()
t0 <- Sys.time()
sites <- read_sites(resolve_path(cfg$sites_file))
meta_site <- sites[sites$site_id == site_id, ]
if (nrow(meta_site) != 1) cli::cli_abort("{site_id}: not in the site table.")
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
meas_excl <- read_measured_exclusions(resolve_path(cfg$measured_exclusions_file))
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$icos_site == site_id & ms$provider == "TU Dresden", ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) TU Dresden SkyVUE row in {.file measured_sites.csv}.")
d_km <- haversine_km(meta_site$lat, meta_site$lon, ms$provider_lat, ms$provider_lon)
cli::cli_inform("SkyVUE PRO at {site_id}: {round(d_km * 1000)} m from the tower ({ms$distance_class}, status {ms$status}).")
meas_dir <- opt("dir", resolve_path(cfg$measured_dir, site = site_id))
std_dir <- file.path(resolve_path(cfg$standard_dir, must_exist = FALSE), site_id)
tab_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables", "internal")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

days <- find_skyvue_days(meas_dir)
if (!is.null(opt("from"))) days <- days[days$day >= as.Date(opt("from")), ]
if (!is.null(opt("to"))) days <- days[days$day <= as.Date(opt("to")), ]
if (nrow(days) == 0) cli::cli_abort("{.path {meas_dir}}: no SkyVUE day file pair in the selected period.")
cli::cli_inform("{nrow(days)} day file pair{?s}, {format(min(days$day))} to {format(max(days$day))}, under {.path {meas_dir}}.")

keep <- c("time_end_utc", "crc_ok", "detection_status", "window_tr_pct", "mlh1_m", "mlh1_q", "mlh2_m", "mlh2_q",
          "mlh3_m", "mlh3_q")
recs <- vector("list", nrow(days))
info <- vector("list", nrow(days))
for (i in seq_len(nrow(days))) {
  fname <- basename(if (is.na(days$csv[i])) days$txt[i] else days$csv[i])
  if (!is.na(days$skip_reason[i])) {
    # fix-14c item 9: a day with only one of the two files is skipped and listed (not an error)
    info[[i]] <- dplyr::tibble(day = days$day[i], csv = fname, skip_reason = days$skip_reason[i], error = NA_character_)
    next
  }
  r <- tryCatch(read_skyvue_day(days$csv[i], days$txt[i]), error = function(e) e)
  if (inherits(r, "skyvue_count_mismatch")) {
    # fix-14a item 10: different CSV / TXT record counts = skipped day (no alignment), not a failed check
    info[[i]] <- dplyr::tibble(day = days$day[i], csv = fname, n_csv = r$n_csv, n_txt = r$n_txt,
                               skip_reason = paste0("record count mismatch: CSV ", r$n_csv, ", TXT ", r$n_txt),
                               error = NA_character_)
    next
  }
  if (inherits(r, "error")) {
    msg <- gsub("\\s+", " ", conditionMessage(r))
    cli::cli_alert_danger("{format(days$day[i])}: {msg}")
    info[[i]] <- dplyr::tibble(day = days$day[i], csv = fname, skip_reason = NA_character_, error = msg)
    next
  }
  m <- attr(r, "skyvue_meta")
  info[[i]] <- dplyr::tibble(
    day = days$day[i], csv = m$csv, csv_layout = m$csv_layout, n_csv = m$n_records, n_txt = m$n_records,
    n_records = m$n_records, n_crc_fail = m$n_crc_fail,
    n_mlh1 = sum(!is.na(r$mlh1_m)), n_mlh2 = sum(!is.na(r$mlh2_m)), n_mlh3 = sum(!is.na(r$mlh3_m)),
    profile_check = m$profile_check, n_profile_records_checked = m$n_profile_records_checked,
    max_rel_diff_profile = m$max_rel_diff_profile,
    scale_pct = paste(m$scale_values, collapse = " "), n_gates = paste(unique(r$n_gates), collapse = " "),
    resolution_m = paste(unique(r$resolution_m), collapse = " "), tilt_deg = paste(unique(r$tilt_deg), collapse = " "),
    interval_s_median = m$interval_s_median, first_utc = min(r$time_end_utc), last_utc = max(r$time_end_utc),
    skip_reason = NA_character_, error = NA_character_)
  recs[[i]] <- r[keep]
  if (i %% 20 == 0) cli::cli_inform("{i} / {nrow(days)} days ({format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}).")
}
info <- dplyr::bind_rows(info)
# fixed column order (bind_rows follows the first day, which may be a skipped one)
info <- dplyr::relocate(info, dplyr::any_of(c(
  "day", "csv", "csv_layout", "n_csv", "n_txt", "n_records", "n_crc_fail", "n_mlh1", "n_mlh2", "n_mlh3", "profile_check",
  "n_profile_records_checked", "max_rel_diff_profile", "scale_pct", "n_gates", "resolution_m", "tilt_deg",
  "interval_s_median", "first_utc", "last_utc", "skip_reason", "error")))
readr::write_csv(info, file.path(tab_dir, paste0("14_skyvue_days_", site_id, ".csv")))
recs <- dplyr::bind_rows(recs)
n_fail <- sum(!is.na(info$error))
n_skip <- sum(!is.na(info$skip_reason))
if (n_skip > 0) cli::cli_alert_warning("{n_skip} day{?s} skipped: {paste(format(info$day[!is.na(info$skip_reason)]), info$skip_reason[!is.na(info$skip_reason)], collapse = '; ')}.")
is_mm <- !is.na(info$skip_reason) & startsWith(info$skip_reason, "record count mismatch")
n_mm <- sum(is_mm)
mm_month <- if (n_mm > 0) {
  tm <- table(format(info$day[is_mm], "%Y-%m"))
  paste(names(tm), tm, sep = ": ", collapse = ", ")
} else "none"
if (nrow(recs) == 0) cli::cli_abort("No day passed the checks.")
recs <- recs[order(recs$time_end_utc), ]
if (anyDuplicated(recs$time_end_utc)) {
  cli::cli_alert_warning("{sum(duplicated(recs$time_end_utc))} duplicated record time{?s} across day files: first kept.")
  recs <- recs[!duplicated(recs$time_end_utc), ]
}

cli::cli_h2("Checks over {nrow(info)} day{?s}")
ok <- info[is.na(info$error) & is.na(info$skip_reason), ]
mr <- if (any(!is.na(ok$max_rel_diff_profile))) signif(max(ok$max_rel_diff_profile, na.rm = TRUE), 2) else "not available"
cli::cli_inform(c(
  "Records {sum(ok$n_records)}; CRC mismatches {sum(ok$n_crc_fail)}; days failing a check {n_fail}; days skipped {n_skip} (record count mismatch {n_mm}, only one file {n_skip - n_mm}).",
  "Record count mismatch days per month: {mm_month}.",
  "CSV layout: full {sum(ok$csv_layout == 'full')} days, short (7 columns, before the MLH option) {sum(ok$csv_layout == 'short')} days.",
  "Profile decode vs CSV (first 20 gates, unambiguous values; full layout only): {sum(ok$n_profile_records_checked, na.rm = TRUE)} records, max relative difference {mr}.",
  "Scale parameter values: {paste(unique(ok$scale_pct), collapse = ', ')} %; gates {paste(unique(ok$n_gates), collapse = ', ')} x {paste(unique(ok$resolution_m), collapse = ', ')} m; tilt {paste(unique(ok$tilt_deg), collapse = ', ')} deg.",
  "Record step: median of the daily medians {stats::median(ok$interval_s_median)} s.",
  "First record with an MLH layer: {if (any(!is.na(recs$mlh1_m))) format(min(recs$time_end_utc[!is.na(recs$mlh1_m)]), '%Y-%m-%d %H:%M:%S') else 'none'} UTC (site: MLH option from 2025-06-25)."
))
cli::cli_inform("Records with layer 1 / 2 / 3: {sum(!is.na(recs$mlh1_m))} / {sum(!is.na(recs$mlh2_m))} / {sum(!is.na(recs$mlh3_m))}; confidence of layer 1: {paste(names(table(recs$mlh1_q[!is.na(recs$mlh1_m)])), table(recs$mlh1_q[!is.na(recs$mlh1_m)]), sep = ': ', collapse = ', ')}.")

h <- skyvue_halfhour(recs)
std <- standardize_skyvue_mlh(h, site_id, exclusions, meta_site$utc_offset_h, meas_excl,
                              days_meta = ok[c("day", "n_records", "n_crc_fail", "max_rel_diff_profile")])
p <- write_standard_pblh(std$data, std$meta, std_dir)
d <- std$data
cli::cli_h2("{std$meta$source_id}: {nrow(d)} half-hours; with layer 1 / 2 / 3: {sum(!is.na(d$x_layer1_m))} / {sum(!is.na(d$x_layer2_m))} / {sum(!is.na(d$x_layer3_m))}; excluded {std$n_excluded} (site periods), {std$n_excluded_measured} (measured-source periods)")
cli::cli_alert_success("{.path {p[['data']]}}")
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
if (n_fail > 0) {
  cli::cli_alert_danger("{n_fail} day{?s} failed a check and {?is/are} not in the standard file: see {.file 14_skyvue_days_{site_id}.csv}.")
  quit(save = "no", status = 1)
}
