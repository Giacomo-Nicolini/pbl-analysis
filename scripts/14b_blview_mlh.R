# Spec 14 Part B: Vaisala BL-View L3 mixing-layer height (CZ-Lnz, CL61 at the Lanzhot grassland tower) as the
# secondary standard source `vaisala_blview_mlh` (hourly, all quality values kept, quality in x_quality).
# Input: every .nc under <measured_dir> of the site (recursive; BL-View 3.1.1, 3.1.2, 3.2.2 detected from the
#   variables, never from the folder name), or under --dir. Stored time read with --time-offset (default 0), then
#   corrected (fix-14a item 11, below).
# Pairing: provider "CzechGlobe" in config/measured_sites.csv (status secondary).
# Output: <standard_dir>/<site>/<site>__vaisala_blview_mlh.csv.gz + .yml; <output_dir>/tables/internal/
#   14_blview_files_<site>.csv (one row per file: version, firmware, profiles, hours, valid hours, inconsistent
#   hours, date_stamp check, UTC offset variable, error), 14_blview_versions_<site>.csv (per version: date range,
#   days with a valid MLH, share of valid hours, quality distribution, files without MLH),
#   14_blview_overlap_<site>.csv (hours present in more than one file: version kept / dropped, how many differ).
#   Internal use only (measured data are not released).
# fix-14a item 11 (G 2026-11-12, from the 14c lag tests): time correction by default (blview_time_correction()):
#   (a) 3.2.2 stamps read as Europe/Prague local time (with summer time) and converted to UTC (nonexistent
#   spring-forward hours dropped, autumn collisions counted); (b) all versions shifted by -1 h, an empirical
#   alignment (mechanism unknown). Both, the lag evidence and the DST hours go to the sidecar; x_time_shift_h per
#   row. The version, file and overlap tables use the corrected hours (the 3.1.2 calculation-time check the stored).
#
# Usage: Rscript scripts/14b_blview_mlh.R [SITE] [--dir=PATH] [--time-offset=H] [--blview-322-tz=TZ|none]
#          [--shift-h=H]
#   Defaults: site CZ-Lnz, --blview-322-tz=Europe/Prague, --shift-h=-1 (--blview-322-tz=none --shift-h=0 keeps the
#   stored time, for checks only). --time-offset (default 0) must stay 0 with the 3.2.2 conversion.
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a) == 1) sub(paste0("^--", name, "="), "", a) else default
}
pos <- args[!startsWith(args, "--")]
bad <- setdiff(sub("=.*$", "", args[startsWith(args, "--")]), c("--dir", "--time-offset", "--blview-322-tz", "--shift-h"))
if (length(bad) > 0) cli::cli_abort("Unknown option{?s} {.val {bad}}.")
site_id <- if (length(pos) >= 1) pos[1] else "CZ-Lnz"
time_offset_h <- as.numeric(opt("time-offset", "0"))
tz_322 <- opt("blview-322-tz", "Europe/Prague")
if (identical(tz_322, "none")) tz_322 <- NULL
shift_h <- as.numeric(opt("shift-h", "-1"))
if (!is.null(tz_322) && time_offset_h != 0) cli::cli_abort("--time-offset must be 0 with the 3.2.2 local-time conversion (use --blview-322-tz=none).")
cfg <- get_cfg()
t0 <- Sys.time()
sites <- read_sites(resolve_path(cfg$sites_file))
meta_site <- sites[sites$site_id == site_id, ]
if (nrow(meta_site) != 1) cli::cli_abort("{site_id}: not in the site table.")
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
meas_excl <- read_measured_exclusions(resolve_path(cfg$measured_exclusions_file))
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
ms <- ms[ms$icos_site == site_id & ms$provider == "CzechGlobe", ]
if (nrow(ms) != 1) cli::cli_abort("{site_id}: no (or more than one) CzechGlobe BL-View row in {.file measured_sites.csv}.")
meas_dir <- opt("dir", resolve_path(cfg$measured_dir, site = site_id))
std_dir <- file.path(resolve_path(cfg$standard_dir, must_exist = FALSE), site_id)
tab_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables", "internal")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

cli::cli_h1("BL-View L3 at {site_id} ({ms$distance_km} km, status {ms$status}); stored time read with offset {time_offset_h} h")
r <- read_blview_dir(meas_dir, time_offset_h)
files <- r$files
files$path <- NULL
readr::write_csv(files, file.path(tab_dir, paste0("14_blview_files_", site_id, ".csv")))
vt <- blview_version_table(r)
readr::write_csv(vt, file.path(tab_dir, paste0("14_blview_versions_", site_id, ".csv")))
cli::cli_h2("Per version (spec 14 B3)")
print_tbl(vt)
nm <- files[files$version %in% "no_mlh", ]
if (nrow(nm) > 0) cli::cli_inform("Files without MLH variables: {nrow(nm)} (listed in {.file 14_blview_files_{site_id}.csv}, version no_mlh).")
er <- files[!is.na(files$error), ]
if (nrow(er) > 0) {
  cli::cli_alert_danger("{nrow(er)} file{?s} could not be read:")
  print_tbl(er[c("file", "error")])
}
fw <- files[!is.na(files$version) & files$version != "no_mlh", ]
cli::cli_inform("Firmware (sw_version) by version: {paste(unique(paste(fw$version, fw$sw_version)), collapse = '; ')}.")
cli::cli_inform("3.1.x date_stamp equal to time in UTC: {sum(fw$date_stamp_matches_utc %in% TRUE)} of {sum(fw$version != '3.2.2')} files; location_utc_offset values: {paste(unique(fw$location_utc_offset[fw$version == '3.1.2']), collapse = ', ')}.")
hrs <- r$hours
tz <- unique(hrs$time_zone_h[!is.na(hrs$time_zone_h)])
if (length(tz) > 0) cli::cli_inform("3.2.2 variable time_zone: {paste(tz, collapse = ', ')} h (stored time read as UTC + {time_offset_h} h).")
cli::cli_inform("Hours with a value changing inside the hour: {sum(!hrs$consistent)} (kept as NA, x_consistent FALSE).")
x312 <- hrs[hrs$version == "3.1.2" & !is.na(hrs$calc_time_utc), ]
if (nrow(x312) > 0) {
  lag <- as.numeric(x312$calc_time_utc) - as.numeric(x312$hour_end_utc)
  cli::cli_inform("3.1.2 calculation time minus hour end: median {stats::median(lag)} s, range {min(lag)} to {max(lag)} s (supports: value = the hour ending at hh:00).")
}

tc <- blview_time_correction(hrs, tz_322, shift_h)
hrs <- tc$hours
cli::cli_h2("Time correction (fix-14a item 11)")
cli::cli_inform(c(
  "3.2.2 local time: {if (is.null(tz_322)) 'not converted' else tz_322}; all versions shifted by {shift_h} h.",
  "DST: {tc$n_dropped_dst} stored 3.2.2 hour end{?s} not existing in local time dropped{if (tc$n_dropped_dst > 0) paste0(' (', paste(tc$dropped_dst, collapse = ', '), ')') else ''}; {tc$n_duplicated_dst} hour{?s} duplicated after the conversion{if (tc$n_duplicated_dst > 0) paste0(' (', paste(tc$duplicated_dst, collapse = ', '), ' UTC)') else ''}.",
  "Total shift by version (h): {paste(unique(paste(hrs$version, hrs$time_shift_h)), collapse = '; ')}."))
r$hours <- hrs
vt <- blview_version_table(r)
readr::write_csv(vt, file.path(tab_dir, paste0("14_blview_versions_", site_id, ".csv")))
dd <- blview_dedup_hours(hrs)
readr::write_csv(dd$overlap, file.path(tab_dir, paste0("14_blview_overlap_", site_id, ".csv")))
if (nrow(dd$overlap) > 0) {
  cli::cli_h2("Hours in more than one file (newest version kept)")
  print_tbl(dd$overlap)
}
std <- standardize_blview_mlh(dd$hours, site_id, exclusions, meta_site$utc_offset_h, time_offset_h, meas_excl, vt,
                              time_correction = tc[names(tc) != "hours"])
p <- write_standard_pblh(std$data, std$meta, std_dir)
d <- std$data
cli::cli_h2("{std$meta$source_id}: {nrow(d)} hours, {sum(!is.na(d$pblh_m))} with a value; excluded {std$n_excluded} (site periods), {std$n_excluded_measured} (measured-source periods)")
cli::cli_alert_success("{.path {p[['data']]}}")
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
