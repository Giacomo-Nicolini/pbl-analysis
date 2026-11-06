# Spec 13 item 3: additional QC of a MeteoSwiss E-PROFILE MLH file (asked by MeteoSwiss, especially Kenttarova).
# Internal report, not committed: tables and figures under <output_dir>. Statistics on the native steps (good =
# quality level 0); day / night by apparent solar time (midday 10-16 h, night 22-04 h). Nothing is excluded here:
# the proposed exclusion periods are a table for G (never applied without approval).
# Input: the MeteoSwiss_STRATfinder_*.nc or MeteoSwiss_CABAM_*.nc file under <measured_dir> (as scripts/11b).
# Output (<output_dir>):
#   tables/13_qc_month_hour_<site>.csv   month x hour: steps, share good, median / IQR of good MLH, share at the
#                                        lowest detectable height (bit 8), share window-flagged
#   tables/13_qc_monthly_<site>.csv      per month: days with data, shares, midday / night medians, segment
#   tables/13_qc_days_<site>.csv         days of the record missing from the file, by month and sun class
#   tables/13_qc_sun_<site>.csv          share good and median MLH by hour since sunrise / until sunset, per season
#   tables/13_qc_breakpoints_<site>.csv  60 days before / after each instrument breakpoint
#   tables/13_qc_temp_overlap_<site>.csv STRATfinder bit 8192 vs temperature overlap model (Juelich)
#   tables/13_qc_events_<site>.csv       STRATfinder: window-transmission days of qc_events_stratfinder.csv (when
#                                        present in <measured_dir of _shared>) vs the file's bit 4
#   tables/13_qc_proposed_exclusions_<site>.csv  months with >= 30 % window-flagged steps (proposal only)
#   figures/13_qc_diurnal_<site>.png, 13_qc_share_good_<site>.png, 13_qc_lowest_<site>.png, 13_qc_monthly_<site>.png
#
# Usage: Rscript scripts/13a_meteoswiss_qc.R SITE
#   Profile: env var R_CONFIG_ACTIVE (default | cluster). Light: one file, ~15 s to read at most.

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
site_id <- commandArgs(trailingOnly = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val FI-Lom}.")
meas_dir <- resolve_path(cfg$measured_dir, site = site_id)
f_nc <- list.files(meas_dir, pattern = "^MeteoSwiss_(STRATfinder|CABAM)_.*\\.nc$", recursive = TRUE, full.names = TRUE)
if (length(f_nc) != 1) cli::cli_abort("{.path {meas_dir}}: expected one MeteoSwiss file, found {length(f_nc)}.")
out <- resolve_path(cfg$output_dir, must_exist = FALSE)
tdir <- file.path(out, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
wt <- function(x, what) readr::write_csv(x, file.path(tdir, paste0("13_qc_", what, "_", site_id, ".csv")))
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 3)

r <- read_meteoswiss_mlh(f_nc)
mm <- attr(r, "meteoswiss_meta")
st <- meteoswiss_qc_steps(r, mm$station_lon)
cli::cli_h1("QC of {mm$file} ({mm$station_code}, {mm$method}) for {site_id}")
cli::cli_inform("{nrow(st)} steps, {format(min(st$date))} to {format(max(st$date))}; good {round(100 * mean(st$good), 1)} %; window flag {round(100 * mean(st$window), 1)} %; bit 8 among good {round(100 * mean(st$lowest[st$good]), 1)} %.")

# 1-3: monthly diurnal cycle, share good by hour, lowest detectable height
mh <- qc_month_hour(st)
wt(mh, "month_hour")
cli::cli_h2("Share good and median good MLH by hour (all months)")
hh <- dplyr::summarise(dplyr::group_by(st, hour_utc), share_good = mean(good), median_good_m = qc_med(mlh_m[good]),
                       share_lowest_good = mean(lowest[good]), .groups = "drop")
print_tbl(hh)
low <- st$mlh_m[st$good]
cli::cli_inform("Good values <= 50 m: {sum(low <= 50)} ({round(100 * mean(low <= 50), 2)} %); <= 100 m: {round(100 * mean(low <= 100), 2)} %; minimum {min(low)} m.")

# 4: sunrise / sunset, polar day and night, days missing from the file
dy <- qc_days(st, mm$station_lat, mm$station_lon)
dm <- dplyr::summarise(dplyr::group_by(dplyr::mutate(dy, year_month = format(date, "%Y-%m")), year_month),
                       days = dplyr::n(), missing = sum(!in_file), polar_day = sum(sun == "polar_day"),
                       polar_night = sum(sun == "polar_night"), missing_polar_day = sum(!in_file & sun == "polar_day"),
                       missing_polar_night = sum(!in_file & sun == "polar_night"), .groups = "drop")
wt(dm, "days")
cli::cli_h2("Days of the record missing from the file: {sum(!dy$in_file)} of {nrow(dy)}")
print_tbl(dplyr::count(dy[!dy$in_file, ], sun, name = "missing_days"))
print_tbl(dm[dm$missing > 0, c("year_month", "days", "missing", "missing_polar_day", "missing_polar_night")])
if (any(dy$sun == "polar_day")) {
  pd <- dy[format(dy$date, "%m") %in% c("05", "06", "07", "08"), ]
  cli::cli_inform("May-August: {sum(!pd$in_file)} missing days of {nrow(pd)} ({sum(!pd$in_file & pd$sun == 'polar_day')} in polar day) (MeteoSwiss: 49 days stopped in the CABAM noon search).")
}
sr <- qc_sun_relative(st, mm$station_lat, mm$station_lon)
wt(sr, "sun")
cli::cli_h2("Around sunrise and sunset (normal days; all seasons pooled)")
print_tbl(dplyr::summarise(dplyr::group_by(sr, event, h_rel), share_good = sum(share_good * n) / sum(n), n = sum(n),
                           .groups = "drop"))

# 5: breakpoints
bp <- qc_breakpoints(st)
wt(bp, "breakpoints")
if (nrow(bp) > 0) {
  cli::cli_h2("60 days before / after each breakpoint")
  print_tbl(bp)
}

# 6: window flag by month (dirty window at Melun), proposals
mon <- qc_monthly(st)
wt(mon, "monthly")
cli::cli_h2("Months with >= 10 % window-flagged steps")
print_tbl(mon[mon$share_window >= 0.1, c("year_month", "days_with_data", "share_good", "share_window", "median_midday_m", "median_night_m")])
cmp <- dplyr::summarise(dplyr::group_by(dplyr::mutate(mon, flagged = share_window >= 0.3), flagged),
                        months = dplyr::n(), median_midday_m = stats::median(median_midday_m, na.rm = TRUE),
                        median_night_m = stats::median(median_night_m, na.rm = TRUE), .groups = "drop")
print_tbl(cmp)
px <- qc_proposed_exclusions(mon)
wt(px, "proposed_exclusions")
cli::cli_h2("Proposed exclusion periods (NOT applied; for G): {nrow(px)} month{?s}")
if (nrow(px) > 0) print_tbl(px)

# 7: Juelich without temperature overlap model (bit 8192)
to <- qc_temp_overlap(st)
if (nrow(to) > 0) {
  wt(to, "temp_overlap")
  cli::cli_h2("Bit 8192 (no temperature overlap model) vs temperature model, same calendar months")
  print_tbl(to)
}

# STRATfinder: window-transmission days of the shared qc_events file vs bit 4 in the file
if (mm$method == "STRATfinder") {
  qe <- file.path(resolve_path(cfg$measured_dir, site = "_shared", must_exist = FALSE), "meteoswiss_eprofile_2026-11-05",
                  "qc_events_stratfinder.csv")
  if (file.exists(qe)) {
    ev <- readr::read_csv(qe, col_types = readr::cols(.default = "c"), progress = FALSE)
    ev <- ev[ev$stn == mm$station_code & ev$category == "intervention", ]
    days <- if (nrow(ev)) as.Date(unlist(lapply(seq_len(nrow(ev)), function(i) {
      seq(as.Date(substr(ev$start_utc[i], 1, 10)), as.Date(substr(ev$end_utc[i], 1, 10)), by = "day")
    })), origin = "1970-01-01") else as.Date(character())
    b4 <- tapply(st$window, st$date, mean)
    chk <- dplyr::tibble(date = days, in_file = days %in% st$date, share_bit4 = as.numeric(b4[as.character(days)]))
    wt(chk, "events")
    b4_days <- as.Date(names(b4)[b4 > 0])
    cli::cli_inform("qc_events_stratfinder.csv: {length(days)} intervention days; with bit 4 on >= 99.9 % of their profiles: {sum(chk$share_bit4 >= 0.999, na.rm = TRUE)} (the events end at 23:59:00, so the last profile of the day can be unflagged); days with bit 4 not listed: {length(setdiff(b4_days, days))}.")
  } else {
    cli::cli_alert_warning("{.path {qe}} not found: events check skipped.")
  }
}

# Figures
br <- if (nrow(bp)) unique(bp$breakpoint) else as.Date(character())
ttl <- paste0(site_id, ": ", mm$station_code, " ", mm$method, " (quality 0)")
figs <- list(
  list(plot_qc_month_hour(mh, "median_good_m", "Median MLH [m]", ttl), "13_qc_diurnal", 16, 10),
  list(plot_qc_month_hour(mh, "share_good", "Share good [-]", ttl), "13_qc_share_good", 16, 10),
  list(plot_qc_month_hour(mh, "share_lowest_good", "Share bit 8 [-]", ttl), "13_qc_lowest", 16, 10),
  list(plot_qc_monthly(mon, br, ttl), "13_qc_monthly", 20, 12)
)
for (fg in figs) {
  pth <- save_figure(fg[[1]], paste0(fg[[2]], "_", site_id), fg[[3]], fg[[4]])
  cli::cli_alert_success("{.path {pth[['png']]}}")
}
