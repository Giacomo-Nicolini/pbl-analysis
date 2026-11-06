# Additional QC of the MeteoSwiss E-PROFILE MLH files (spec 13 item 3; asked by MeteoSwiss, especially Kenttarova).
# Internal report only (scripts/13a_meteoswiss_qc.R): tables and figures under <output_dir>, never committed.
# All statistics on the native steps (STRATfinder ~60 s profiles, CABAM 15-min steps); "good" = quality level 0
# (provider recommendation). Day / night by apparent solar time (AST): midday 10-16 h (the aerosol window of spec 10,
# measured_compare_options()), night 22-04 h.

#' Native MeteoSwiss steps in one table
#'
#' @param r Output of [read_meteoswiss_stratfinder()] or [read_meteoswiss_cabam()].
#' @param lon Ceilometer longitude [deg E] (apparent solar time).
#' @return Tibble: `time_utc` [POSIXct UTC; profile time, or the start of the CABAM step], `mlh_m` [m agl], `good`
#'   [logical, quality 0 with a value], `quality`, `qc_flag_mlh`, `segment` [integer], `window` [logical: CABAM bit 2
#'   dirty window, STRATfinder bit 4 window transmission < 60 % (no maintenance log)], `lowest` [logical, bit 8],
#'   `no_temp_overlap` [logical, STRATfinder bit 8192], `ast_h` [h, apparent solar time at the step midpoint],
#'   `month`, `hour_utc` [integer], `date` [Date, UTC]. Attribute `step_s` [s, nominal step].
meteoswiss_qc_steps <- function(r, lon) {
  cabam <- identical(attr(r, "meteoswiss_meta")$method, "CABAM")
  t <- if (cabam) r$interval_start_utc else r$time_utc
  step <- if (cabam) 900 else 60
  q <- r$qc_flag_mlh
  out <- dplyr::tibble(time_utc = t, mlh_m = r$mlh_m, good = r$quality_mlh %in% 0L & !is.na(r$mlh_m),
                       quality = r$quality_mlh, qc_flag_mlh = q, segment = r$segment,
                       window = bitwAnd(q, if (cabam) 2L else 4L) > 0, lowest = bitwAnd(q, 8L) > 0,
                       no_temp_overlap = if (cabam) FALSE else bitwAnd(q, 8192L) > 0,
                       ast_h = apparent_solar_hour(t + step / 2, lon),
                       month = as.integer(format(t, "%m", tz = "UTC")), hour_utc = as.integer(format(t, "%H", tz = "UTC")),
                       date = as.Date(t, tz = "UTC"))
  attr(out, "step_s") <- step
  out
}

# AST classes used throughout: midday 10-16 h, night 22-04 h.
qc_midday <- function(ast) ast >= 10 & ast < 16
qc_night <- function(ast) ast >= 22 | ast < 4
qc_med <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)

#' Month x hour summary of the MeteoSwiss steps
#'
#' @param st Output of [meteoswiss_qc_steps()].
#' @return Tibble per `month` x `hour_utc`: `n` (steps in the file), `share_good` [0-1], `median_good_m`,
#'   `p25_good_m`, `p75_good_m` [m agl], `share_lowest_good` [0-1, bit 8 among the good steps], `share_window` [0-1].
qc_month_hour <- function(st) {
  dplyr::summarise(dplyr::group_by(st, month, hour_utc), n = dplyr::n(), share_good = mean(good),
                   median_good_m = qc_med(mlh_m[good]),
                   p25_good_m = if (any(good)) unname(stats::quantile(mlh_m[good], 0.25)) else NA_real_,
                   p75_good_m = if (any(good)) unname(stats::quantile(mlh_m[good], 0.75)) else NA_real_,
                   share_lowest_good = if (any(good)) mean(lowest[good]) else NA_real_,
                   share_window = mean(window), .groups = "drop")
}

#' Monthly series of the MeteoSwiss steps
#'
#' @param st Output of [meteoswiss_qc_steps()].
#' @return Tibble per `year_month` (`YYYY-MM`): `days_with_data`, `days_in_month`, `n`, `share_good`, `share_window`
#'   [0-1], `share_no_temp_overlap` [0-1], `median_midday_m`, `median_night_m` [m agl, good steps, AST classes],
#'   `share_lowest_night` [0-1, bit 8 among the good night steps], `segment` [largest index].
qc_monthly <- function(st) {
  st$year_month <- format(st$time_utc, "%Y-%m", tz = "UTC")
  out <- dplyr::summarise(dplyr::group_by(st, year_month), days_with_data = dplyr::n_distinct(date), n = dplyr::n(),
                          share_good = mean(good), share_window = mean(window), share_no_temp_overlap = mean(no_temp_overlap),
                          median_midday_m = qc_med(mlh_m[good & qc_midday(ast_h)]),
                          median_night_m = qc_med(mlh_m[good & qc_night(ast_h)]),
                          share_lowest_night = if (any(good & qc_night(ast_h))) mean(lowest[good & qc_night(ast_h)]) else NA_real_,
                          segment = max(segment), .groups = "drop")
  d1 <- as.Date(paste0(out$year_month, "-01"))
  out$days_in_month <- as.integer(as.Date(format(d1 + 31, "%Y-%m-01")) - d1)
  out[c("year_month", "days_with_data", "days_in_month", setdiff(names(out), c("year_month", "days_with_data", "days_in_month")))]
}

#' Days without data and the sun class of every day
#'
#' Days between the first and the last day of the file that are absent from it (E-PROFILE files: days without any
#' data are not in the file), with the sun class from [sun_times()] (polar day / polar night / normal).
#'
#' @param st Output of [meteoswiss_qc_steps()].
#' @param lat,lon Ceilometer position [deg].
#' @return Tibble per day of the record: `date`, `in_file` [logical], `sun` (`"polar_day"`, `"polar_night"`,
#'   `"normal"`), `n_good` [steps].
qc_days <- function(st, lat, lon) {
  days <- seq(min(st$date), max(st$date), by = "day")
  s <- sun_times(days, lat, lon)
  ng <- tapply(st$good, st$date, sum)
  dplyr::tibble(date = days, in_file = days %in% st$date,
                sun = ifelse(s$polar %in% "day", "polar_day", ifelse(s$polar %in% "night", "polar_night", "normal")),
                n_good = as.integer(ifelse(is.na(ng[as.character(days)]), 0L, ng[as.character(days)])))
}

#' MLH around sunrise and sunset
#'
#' Steps of the normal days (not polar) binned by whole hours since sunrise and until sunset (geometric, [sun_times()]):
#' share good and median good MLH per season. The transitions are where aerosol retrievals are least reliable (the
#' residual layer is still the strongest gradient after sunrise; spec 10, Kotthaus2020).
#'
#' @param st Output of [meteoswiss_qc_steps()].
#' @param lat,lon Ceilometer position [deg].
#' @param before,after Hours before and after the event [h].
#' @return Tibble: `event` (`sunrise`, `sunset`), `season`, `h_rel` [h, start of the bin relative to the event],
#'   `n`, `share_good`, `median_good_m`.
qc_sun_relative <- function(st, lat, lon, before = 3, after = 6) {
  s <- sun_times(sort(unique(st$date)), lat, lon)
  i <- match(st$date, s$date)
  mid <- st$time_utc + attr(st, "step_s") / 2
  res <- list()
  for (ev in c("sunrise", "sunset")) {
    e <- s[[paste0(ev, "_utc")]][i]
    h <- floor(as.numeric(difftime(mid, e, units = "hours")))
    keep <- !is.na(h) & h >= -(if (ev == "sunrise") before else after) & h < (if (ev == "sunrise") after else before)
    d <- dplyr::tibble(season = season_of(st$time_utc[keep]), h_rel = h[keep], good = st$good[keep], mlh = st$mlh_m[keep])
    res[[ev]] <- dplyr::mutate(dplyr::summarise(dplyr::group_by(d, season, h_rel), n = dplyr::n(), share_good = mean(good),
                                                median_good_m = qc_med(mlh[good]), .groups = "drop"),
                               event = ev, .before = 1)
  }
  dplyr::bind_rows(res)
}

#' MLH before and after the instrument breakpoints
#'
#' For every breakpoint (start of a new `segment`), the good steps in the `window_d` days before and after:
#' midday and night medians, share of night values at the lowest detectable height, share good.
#'
#' @param st Output of [meteoswiss_qc_steps()].
#' @param window_d Days on each side [d].
#' @return Tibble per breakpoint x side (`before`, `after`): `breakpoint` [Date], `segment_after`, `side`, `n`,
#'   `share_good`, `median_midday_m`, `median_night_m`, `share_lowest_night`.
qc_breakpoints <- function(st, window_d = 60) {
  seg_start <- tapply(st$time_utc, st$segment, min)
  seg_start <- seg_start[order(as.integer(names(seg_start)))][-1]
  if (length(seg_start) == 0) return(dplyr::tibble())
  dplyr::bind_rows(lapply(seq_along(seg_start), function(k) {
    b <- as.POSIXct(seg_start[[k]], origin = "1970-01-01", tz = "UTC")
    b <- as.POSIXct(as.Date(b, tz = "UTC"), tz = "UTC")
    dplyr::bind_rows(lapply(c("before", "after"), function(side) {
      sel <- if (side == "before") st$time_utc >= b - window_d * 86400 & st$time_utc < b else
        st$time_utc >= b & st$time_utc < b + window_d * 86400
      x <- st[sel, ]
      ng <- x$good & qc_night(x$ast_h)
      dplyr::tibble(breakpoint = as.Date(b), segment_after = as.integer(names(seg_start)[k]), side = side, n = nrow(x),
                    share_good = if (nrow(x)) mean(x$good) else NA_real_,
                    median_midday_m = qc_med(x$mlh_m[x$good & qc_midday(x$ast_h)]),
                    median_night_m = qc_med(x$mlh_m[ng]), share_lowest_night = if (any(ng)) mean(x$lowest[ng]) else NA_real_)
    }))
  }))
}

#' Night MLH with and without a temperature overlap model (STRATfinder bit 8192)
#'
#' Juelich: TUB120017 had no temperature-dependent overlap model (bit 8192, 2021-03 to 2025-11); TUB160062 has one
#' (data from 2026-06-18). Compared on the same calendar months (those with good steps in both groups), at night
#' (AST 22-04 h) and midday.
#'
#' @param st Output of [meteoswiss_qc_steps()].
#' @return Tibble per group (`no_temperature_model`, `temperature_model`): `months`, `n_night`, `median_night_m`,
#'   `p10_night_m`, `p90_night_m`, `mode_night_m` [m, 10 m bins], `share_lowest_night`, `median_midday_m`; empty
#'   without bit 8192.
qc_temp_overlap <- function(st) {
  if (!any(st$no_temp_overlap)) return(dplyr::tibble())
  g <- ifelse(st$no_temp_overlap, "no_temperature_model", "temperature_model")
  mo <- intersect(unique(st$month[g == "no_temperature_model" & st$good]), unique(st$month[g == "temperature_model" & st$good]))
  dplyr::bind_rows(lapply(c("no_temperature_model", "temperature_model"), function(k) {
    x <- st[g == k & st$month %in% mo, ]
    ng <- x$good & qc_night(x$ast_h)
    v <- x$mlh_m[ng]
    tb <- table(floor(v / 10) * 10)
    dplyr::tibble(group = k, months = paste(sort(mo), collapse = ","), n_night = length(v), median_night_m = qc_med(v),
                  p10_night_m = if (length(v)) unname(stats::quantile(v, 0.1)) else NA_real_,
                  p90_night_m = if (length(v)) unname(stats::quantile(v, 0.9)) else NA_real_,
                  mode_night_m = if (length(v)) as.numeric(names(tb)[which.max(tb)]) else NA_real_,
                  share_lowest_night = if (any(ng)) mean(x$lowest[ng]) else NA_real_,
                  median_midday_m = qc_med(x$mlh_m[x$good & qc_midday(x$ast_h)]))
  }))
}

#' Proposed exclusion periods (not applied)
#'
#' Spec 13: proposals for G, never applied without approval. Rule: calendar months in which at least `min_window` of
#' the steps carry the window flag (CABAM bit 2 dirty window, STRATfinder bit 4): the remaining good values of those
#' months are suspect (MeteoSwiss: spurious deep layers from a contaminated window).
#'
#' @param mon Output of [qc_monthly()].
#' @param min_window Minimum share of window-flagged steps [0-1].
#' @return Tibble: `year_month`, `share_window`, `share_good`, `reason`.
qc_proposed_exclusions <- function(mon, min_window = 0.3) {
  x <- mon[mon$share_window >= min_window, c("year_month", "share_window", "share_good")]
  x$reason <- sprintf("window flag on %.0f %% of the steps (>= %.0f %%): remaining values suspect", 100 * x$share_window,
                      100 * min_window)
  x
}
