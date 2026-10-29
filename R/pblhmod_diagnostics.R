# Diagnostics of PBLHmod runs from their saved output (spec 10, G 2026-10-29: PBLHmod deficit at FR-Gri). No model
# change: the functions only read the output of run_pblh_model() (or its derived csv).

#' Batchvarova-Gryning start height per day
#'
#' Rebuilds the initial height `h0` of the slab model with the rule of [run_pblh_model()] (spec 02 section 6,
#' Amendment 8): the day starts at the first half-hour with `is_day` after a non-day half-hour; `h0` is the
#' diagnostic height `h_diag_m` of the half-hour before, else the nearest valid one within `h0_search` half-hours
#' before.
#'
#' @param out Output of [run_pblh_model()] or its derived csv: complete, sorted 30-min grid with
#'   `timestamp_start_utc` [POSIXct UTC or ISO text], `is_day`, `day_id`, `regime`, `method`, `h_diag_m` [m],
#'   `h_bg_m` [m], `flag_no_h0`, `flag_bg_stopped_gap`.
#' @param h0_search Search window before sunrise [half-hours] (`pblh_constants()$h0_search_halfhours`).
#' @return Tibble per day start: `day_id` [Date], `sunrise_utc` [POSIXct UTC, start of the first daylight
#'   half-hour], `h0_m` [m, `NA` if none], `h0_method` (`B1`-`B4`), `h0_lag_halfhours` [1..h0_search],
#'   `n_convective` [half-hours of the day in the convective window], `first_h_bg_m` [m, `h_bg_m` of the first
#'   window half-hour], `flag_no_h0`, `flag_bg_stopped_gap` [logical, any half-hour of the day].
bg_start_heights <- function(out, h0_search = pblh_constants()$h0_search_halfhours) {
  ts <- as_utc_time(out$timestamp_start_utc)
  n <- length(ts)
  if (n > 1 && any(diff(as.numeric(ts)) != 1800)) cli::cli_abort("{.arg out} must be a complete, sorted 30-min grid.")
  dd <- dplyr::coalesce(as.logical(out$is_day), FALSE)
  s <- which(dd & c(TRUE, !dd[-n]))
  day <- as.Date(out$day_id)
  lag <- vapply(s, function(i) {
    cand <- seq.int(i - 1, i - h0_search)
    cand <- cand[cand >= 1 & !is.na(out$h_diag_m[pmax(cand, 1)])]
    if (length(cand) == 0) NA_integer_ else as.integer(i - cand[1])
  }, integer(1))
  j <- s - lag
  conv <- out$regime %in% "convective"
  by_day <- function(v) tapply(v, day, any)
  first_conv <- tapply(seq_len(n)[conv], day[conv], min)
  dk <- as.character(day[s])
  dplyr::tibble(
    day_id = day[s],
    sunrise_utc = ts[s],
    h0_m = ifelse(is.na(j), NA_real_, out$h_diag_m[pmax(j, 1)]),
    h0_method = ifelse(is.na(j), NA_character_, as.character(out$method[pmax(j, 1)])),
    h0_lag_halfhours = lag,
    n_convective = as.integer(tapply(conv, day, sum)[dk]),
    first_h_bg_m = out$h_bg_m[first_conv[dk]],
    flag_no_h0 = as.logical(by_day(as.logical(out$flag_no_h0) %in% TRUE)[dk]),
    flag_bg_stopped_gap = as.logical(by_day(as.logical(out$flag_bg_stopped_gap) %in% TRUE)[dk])
  )
}

#' Daily clearness index
#'
#' `kt = sum(SW_IN) / sum(SW_IN_POT)` over the half-hours of the day with `SW_IN_POT > 0` and both values present.
#'
#' @param day Day of each half-hour (e.g. `day_id`) [Date].
#' @param sw_in Incoming shortwave radiation [W m-2].
#' @param sw_in_pot Potential (top-of-atmosphere) shortwave radiation [W m-2].
#' @return Tibble per day: `day_id` [Date], `kt` [-], `coverage` [0-1, share of the half-hours with
#'   `SW_IN_POT > 0` that have `SW_IN`], `n_light` [-].
daily_clearness <- function(day, sw_in, sw_in_pot) {
  light <- !is.na(sw_in_pot) & sw_in_pot > 0
  d <- dplyr::tibble(day_id = as.Date(day), sw = sw_in, pot = sw_in_pot)[light, ]
  d <- dplyr::summarise(dplyr::group_by(d, day_id), n_light = dplyr::n(), coverage = mean(!is.na(sw)),
                        kt = sum(sw[!is.na(sw)]) / sum(pot[!is.na(sw)]), .groups = "drop")
  d[c("day_id", "kt", "coverage", "n_light")]
}
