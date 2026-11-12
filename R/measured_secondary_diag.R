# Diagnostics of the secondary measured sources (spec 14 Part C): time-convention lag test on clear-day mornings,
# hourly layer summaries. Internal outputs only (tables/internal/14_*, figures/internal/14_*).

#' Hourly UTC bins of an interval-mean series
#'
#' Each value goes to the UTC hour [h, h + 1 h) containing its interval midpoint; the bin value is the mean of the
#' values in it (one per hour for hourly sources, the two half-hours for half-hourly ones).
#'
#' @param t_start,t_end Interval bounds [POSIXct UTC].
#' @param value Values [m].
#' @return Tibble: `bin_utc` [POSIXct UTC, start of the hour], `value` [m, mean of the non-missing values], `n` [-].
hour_bins <- function(t_start, t_end, value) {
  mid <- as.numeric(t_start) + (as.numeric(t_end) - as.numeric(t_start)) / 2
  b <- floor(mid / 3600) * 3600
  ok <- !is.na(value)
  if (!any(ok)) {
    return(dplyr::tibble(bin_utc = as.POSIXct(character(), tz = "UTC"), value = numeric(), n = integer()))
  }
  s <- tapply(value[ok], b[ok], mean)
  n <- tapply(value[ok], b[ok], length)
  dplyr::tibble(bin_utc = as.POSIXct(as.numeric(names(s)), origin = "1970-01-01", tz = "UTC"), value = as.numeric(s),
                n = as.integer(n))
}

#' Hourly bins of an instantaneous series (ERA5)
#'
#' Bin [h, h + 1 h) = mean of the instantaneous values at h and h + 1 h (both required), the same rule as the
#' comparisons with hourly interval means (`measured_time_base()`).
#'
#' @param t Times of the instantaneous values [POSIXct UTC].
#' @param value Values [m].
#' @return Tibble: `bin_utc` [POSIXct UTC], `value` [m], `n` [-, 2].
instant_hour_bins <- function(t, value) {
  tt <- as.numeric(t)
  nxt <- value[match(tt + 3600, tt)]
  ok <- !is.na(value) & !is.na(nxt)
  dplyr::tibble(bin_utc = as.POSIXct(tt[ok], origin = "1970-01-01", tz = "UTC"), value = (value[ok] + nxt[ok]) / 2,
                n = 2L)
}

#' Morning hour bins of given days (sunrise to solar noon)
#'
#' @param days Days [Date; the UTC date of solar noon, equal to the local date in Europe].
#' @param lat,lon Site coordinates [decimal degrees].
#' @return Tibble: `day` [Date], `bin_utc` [POSIXct UTC, start of an hour whose midpoint lies between sunrise and
#'   solar noon], `season` (DJF, MAM, JJA, SON). Polar days and nights give no bins.
morning_bins <- function(days, lat, lon) {
  days <- unique(as.Date(days))
  if (length(days) == 0) {
    return(dplyr::tibble(day = as.Date(character()), bin_utc = as.POSIXct(character(), tz = "UTC"), season = character()))
  }
  st <- sun_times(days, lat, lon)
  st <- st[!is.na(st$sunrise_utc), ]
  out <- lapply(seq_len(nrow(st)), function(i) {
    b0 <- floor(as.numeric(st$sunrise_utc[i]) / 3600) * 3600
    b <- seq(b0, as.numeric(st$solar_noon_utc[i]), by = 3600)
    mid <- b + 1800
    b <- b[mid >= as.numeric(st$sunrise_utc[i]) & mid <= as.numeric(st$solar_noon_utc[i])]
    if (length(b) == 0) return(NULL)
    dplyr::tibble(day = st$date[i], bin_utc = as.POSIXct(b, origin = "1970-01-01", tz = "UTC"))
  })
  out <- dplyr::bind_rows(out)
  if (nrow(out) == 0) {
    return(dplyr::tibble(day = as.Date(character()), bin_utc = as.POSIXct(character(), tz = "UTC"), season = character()))
  }
  out$season <- as.character(season_of(out$bin_utc))
  out
}

#' Lag correlation of the morning growth (time-convention test)
#'
#' Spec 14 Part C: for each lag L (h), the Pearson correlation between the reference in the morning bins (sunrise to
#' solar noon of the clear days) and the measured value in the bin L hours later. A measured series stamped in local
#' time (CET, UTC + 1) instead of UTC shows its best lag at L = +1 (its value for UTC hour h sits in bin h + 1);
#' CEST, +2. Reported, never corrected automatically.
#'
#' @param meas,ref Hour bins of the measured series and of the reference ([hour_bins()], [instant_hour_bins()]).
#' @param mb Morning bins ([morning_bins()]).
#' @param lags Lags [h].
#' @param min_n Minimum number of pairs for a correlation [-].
#' @return Tibble: `season` (each season and `all`), `lag_h`, `n` [-, pairs], `n_days` [-], `r` [-, Pearson; `NA`
#'   below `min_n`].
lag_correlation <- function(meas, ref, mb, lags = -3:3, min_n = 10L) {
  rt <- as.numeric(ref$bin_utc)
  mt <- as.numeric(meas$bin_utc)
  base <- mb
  base$ref <- ref$value[match(as.numeric(mb$bin_utc), rt)]
  groups <- c(list(all = rep(TRUE, nrow(base))),
              stats::setNames(lapply(c("DJF", "MAM", "JJA", "SON"), function(s) base$season == s),
                              c("DJF", "MAM", "JJA", "SON")))
  dplyr::bind_rows(lapply(names(groups), function(g) {
    dplyr::bind_rows(lapply(lags, function(L) {
      x <- base[groups[[g]], ]
      m <- meas$value[match(as.numeric(x$bin_utc) + L * 3600, mt)]
      ok <- !is.na(x$ref) & !is.na(m)
      n <- sum(ok)
      r <- if (n >= min_n && stats::sd(x$ref[ok]) > 0 && stats::sd(m[ok]) > 0) stats::cor(x$ref[ok], m[ok]) else NA_real_
      dplyr::tibble(season = g, lag_h = L, n = n, n_days = length(unique(x$day[ok])), r = r)
    }))
  }))
}

#' Best lag per season
#'
#' @param lc Output of [lag_correlation()] (possibly with grouping columns).
#' @param by Grouping columns besides `season` [character].
#' @return Tibble: `by` columns, `season`, `best_lag_h`, `r_best`, `r_lag0` [-], `n_best` [-]; seasons without any
#'   correlation are dropped.
best_lag <- function(lc, by = character()) {
  lc <- lc[!is.na(lc$r), ]
  if (nrow(lc) == 0) return(dplyr::tibble())
  g <- dplyr::group_by(lc, dplyr::across(dplyr::all_of(c(by, "season"))))
  dplyr::summarise(g, best_lag_h = lag_h[which.max(r)], r_best = max(r), r_lag0 = if (any(lag_h == 0)) r[lag_h == 0] else NA_real_,
                   n_best = n[which.max(r)], .groups = "drop")
}
