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
#' CEST, +2. Reported, never corrected automatically. Lags -3 ... +6 h (fix-14c: the 3.2.2 peak sat at the +3 h edge).
#'
#' @param meas,ref Hour bins of the measured series and of the reference ([hour_bins()], [instant_hour_bins()]).
#' @param mb Morning bins ([morning_bins()]).
#' @param lags Lags [h].
#' @param min_n Minimum number of pairs for a correlation [-].
#' @return Tibble: `season` (each season and `all`), `lag_h`, `n` [-, pairs], `n_days` [-], `r` [-, Pearson; `NA`
#'   below `min_n`].
lag_correlation <- function(meas, ref, mb, lags = -3:6, min_n = 10L) {
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


#' Lag table of measured series against named references
#'
#' Builds the full lag table without data masking (fix-14c: inside `dplyr::mutate()` the correlation column `r`
#' masked the reference name, so the reference column held correlations and nothing was reduced).
#'
#' @param meas_bins Named list (names = series) of measured hour bins ([hour_bins()]).
#' @param ref_bins Named list (names = reference ids, e.g. `pblhmod_default`, `era5_blh_center`) of hour bins.
#' @param mb Morning bins ([morning_bins()]).
#' @param site_id Station code [character].
#' @param lags Lags [h].
#' @return Tibble: `site_id`, `series`, `reference`, `season`, `lag_h`, `n`, `n_days`, `r`.
lag_table <- function(meas_bins, ref_bins, mb, site_id, lags = -3:6) {
  out <- list()
  for (s in names(meas_bins)) for (ref_id in names(ref_bins)) {
    lc <- lag_correlation(meas_bins[[s]], ref_bins[[ref_id]], mb, lags = lags)
    out[[length(out) + 1]] <- dplyr::bind_cols(
      dplyr::tibble(site_id = rep(site_id, nrow(lc)), series = rep(s, nrow(lc)), reference = rep(ref_id, nrow(lc))), lc)
  }
  dplyr::bind_rows(out)
}

#' Best lag per group and season
#'
#' One row per combination of the `by` columns and `season` (fix-14c: written without `dplyr::summarise()`, so that a
#' group can never return several rows).
#'
#' @param lc Output of [lag_table()] or [lag_correlation()].
#' @param by Grouping columns besides `season` [character].
#' @return Tibble: `by` columns, `season`, `best_lag_h` [h], `r_best` [-], `r_lag0` [-, `NA` if lag 0 has no
#'   correlation], `n_best` [-, pairs at the best lag]; groups without any correlation are dropped.
best_lag <- function(lc, by = character()) {
  lc <- lc[!is.na(lc$r), ]
  keys <- c(by, "season")
  if (nrow(lc) == 0) {
    return(dplyr::as_tibble(c(stats::setNames(rep(list(character()), length(keys)), keys),
                              list(best_lag_h = integer(), r_best = numeric(), r_lag0 = numeric(), n_best = integer()))))
  }
  k <- do.call(paste, c(lc[keys], sep = "\r"))
  dplyr::bind_rows(lapply(split(seq_len(nrow(lc)), factor(k, levels = unique(k))), function(i) {
    x <- lc[i, ]
    j <- which.max(x$r)
    dplyr::bind_cols(x[1, keys], dplyr::tibble(best_lag_h = as.integer(x$lag_h[j]), r_best = x$r[j],
                                               r_lag0 = if (any(x$lag_h == 0)) x$r[x$lag_h == 0][1] else NA_real_,
                                               n_best = as.integer(x$n[j])))
  }))
}

#' Labels of the lag-test references
#'
#' @param ids Reference ids [character].
#' @return Character: `PBLHmod`, `ERA5` (from [pblh_sources_style]), else the id.
reference_label <- function(ids) {
  lab <- c(pblhmod_default = "PBLHmod", era5_blh_center = "ERA5")
  ifelse(ids %in% names(lab), lab[ids], ids)
}

#' Window metrics of a measured series against PBLHmod and ERA5
#'
#' Spec 14 Part C (moved from scripts/14c_secondary_diagnostics.R, fix-14c): aerosol window (apparent solar time,
#' `opts$aerosol_window`) and PBLHmod convective, on the measurement time base ([measured_time_base()]: an hourly
#' interval [t, t + 1 h) is compared with the mean of the two PBLHmod half-hours inside it and with the mean of ERA5 at
#' t and t + 1 h, all labelled by the interval start). By quality class where the series carries `x_quality`
#' (BL-View); one row per source otherwise (SkyVUE layers: `pblh_m` set to the layer before the call).
#'
#' @param s Measured standard rows with `pblh_m` set to the series to evaluate.
#' @param label Series name [character].
#' @param mod PBLHmod standard table (with `x_regime`).
#' @param era Named list with `era5_blh_center`.
#' @param lon Site longitude [decimal degrees E].
#' @param opts Output of [measured_compare_options()].
#' @param site_id Station code [character].
#' @return Tibble: `site_id`, `series`, `source` (`pblhmod_default`, `era5_blh_center`), `quality` (`all`, `q1`-`q3`),
#'   `n_days`, the metrics of [pblh_metrics()] (ln(source / measured)); `NULL` without any measured value.
secondary_window_metrics <- function(s, label, mod, era, lon, opts, site_id) {
  s <- s[!is.na(s$timestamp_start_utc), ]
  if (sum(!is.na(s$pblh_m)) == 0) return(NULL)
  step <- stats::median(as.numeric(s$timestamp_end_utc) - as.numeric(s$timestamp_start_utc))
  tb <- measured_time_base(s, mod, era, step_s = step)
  ex <- function(v) v[v$kind %in% "exact", c("timestamp_utc", "pblh_m")]
  prs <- list(pblhmod_default = list(x = tb$mod, ref = ex(tb$meas), mid = tb$mid_s),
              era5_blh_center = list(x = tb$era$era5_blh_center, ref = ex(tb$era_meas), mid = tb$era_mid_s))
  dplyr::bind_rows(lapply(names(prs), function(src_id) {
    pt <- pair_table(prs[[src_id]]$x, prs[[src_id]]$ref, src_id)
    a <- apparent_solar_hour(pt$timestamp_utc + prs[[src_id]]$mid, lon)
    pt <- pt[pt$regime %in% "convective" & a >= opts$aerosol_window[1] & a < opts$aerosol_window[2], ]
    q <- if ("x_quality" %in% names(s)) s$x_quality[match(as.numeric(pt$timestamp_utc), as.numeric(s$timestamp_start_utc))]
    cls <- if (is.null(q)) list(all = rep(TRUE, nrow(pt))) else
      c(list(all = rep(TRUE, nrow(pt))), stats::setNames(lapply(1:3, function(k) q %in% k), paste0("q", 1:3)))
    dplyr::bind_rows(lapply(names(cls), function(k) {
      x <- pt[cls[[k]], ]
      dplyr::bind_cols(dplyr::tibble(site_id = site_id, series = label, source = src_id, quality = k,
                                     n_days = length(unique(as.Date(x$timestamp_utc)))),
                       pblh_metrics(x$mod, x$ref))
    }))
  }))
}
