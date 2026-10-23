# Timing of the diurnal cycle of the sensible heat flux: EC vs ERA5 (spec 07, G / Cowork 2026-10-23).
# Each value is placed at the midpoint of its averaging interval: EC half-hour = start + 15 min; ERA5 hourly
# accumulation (hour ending at the timestamp) = timestamp - 30 min.

#' Mean diurnal cycle by season, on interval midpoints
#'
#' @param t_mid Interval midpoints [POSIXct UTC].
#' @param h Sensible heat flux [W m-2, positive upward].
#' @return Tibble: `season`, `minute` [minutes after 00:00 UTC of the midpoint],
#'   `H` [W m-2, mean], `n` [-].
h_mean_cycle <- function(t_mid, h) {
  ok <- !is.na(h) & !is.na(t_mid)
  d <- dplyr::tibble(season = season_of(t_mid[ok]),
                     minute = as.integer(round((as.numeric(t_mid[ok]) %% 86400) / 60)), h = h[ok])
  out <- dplyr::summarise(d, .by = c(season, minute), H = mean(h), n = dplyr::n())
  out[order(out$season, out$minute), ]
}

#' Morning and afternoon crossing times of a mean diurnal cycle
#'
#' For each level `frac * max(H)`: the morning crossing is the last upward
#' crossing before the time of the maximum, the afternoon crossing the first
#' downward crossing after it; linear interpolation between neighbouring values.
#'
#' @param minute Minutes after 00:00 UTC [0, 1440), one cycle.
#' @param H Mean flux at those minutes [W m-2].
#' @param frac Levels as fractions of the maximum [-]; 0 = sign change.
#' @return Tibble: `frac`, `level_Wm2`, `morning_min`, `afternoon_min` [minutes
#'   after 00:00 UTC; `NA` if no crossing].
h_crossings <- function(minute, H, frac = c(0, 0.2)) {
  o <- order(minute)
  m <- minute[o]
  h <- H[o]
  imax <- which.max(h)
  cross <- function(L, idx, up) {
    for (k in idx) {
      a <- h[k]
      b <- h[k + 1]
      if ((up && a < L && b >= L) || (!up && a >= L && b < L)) return(m[k] + (L - a) / (b - a) * (m[k + 1] - m[k]))
    }
    NA_real_
  }
  dplyr::bind_rows(lapply(frac, function(f) {
    L <- f * h[imax]
    mo <- if (imax > 1) cross(L, rev(seq_len(imax - 1)), TRUE) else NA_real_
    af <- if (imax < length(h)) cross(L, imax:(length(h) - 1), FALSE) else NA_real_
    dplyr::tibble(frac = f, level_Wm2 = L, morning_min = mo, afternoon_min = af)
  }))
}

#' Best-fit time lag between two mean diurnal cycles
#'
#' Both cycles are interpolated (periodic, 1-min grid); the lag minimises the
#' mean squared difference `ref(t) - mod(t - lag)`. A positive lag means `mod`
#' is earlier than `ref` (it must be shifted later by `lag` to match).
#'
#' @param ref_min,ref_H Reference cycle: minutes after 00:00 UTC, flux [W m-2].
#' @param mod_min,mod_H Compared cycle.
#' @param lags Candidate lags [min].
#' @return Best lag [min].
h_best_lag <- function(ref_min, ref_H, mod_min, mod_H, lags = -180:180) {
  grid <- 0:1439
  per <- function(x, y) {
    o <- order(x)
    x <- x[o]
    y <- y[o]
    stats::approx(c(x - 1440, x, x + 1440), c(y, y, y), xout = grid)$y
  }
  r <- per(ref_min, ref_H)
  m <- per(mod_min, mod_H)
  sse <- vapply(lags, function(l) mean((r - m[((grid - l) %% 1440) + 1])^2), numeric(1))
  lags[which.min(sse)]
}

#' EC vs ERA5 sensible-heat timing table
#'
#' @param series Named list of tibbles with `t_mid` [POSIXct UTC] and `h`
#'   [W m-2]; must contain `ERA5_H`; the other entries are EC series.
#' @param frac Crossing levels (fractions of the maximum).
#' @return Tibble per season and series: `max_H` [W m-2], `t_max_min`, and for
#'   each level `morning_<p>`, `afternoon_<p>` [min after 00:00 UTC, midpoint
#'   placement], `lag_era5_min` [min, best lag of ERA5 relative to the EC series;
#'   positive = ERA5 earlier; `NA` on the ERA5 row].
h_timing_table <- function(series, frac = c(0, 0.2)) {
  cyc <- lapply(series, function(s) h_mean_cycle(s$t_mid, s$h))
  seasons <- levels(season_of(Sys.time()))
  dplyr::bind_rows(lapply(seasons, function(se) {
    era <- cyc$ERA5_H[cyc$ERA5_H$season == se, ]
    dplyr::bind_rows(lapply(names(cyc), function(nm) {
      cy <- cyc[[nm]][cyc[[nm]]$season == se, ]
      if (nrow(cy) < 4) return(NULL)
      cr <- h_crossings(cy$minute, cy$H, frac)
      row <- dplyr::tibble(season = se, series = nm, max_H = max(cy$H), t_max_min = cy$minute[which.max(cy$H)])
      for (i in seq_len(nrow(cr))) {
        p <- paste0(round(100 * cr$frac[i]), "pct")
        row[[paste0("morning_", p)]] <- cr$morning_min[i]
        row[[paste0("afternoon_", p)]] <- cr$afternoon_min[i]
      }
      row$lag_era5_min <- if (nm == "ERA5_H" || nrow(era) < 4) NA_real_ else h_best_lag(cy$minute, cy$H, era$minute, era$H)
      row
    }))
  }))
}
