# Free-atmosphere potential-temperature gradient gamma from ERA5 pressure levels (spec 11 Part A; project plan,
# phase 1, gamma options 2-3). gamma enters the Batchvarova & Gryning (1991) slab model (Eq. B5) and N =
# sqrt(g gamma / T) of the neutral Eq. B3 and the B2/B3 switch (Kljun et al. 2015 App. B).
# Profiles: 6-hourly T and geopotential on 13 levels (read_era5_pressure_levels(): theta = T (1000 / p)^(R_d / c_p),
# height = Phi / g), tower cell. The morning profile is used: in summer its lower part is often the residual layer
# of the previous day, i.e. what the growing CBL meets (G 2026-10-30).

#' Options of the ERA5 gamma
#'
#' @param zmin_agl,zmax_agl Fit layer [m above the site; levels with height in `[zmin_agl, zmax_agl]`]. Default
#'   500-3500 m (G 2026-10-30: 3500 m keeps 700 hPa and excludes 600 hPa at the lowland sites, where 700 hPa lies
#'   near 3000 m; with 13 levels and no 775 hPa the layer holds 2-3 levels).
#' @param gamma_min,gamma_max Plausible range [K m-1]; values outside become `NA` (flag `implausible`).
#' @param profile_hour `NULL` (default: the profile closest to and not after sunrise) or a UTC hour (0, 6, 12, 18)
#'   forcing the profile time [h].
#' @param min_levels Minimum number of levels in the fit [-].
#' @return Named list.
gamma_options <- function(zmin_agl = 500, zmax_agl = 3500, gamma_min = 0.001, gamma_max = 0.02,
                          profile_hour = NULL, min_levels = 2L) {
  if (!(is.numeric(zmin_agl) && is.numeric(zmax_agl) && zmin_agl >= 0 && zmax_agl > zmin_agl)) {
    cli::cli_abort("{.field zmin_agl} and {.field zmax_agl} must satisfy 0 <= zmin_agl < zmax_agl [m].")
  }
  if (!is.null(profile_hour) && !(profile_hour %in% c(0, 6, 12, 18))) {
    cli::cli_abort("{.field profile_hour} must be NULL or one of 0, 6, 12, 18 [h UTC].")
  }
  list(zmin_agl = zmin_agl, zmax_agl = zmax_agl, gamma_min = gamma_min, gamma_max = gamma_max,
       profile_hour = profile_hour, min_levels = as.integer(min_levels))
}

#' gamma of one profile
#'
#' Least-squares slope of potential temperature on height over the levels in the fit layer; levels below the
#' ground (height above the site < 0) are always excluded (e.g. 850 hPa at Davos).
#'
#' @param z_agl Level heights above the site [m].
#' @param theta Potential temperature [K].
#' @param p_hPa Pressure of the levels [hPa] (to report the levels used).
#' @param opts [gamma_options()].
#' @return List: `gamma` [K m-1, `NA` if flagged], `gamma_raw` [K m-1, before the plausibility test], `n_levels`
#'   [-], `r2` [-], `levels` [character, e.g. `"925;850;700"`], `flag` (`ok`, `too_few_levels`, `implausible`).
gamma_fit_profile <- function(z_agl, theta, p_hPa, opts = gamma_options()) {
  use <- !is.na(z_agl) & !is.na(theta) & z_agl >= 0 & z_agl >= opts$zmin_agl & z_agl <= opts$zmax_agl
  n <- sum(use)
  lev <- paste(sort(p_hPa[use], decreasing = TRUE), collapse = ";")
  if (n < opts$min_levels) {
    return(list(gamma = NA_real_, gamma_raw = NA_real_, n_levels = n, r2 = NA_real_, levels = lev,
                flag = "too_few_levels"))
  }
  z <- z_agl[use]
  th <- theta[use]
  b <- sum((z - mean(z)) * (th - mean(th))) / sum((z - mean(z))^2)
  r2 <- if (n > 2) stats::cor(z, th)^2 else 1
  ok <- b >= opts$gamma_min && b <= opts$gamma_max
  list(gamma = if (ok) b else NA_real_, gamma_raw = b, n_levels = n, r2 = r2, levels = lev,
       flag = if (ok) "ok" else "implausible")
}

#' gamma of every profile
#'
#' @param pl Output of [read_era5_pressure_levels()] (the tower cell is used: `cell_row = cell_col = 0`).
#' @param altitude_m Site altitude [m a.s.l.] (heights above the site).
#' @param orography_m ERA5 cell orography [m a.s.l.] (`z / g` of the static field), for a second fit reported as
#'   `gamma_orog` (heights above the model surface); `NA` to skip.
#' @param opts [gamma_options()].
#' @return Tibble per profile time: `timestamp_utc` [POSIXct UTC], `gamma`, `gamma_raw` [K m-1], `n_levels`, `r2`,
#'   `levels`, `flag`, `gamma_orog` [K m-1, same layer above the ERA5 orography], `levels_orog`.
gamma_profiles <- function(pl, altitude_m, orography_m = NA_real_, opts = gamma_options()) {
  pl <- pl[pl$cell_row == 0 & pl$cell_col == 0, ]
  pl <- pl[order(pl$timestamp_utc, -pl$p_hPa), ]
  idx <- split(seq_len(nrow(pl)), as.numeric(pl$timestamp_utc))
  res <- lapply(idx, function(i) {
    f <- gamma_fit_profile(pl$zg_m[i] - altitude_m, pl$theta_K[i], pl$p_hPa[i], opts)
    g <- if (is.na(orography_m)) list(gamma = NA_real_, levels = NA_character_) else
      gamma_fit_profile(pl$zg_m[i] - orography_m, pl$theta_K[i], pl$p_hPa[i], opts)
    c(f, gamma_orog = g$gamma, levels_orog = g$levels)
  })
  res <- unname(res)
  dplyr::tibble(
    timestamp_utc = as.POSIXct(as.numeric(names(idx)), origin = "1970-01-01", tz = "UTC"),
    gamma = vapply(res, `[[`, 0, "gamma"), gamma_raw = vapply(res, `[[`, 0, "gamma_raw"),
    n_levels = vapply(res, `[[`, 0L, "n_levels"), r2 = vapply(res, `[[`, 0, "r2"),
    levels = vapply(res, `[[`, "", "levels"), flag = vapply(res, `[[`, "", "flag"),
    gamma_orog = vapply(res, `[[`, 0, "gamma_orog"), levels_orog = vapply(res, `[[`, "", "levels_orog")
  )
}

#' Daily morning gamma
#'
#' For each day, the 6-hourly profile closest to and not after the geometric sunrise (the start of the B&G
#' integration), or the profile at `opts$profile_hour` UTC of that date.
#'
#' @param prof Output of [gamma_profiles()].
#' @param dates Days [Date; local standard time dates, as `day_id` of [run_pblh_model()]].
#' @param lat,lon Site coordinates [decimal degrees] (sunrise, [sun_times()]).
#' @param opts [gamma_options()].
#' @return Tibble per date: `date` [Date], `sunrise_utc`, `profile_utc` [POSIXct UTC], the columns of
#'   [gamma_profiles()] for that profile; `flag = "no_profile"` where none is available.
gamma_daily <- function(prof, dates, lat, lon, opts = gamma_options()) {
  dates <- as.Date(dates)
  st <- sun_times(dates, lat, lon)
  t_prof <- as.numeric(prof$timestamp_utc)
  if (is.null(opts$profile_hour)) {
    sr <- as.numeric(st$sunrise_utc)
    j <- findInterval(sr, t_prof) # last profile <= sunrise
    j[j == 0 | is.na(sr)] <- NA
    # the profile must be the one of the same 6-hour slot: not more than 6 h before sunrise
    j[!is.na(j) & sr - t_prof[pmax(j, 1)] >= 6 * 3600] <- NA
  } else {
    want <- as.numeric(as.POSIXct(dates, tz = "UTC")) + opts$profile_hour * 3600
    j <- match(want, t_prof)
  }
  p <- prof[ifelse(is.na(j), NA_integer_, j), ]
  out <- dplyr::bind_cols(dplyr::tibble(date = dates, sunrise_utc = st$sunrise_utc,
                                        profile_utc = p$timestamp_utc), p[setdiff(names(p), "timestamp_utc")])
  out$flag[is.na(j)] <- "no_profile"
  out
}

#' Monthly gamma climatology
#'
#' @param daily Output of [gamma_daily()].
#' @return Tibble per calendar month 1-12: `month`, `median`, `p25`, `p75` [K m-1], `n_days` [-, days with a
#'   valid gamma], `n_flagged` [-, days with a profile but no valid gamma].
gamma_climatology <- function(daily) {
  m <- as.integer(format(daily$date, "%m"))
  dplyr::bind_rows(lapply(1:12, function(k) {
    g <- daily$gamma[m == k & !is.na(daily$gamma)]
    q <- if (length(g) > 0) stats::quantile(g, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3)
    dplyr::tibble(month = k, median = q[2], p25 = q[1], p75 = q[3], n_days = length(g),
                  n_flagged = sum(m == k & daily$flag %in% c("too_few_levels", "implausible")))
  }))
}

#' gamma per model row from the ERA5 tables
#'
#' Spec 11 B: `"era5_clim"`, the monthly climatology of the day's month; `"era5_daily"`, the day's morning gamma,
#' else the climatology (`"era5_clim_fallback"`). Rows of a day share one value (`day_id` of [run_pblh_model()]:
#' from the sunrise of the day to the next sunrise, so the morning profile also applies to the following night). A
#' month without a climatology value aborts.
#'
#' @param day_id Day of each row [Date] ([run_pblh_model()]).
#' @param source `"era5_clim"` or `"era5_daily"`.
#' @param daily Output of [gamma_daily()] (`date`, `gamma`).
#' @param clim Output of [gamma_climatology()] (`month`, `median`).
#' @return List: `gamma` [K m-1 per row], `source_used` (`era5_clim`, `era5_daily`, `era5_clim_fallback`; `NA`
#'   where `day_id` is `NA`).
gamma_for_rows <- function(day_id, source = c("era5_clim", "era5_daily"), daily, clim) {
  source <- match.arg(source)
  day_id <- as.Date(day_id)
  m <- as.integer(format(day_id, "%m"))
  g_clim <- clim$median[match(m, clim$month)]
  if (any(!is.na(m) & is.na(g_clim))) {
    cli::cli_abort("ERA5 gamma climatology missing for month{?s} {.val {sort(unique(m[!is.na(m) & is.na(g_clim)]))}}.")
  }
  if (source == "era5_clim") {
    return(list(gamma = g_clim, source_used = ifelse(is.na(day_id), NA_character_, "era5_clim")))
  }
  g_day <- daily$gamma[match(day_id, as.Date(daily$date))]
  list(gamma = ifelse(is.na(g_day), g_clim, g_day),
       source_used = ifelse(is.na(day_id), NA_character_, ifelse(is.na(g_day), "era5_clim_fallback", "era5_daily")))
}
