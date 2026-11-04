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
#'   [gamma_profiles()] for that profile; `flag = "no_profile"` where none is available, `"polar_night"` /
#'   `"polar_day"` on days without sunrise (no morning profile; spec 09 amendment 1: B&G is not run on these
#'   days, spec 02 Amendment 2, so gamma enters only N of Eq. B3 there, from the month's climatology).
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
  if (is.null(opts$profile_hour)) {
    np <- is.na(st$sunrise_utc) & !is.na(st$polar) # no sunrise (a day with sunrise but no sunset keeps its profile)
    out$flag[np] <- paste0("polar_", st$polar[np])
  }
  out
}

#' Monthly gamma climatology
#'
#' A month without a valid daily gamma (e.g. polar night or polar day for the whole month) has `NA` statistics
#' and `n_days = 0`; [gamma_for_rows()] then uses the constant gamma (spec 11 B, `"constant_fallback"`).
#'
#' @param daily Output of [gamma_daily()].
#' @return Tibble per calendar month 1-12: `month`, `median`, `p25`, `p75` [K m-1], `n_days` [-, days with a
#'   valid gamma], `n_flagged` [-, days with a profile but no valid gamma], `n_polar_night`, `n_polar_day`
#'   [-, days without sunrise, hence without a morning profile].
gamma_climatology <- function(daily) {
  m <- as.integer(format(daily$date, "%m"))
  dplyr::bind_rows(lapply(1:12, function(k) {
    g <- daily$gamma[m == k & !is.na(daily$gamma)]
    q <- if (length(g) > 0) stats::quantile(g, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3)
    dplyr::tibble(month = k, median = q[2], p25 = q[1], p75 = q[3], n_days = length(g),
                  n_flagged = sum(m == k & daily$flag %in% c("too_few_levels", "implausible")),
                  n_polar_night = sum(m == k & daily$flag %in% "polar_night"),
                  n_polar_day = sum(m == k & daily$flag %in% "polar_day"))
  }))
}

#' gamma per model row from the ERA5 tables
#'
#' Spec 11 B: `"era5_clim"`, the monthly climatology of the day's month; `"era5_daily"`, the day's morning gamma,
#' else the climatology (`"era5_clim_fallback"`). Rows of a day share one value (`day_id` of [run_pblh_model()]:
#' from the sunrise of the day to the next sunrise, so the morning profile also applies to the following night). A
#' month without a climatology value falls back to the constant `gamma_const` (`"constant_fallback"`, with a
#' warning; G 2026-10-30: no network fallback in the code).
#'
#' @param day_id Day of each row [Date] ([run_pblh_model()]).
#' @param source `"era5_clim"` or `"era5_daily"`.
#' @param daily Output of [gamma_daily()] (`date`, `gamma`).
#' @param clim Output of [gamma_climatology()] (`month`, `median`).
#' @param gamma_const Constant gamma for months without a climatology value [K m-1], scalar or one value per row
#'   (option `gamma` of [pblh_options()]).
#' @return List: `gamma` [K m-1 per row], `source_used` (`era5_clim`, `era5_daily`, `era5_clim_fallback`,
#'   `constant_fallback`; `NA` where `day_id` is `NA`).
gamma_for_rows <- function(day_id, source = c("era5_clim", "era5_daily"), daily, clim, gamma_const = 0.01) {
  source <- match.arg(source)
  day_id <- as.Date(day_id)
  gamma_const <- rep_len(gamma_const, length(day_id))
  m <- as.integer(format(day_id, "%m"))
  g_clim <- clim$median[match(m, clim$month)]
  no_clim <- !is.na(m) & is.na(g_clim)
  if (any(no_clim)) {
    mm <- sort(unique(m[no_clim]))
    # explicit quantity: a numeric vector cannot be the quantity of {?s} (spec 09 amendment 1, SE-Sto)
    cli::cli_warn(c("ERA5 gamma climatology missing for {cli::qty(length(mm))}month{?s} {.val {mm}}.",
                    "i" = "Constant gamma used on {sum(no_clim)} row{?s} ({.field gamma_source_used} = {.val constant_fallback})."))
  }
  g_clim <- ifelse(no_clim, gamma_const, g_clim)
  src_clim <- ifelse(no_clim, "constant_fallback", "era5_clim")
  if (source == "era5_clim") {
    return(list(gamma = g_clim, source_used = ifelse(is.na(day_id), NA_character_, src_clim)))
  }
  g_day <- daily$gamma[match(day_id, as.Date(daily$date))]
  src_fb <- ifelse(no_clim, "constant_fallback", "era5_clim_fallback")
  list(gamma = ifelse(is.na(g_day), g_clim, g_day),
       source_used = ifelse(is.na(day_id), NA_character_, ifelse(is.na(g_day), src_fb, "era5_daily")))
}

#' Constant gamma of one site: median of its daily ERA5 gamma
#'
#' Side result of spec 09 amendment 3 (G 2026-11-04): case `gamma_site` of step 02, one constant per site = the
#' median of the site's valid daily morning gamma ([gamma_daily()], default fit layer, all days of the ERA5 period).
#'
#' @param daily Output of [gamma_daily()] (`<derived_dir>/gamma/<site>_daily.csv`; `gamma` [K m-1], `NA` if flagged).
#' @return Tibble: `gamma` [K m-1, `NA` without a valid day], `n_days` [-, valid days].
gamma_site_constant <- function(daily) {
  g <- daily$gamma[is.finite(daily$gamma)]
  dplyr::tibble(gamma = if (length(g)) stats::median(g) else NA_real_, n_days = length(g))
}

#' Network constant gamma: median of the site medians of the daily ERA5 gamma
#'
#' Side result of spec 09 amendment 3 (G 2026-11-04): case `gamma_med` of step 02, one constant for all stations =
#' the median over the stations of each station's median valid daily morning gamma ([gamma_site_constant()], default
#' fit layer; every station weighs the same). The median of the daily values pooled over every site-day is kept in the
#' table for reference (G 2026-11-04).
#'
#' @param daily_by_site Named list (names = station codes) of [gamma_daily()] tables.
#' @return One-row tibble: `gamma_med` [K m-1, median of the site medians], `gamma_med_pooled` [K m-1, reference],
#'   `n_sites`, `n_days` [-, pooled valid days], `sites` [character, `;`-separated, sorted].
gamma_network_median <- function(daily_by_site) {
  if (length(daily_by_site) == 0 || is.null(names(daily_by_site))) {
    cli::cli_abort("{.arg daily_by_site} must be a named, non-empty list of daily gamma tables.")
  }
  g <- lapply(daily_by_site, function(x) x$gamma[is.finite(x$gamma)])
  g <- g[lengths(g) > 0]
  if (length(g) == 0) cli::cli_abort("No valid daily gamma at any site.")
  dplyr::tibble(gamma_med = stats::median(vapply(g, stats::median, 0)),
                gamma_med_pooled = stats::median(unlist(g, use.names = FALSE)),
                n_sites = length(g), n_days = sum(lengths(g)), sites = paste(sort(names(g)), collapse = ";"))
}

#' Path of the network gamma table
#'
#' Written by `scripts/09_collect.R` ([gamma_network_median()]), read by step 02 (case `gamma_med`) and part of the
#' step-02 inputs (input group `gamma_network`). Spec 09 amendment 3.
#'
#' @param cfg Configuration list with `derived_dir`.
#' @return Path [character]: `<derived_dir>/gamma/network_gamma.csv`.
gamma_network_path <- function(cfg) {
  file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "gamma", "network_gamma.csv")
}

#' Read the network gamma table
#'
#' @param path Output of [gamma_network_path()].
#' @return The one-row table of [gamma_network_median()], or `NULL` if the file does not exist (first pass: case
#'   `gamma_med` not run). Aborts on a malformed file or a non-positive `gamma_med`.
read_gamma_network <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (nrow(x) != 1 || !"gamma_med" %in% names(x) || !(is.finite(x$gamma_med) && x$gamma_med > 0)) {
    cli::cli_abort("{.path {path}}: one row with a positive {.field gamma_med} [K m-1] expected.")
  }
  x
}

#' Write the network gamma table only if its content changed
#'
#' The step-02 inputs hash uses file size and modification time, so an unchanged table is not rewritten (no rerun of
#' step 02 at the next collection).
#'
#' @param x Output of [gamma_network_median()].
#' @param path Output of [gamma_network_path()].
#' @return `TRUE` if written, `FALSE` if unchanged (invisibly).
write_gamma_network <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(x, tmp)
  if (file.exists(path) && identical(unname(tools::md5sum(tmp)), unname(tools::md5sum(path)))) {
    unlink(tmp)
    return(invisible(FALSE))
  }
  ok <- file.copy(tmp, path, overwrite = TRUE)
  unlink(tmp)
  if (!ok) cli::cli_abort("Could not write {.path {path}}.")
  invisible(TRUE)
}
