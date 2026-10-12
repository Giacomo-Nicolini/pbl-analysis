# Solar geometry: elevation, sunrise/sunset, day/night per half-hour.
# Equations: NOAA Global Monitoring Laboratory solar calculator
# (https://gml.noaa.gov/grad/solcalc/calcdetails.html; after Meeus 1991,
# "Astronomical Algorithms"). Spec: docs/specs/02_pblh_model.md section 2 and Amendment 1.

deg2rad <- function(x) x * pi / 180
rad2deg <- function(x) x * 180 / pi

#' Solar declination and equation of time (NOAA solar calculator)
#'
#' @param timestamp_utc Times [POSIXct UTC].
#' @return List: `decl` solar declination [deg], `eot` equation of time
#'   [min, apparent minus mean solar time].
solar_params <- function(timestamp_utc) {
  jd <- as.numeric(timestamp_utc) / 86400 + 2440587.5 # Julian day
  tc <- (jd - 2451545) / 36525                         # Julian century (J2000)
  l0 <- (280.46646 + tc * (36000.76983 + tc * 0.0003032)) %% 360 # geom. mean longitude [deg]
  m <- 357.52911 + tc * (35999.05029 - 0.0001537 * tc)          # geom. mean anomaly [deg]
  e <- 0.016708634 - tc * (0.000042037 + 0.0000001267 * tc)     # orbit eccentricity [-]
  ceq <- sin(deg2rad(m)) * (1.914602 - tc * (0.004817 + 0.000014 * tc)) +
    sin(deg2rad(2 * m)) * (0.019993 - 0.000101 * tc) +
    sin(deg2rad(3 * m)) * 0.000289                              # equation of centre [deg]
  true_long <- l0 + ceq
  omega_n <- 125.04 - 1934.136 * tc                             # lunar node [deg]
  app_long <- true_long - 0.00569 - 0.00478 * sin(deg2rad(omega_n)) # apparent longitude [deg]
  obliq0 <- 23 + (26 + (21.448 - tc * (46.815 + tc * (0.00059 - tc * 0.001813))) / 60) / 60
  obliq <- obliq0 + 0.00256 * cos(deg2rad(omega_n))             # corrected obliquity [deg]
  decl <- rad2deg(asin(sin(deg2rad(obliq)) * sin(deg2rad(app_long))))
  y <- tan(deg2rad(obliq / 2))^2
  eot <- 4 * rad2deg(
    y * sin(2 * deg2rad(l0)) - 2 * e * sin(deg2rad(m)) +
      4 * e * y * sin(deg2rad(m)) * cos(2 * deg2rad(l0)) -
      0.5 * y^2 * sin(4 * deg2rad(l0)) - 1.25 * e^2 * sin(2 * deg2rad(m))
  )
  list(decl = decl, eot = eot)
}

#' Solar elevation angle
#'
#' NOAA solar calculator equations, no atmospheric refraction correction.
#'
#' @param timestamp_utc Times [POSIXct UTC].
#' @param lat Latitude [decimal degrees N].
#' @param lon Longitude [decimal degrees E].
#' @return Solar elevation angle [deg] (geometric, negative below the horizon).
solar_elevation <- function(timestamp_utc, lat, lon) {
  sp <- solar_params(timestamp_utc)
  min_utc <- (as.numeric(timestamp_utc) %% 86400) / 60
  tst <- (min_utc + sp$eot + 4 * lon) %% 1440              # true solar time [min]
  ha <- ifelse(tst / 4 < 0, tst / 4 + 180, tst / 4 - 180) # hour angle [deg]
  cosz <- sin(deg2rad(lat)) * sin(deg2rad(sp$decl)) +
    cos(deg2rad(lat)) * cos(deg2rad(sp$decl)) * cos(deg2rad(ha))
  90 - rad2deg(acos(pmin(1, pmax(-1, cosz))))
}

#' Sunrise, sunset and solar noon
#'
#' Hour-angle method of the NOAA solar calculator. Default `zenith_deg = 90`:
#' geometric sunrise/sunset (centre of the solar disc on the horizon, no
#' refraction), the definition used by PBLHmod. `zenith_deg = 90.833`
#' reproduces the NOAA calculator (refraction and disc radius). Declination
#' and equation of time are re-evaluated at the event time (2 iterations).
#'
#' @param date_utc Dates [Date or coercible, UTC calendar date of solar noon].
#' @param lat Latitude [decimal degrees N].
#' @param lon Longitude [decimal degrees E].
#' @param zenith_deg Solar zenith angle defining sunrise/sunset [deg].
#' @return Tibble: `date` [Date], `solar_noon_utc`, `sunrise_utc`,
#'   `sunset_utc` [POSIXct UTC; `NA` without sunrise/sunset], `polar` [`"day"`
#'   (sun always above the threshold), `"night"` (always below) or `NA`].
sun_times <- function(date_utc, lat, lon, zenith_deg = 90) {
  date <- as.Date(date_utc)
  t0 <- as.numeric(as.POSIXct(format(date), format = "%Y-%m-%d", tz = "UTC")) # 00 UTC [s]
  to_time <- function(min) as.POSIXct(t0 + min * 60, origin = "1970-01-01", tz = "UTC")

  # Solar noon [min after 00 UTC]: 720 - 4 lon - EoT (EoT at noon, 2 iterations)
  noon <- 720 - 4 * lon
  for (i in 1:2) noon <- 720 - 4 * lon - solar_params(to_time(noon))$eot

  # Hour angle of the event: cos H = cos(zen) / (cos(lat) cos(decl)) - tan(lat) tan(decl)
  event <- function(sign) {
    tm <- noon
    polar <- rep(NA_character_, length(tm))
    for (i in 1:3) {
      sp <- solar_params(to_time(tm))
      cosh <- cos(deg2rad(zenith_deg)) / (cos(deg2rad(lat)) * cos(deg2rad(sp$decl))) -
        tan(deg2rad(lat)) * tan(deg2rad(sp$decl))
      polar <- ifelse(cosh > 1, "night", ifelse(cosh < -1, "day", NA_character_))
      ha <- rad2deg(acos(pmin(1, pmax(-1, cosh))))
      noon_i <- 720 - 4 * lon - sp$eot
      tm <- noon_i + sign * 4 * ha
    }
    tm[!is.na(polar)] <- NA
    list(t = tm, polar = polar)
  }
  sr <- event(-1)
  ss <- event(+1)
  dplyr::tibble(
    date = date,
    solar_noon_utc = to_time(noon),
    sunrise_utc = to_time(sr$t),
    sunset_utc = to_time(ss$t),
    polar = dplyr::coalesce(sr$polar, ss$polar)
  )
}

#' Day/night flag per averaging interval
#'
#' @param timestamp_start_utc,timestamp_end_utc Interval bounds [POSIXct UTC].
#' @param lat Latitude [decimal degrees N].
#' @param lon Longitude [decimal degrees E].
#' @return Logical: `TRUE` if the geometric solar elevation at the interval
#'   midpoint is > 0 deg.
is_day <- function(timestamp_start_utc, timestamp_end_utc, lat, lon) {
  mid <- timestamp_start_utc + (as.numeric(timestamp_end_utc) - as.numeric(timestamp_start_utc)) / 2
  solar_elevation(mid, lat, lon) > 0
}
