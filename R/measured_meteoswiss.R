# MeteoSwiss STRATfinder mixed-layer height at Davos (spec 11 Part C).
# Dataset: MeteoSwiss client file "MLH and ABLH from STRATfinder, CHM15k ceilometer Davos (DAV)" (CF-1.8, product
# 1.0, received via the CH-Dav PI 2026-10-29): native ~60 s profiles, `time` in seconds since 1970-01-01 UTC; `mlh`,
# `ablh` in m above ground at the ceilometer (sensor at ground, station 1590 m); `quality_level_mlh` 0 good, 1
# suspect, 2 bad, 9 missing (provider: use 0 only); `qc_flag_mlh` bitmask (8 = at the lowest detectable height,
# <= 30 m, informational). STRATfinder: Kotthaus et al. (2020, Kotthaus2020). Licence CC BY-NC 4.0.

#' Read a MeteoSwiss STRATfinder MLH/ABLH NetCDF file
#'
#' Reads the native profiles. Aborts unless `time` is in seconds since 1970-01-01 UTC, the heights of `mlh` are
#' above ground, and the variables are present.
#'
#' @param path Path of the `.nc` file [character].
#' @return Tibble, one row per profile: `time_utc` [POSIXct UTC], `mlh_m`, `ablh_m` [m above ground at the
#'   ceilometer; `NA` = fill], `quality_mlh`, `quality_ablh` [integer: 0 good, 1 suspect, 2 bad, 9 missing],
#'   `qc_flag_mlh` [integer bitmask], `segment` [integer, homogeneous instrument segment]. Attribute `meteoswiss_meta`:
#'   list with `file`, `title`, `product_version`, `date_created`, `licence`, `acknowledgement`, `references`,
#'   `instrument`, `source`, `institution`, `station_name`, `station_code`, `station_lat` [deg N], `station_lon`
#'   [deg E], `station_alt` [m a.s.l.], `sensor_height` [m agl], `search_upper_limit` [m agl],
#'   `instrument_events` [character].
read_meteoswiss_stratfinder <- function(path) {
  if (!file.exists(path)) cli::cli_abort("MeteoSwiss file {.path {path}} not found.")
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  need <- c("mlh", "ablh", "quality_level_mlh", "quality_level_ablh", "qc_flag_mlh", "homogeneous_segment",
            "station_altitude", "sensor_height", "search_upper_limit", "lat", "lon")
  miss <- setdiff(need, names(nc$var))
  if (length(miss) > 0) cli::cli_abort("{.path {path}}: variable{?s} {.val {miss}} not found.")
  t_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  if (!grepl("^seconds since 1970-01-01( 00:00(:00)?)? UTC$", t_units)) {
    cli::cli_abort("{.path {path}}: unexpected time units {.val {t_units}} (expected seconds since 1970-01-01 UTC).")
  }
  h_ref <- ncdf4::ncatt_get(nc, "mlh", "height_reference")$value
  if (!grepl("above ground", h_ref)) cli::cli_abort("{.path {path}}: mlh height reference {.val {h_ref}} is not above ground.")
  get <- function(v) as.vector(ncdf4::ncvar_get(nc, v))
  num <- function(v) {
    x <- as.numeric(get(v))
    x[is.nan(x)] <- NA_real_
    x
  }
  out <- dplyr::tibble(
    time_utc = as.POSIXct(as.numeric(get("time")), origin = "1970-01-01", tz = "UTC"),
    mlh_m = num("mlh"), ablh_m = num("ablh"),
    quality_mlh = as.integer(get("quality_level_mlh")), quality_ablh = as.integer(get("quality_level_ablh")),
    qc_flag_mlh = as.integer(get("qc_flag_mlh")), segment = as.integer(get("homogeneous_segment"))
  )
  if (is.unsorted(out$time_utc) || anyDuplicated(out$time_utc) > 0) {
    cli::cli_abort("{.path {path}}: times not strictly increasing.")
  }
  ga <- ncdf4::ncatt_get(nc, 0)
  attr(out, "meteoswiss_meta") <- list(
    file = basename(path), title = ga$title, product_version = ga$product_version, date_created = ga$date_created,
    licence = ga$license, acknowledgement = ga$acknowledgement, references = ga$references,
    instrument = ga$instrument, source = ga$source, institution = ga$institution, station_name = ga$station_name,
    station_code = ga$station_code, station_lat = as.numeric(get("lat")), station_lon = as.numeric(get("lon")),
    station_alt = as.numeric(get("station_altitude")), sensor_height = as.numeric(get("sensor_height")),
    search_upper_limit = as.numeric(get("search_upper_limit")), instrument_events = ga$instrument_events
  )
  out
}

#' Aggregate MeteoSwiss STRATfinder profiles to half-hours
#'
#' Spec 11 C: per half-hour [t, t + 30 min) (UTC), the median of the quality-0 MLH values (provider recommendation:
#' quality level 0 only); `n_samples` = number of distinct minutes with a quality-0 MLH (the profile times jitter by
#' about 1 s around 60 s, so two profiles can fall in one minute); `qc_flag` from the coverage `n_samples / 30` as
#' [reobs_qc()] (0 if >= 75 %, 1 if >= 50 %, otherwise `NA` and no value). Rows from the first to the last half-hour
#' with a profile.
#'
#' @param r Output of [read_meteoswiss_stratfinder()].
#' @return Tibble: `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC], `pblh_m` [m agl, median quality-0 MLH
#'   where `qc_flag` is 0 or 1], `qc_flag`, `n_samples` [-, minutes], `x_coverage` [0-1], `x_ablh_m` [m agl, median
#'   quality-0 ABLH, all half-hours with one], `x_share_lowest` [0-1, share of the quality-0 MLH profiles with bit 8,
#'   at the lowest detectable height], `x_segment` [integer, largest homogeneous segment index in the half-hour],
#'   `x_n_profiles` [-, all profiles].
meteoswiss_halfhour <- function(r) {
  t <- as.numeric(r$time_utc)
  hh <- floor(t / 1800) * 1800
  grid <- seq(min(hh), max(hh), by = 1800)
  k <- match(hh, grid)
  n <- length(grid)
  ok <- r$quality_mlh %in% 0L & !is.na(r$mlh_m)
  med <- rep(NA_real_, n)
  m <- tapply(r$mlh_m[ok], k[ok], stats::median)
  med[as.integer(names(m))] <- m
  minute <- floor(t / 60)
  first_min <- ok
  first_min[ok] <- !duplicated(minute[ok]) # one count per minute with a quality-0 value
  n_min <- tabulate(k[first_min], nbins = n)
  ok_a <- r$quality_ablh %in% 0L & !is.na(r$ablh_m)
  ablh <- rep(NA_real_, n)
  a <- tapply(r$ablh_m[ok_a], k[ok_a], stats::median)
  ablh[as.integer(names(a))] <- a
  low <- rep(NA_real_, n)
  s8 <- tapply(bitwAnd(r$qc_flag_mlh[ok], 8L) > 0, k[ok], mean)
  low[as.integer(names(s8))] <- s8
  seg <- rep(NA_integer_, n)
  sg <- tapply(r$segment, k, max)
  seg[as.integer(names(sg))] <- as.integer(sg)
  qc <- reobs_qc(n_min, 30L)
  qc[is.na(med)] <- NA_integer_
  st <- as.POSIXct(grid, origin = "1970-01-01", tz = "UTC")
  dplyr::tibble(timestamp_start_utc = st, timestamp_end_utc = st + 1800,
                pblh_m = ifelse(is.na(qc), NA_real_, med), qc_flag = qc, n_samples = n_min,
                x_coverage = n_min / 30, x_ablh_m = ablh, x_share_lowest = low, x_segment = seg,
                x_n_profiles = tabulate(k, nbins = n))
}

#' Standardize the MeteoSwiss STRATfinder MLH (`meteoswiss_stratfinder_mlh`)
#'
#' Half-hourly interval values [start, start + 30 min) from [meteoswiss_halfhour()], `method_class = aerosol`,
#' heights above ground at the ceilometer (used as is, G 2026-10-30; the station altitude is in the sidecar).
#' Exclusion periods of the ICOS site are applied and counted.
#'
#' @param r Output of [read_meteoswiss_stratfinder()].
#' @param site_id ICOS station paired with the ceilometer [character].
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h ICOS site UTC offset [h] (for the exclusion periods).
#' @return List: `data` (standard table; `uncertainty_m` `NA`; extra columns of [meteoswiss_halfhour()]), `meta`
#'   (sidecar list), `n_excluded` [-].
standardize_meteoswiss_mlh <- function(r, site_id, exclusions, utc_offset_h) {
  mm <- attr(r, "meteoswiss_meta")
  src <- "meteoswiss_stratfinder_mlh"
  h <- meteoswiss_halfhour(r)
  d <- dplyr::tibble(timestamp_start_utc = h$timestamp_start_utc, timestamp_end_utc = h$timestamp_end_utc,
                     pblh_m = h$pblh_m, source_id = src, qc_flag = h$qc_flag, uncertainty_m = NA_real_,
                     n_samples = h$n_samples)
  d <- dplyr::bind_cols(d, h[c("x_coverage", "x_ablh_m", "x_share_lowest", "x_segment", "x_n_profiles")])
  d <- flag_exclusions(d, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")
  ex <- drop_excluded_rows(d, src)
  d <- ex$data
  d$excl_flag <- NULL
  d$excl_reason <- NULL

  meta <- standard_meta(
    site_id, src, "aerosol",
    instrument = paste0(mm$instrument, " Station ", mm$station_name, " (", mm$station_code, "), ", mm$station_lat,
                        " N, ", mm$station_lon, " E, ", mm$station_alt, " m a.s.l.; sensor ", mm$sensor_height,
                        " m above ground."),
    retrieval_algorithm = paste0("STRATfinder (Kotthaus et al. 2020), ", mm$source, "; computed by MeteoSwiss ",
                                 "(not retrieved in this project); daytime search up to ", mm$search_upper_limit,
                                 " m agl, night search 0-700 m"),
    height_reference = "agl",
    instrument_altitude_m = mm$station_alt,
    time_reference = "interval_mean",
    original_resolution = "about 60 s profiles (time_coverage_resolution PT60S, steps 59-61 s)",
    original_qc = paste("MeteoSwiss quality_level_mlh (0 good, 1 suspect, 2 bad, 9 missing; from the bitmask",
                        "qc_flag_mlh); only quality 0 used (provider recommendation)"),
    qc_flag_meaning = list(`0` = "quality-0 MLH in >= 75 % of the 30 minutes",
                           `1` = "quality-0 MLH in >= 50 % and < 75 % of the 30 minutes",
                           `NA` = "no quality-0 MLH, or in < 50 % of the minutes (pblh_m set to NA)"),
    processing_history = paste0(
      "MeteoSwiss file ", mm$file, ", product_version ", mm$product_version, ", created ", mm$date_created,
      "; read by read_meteoswiss_stratfinder(); half-hourly median of the quality-0 MLH profiles ",
      "(meteoswiss_halfhour()); instrument history (file attribute instrument_events, segments in x_segment): ",
      mm$instrument_events, "; rows in exclusion periods dropped (", ex$n_excluded,
      "); standardized by standardize_meteoswiss_mlh() (spec 11)."
    ),
    licence = paste0(mm$licence, " (file attribute 'license'); measured data are not released (project plan, ",
                     "phase 6)"),
    citation = paste0("Kotthaus2020 (docs/references/pblh_literature_2026-10.bib); data provider MeteoSwiss (",
                      mm$institution, "). Acknowledgement: ", mm$acknowledgement)
  )
  list(data = d, meta = meta, n_excluded = ex$n_excluded)
}
