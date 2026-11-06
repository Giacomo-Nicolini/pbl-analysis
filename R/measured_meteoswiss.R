# MeteoSwiss client files of mixed-layer height (spec 11 Part C: Davos; spec 13: five E-PROFILE stations).
# Davos: "MLH and ABLH from STRATfinder, CHM15k ceilometer Davos (DAV)" (CF-1.8, product 1.0, received via the CH-Dav
# PI 2026-10-29). E-PROFILE (spec 13, product 1.0 of 2026-11-05, same layout): STRATfinder files (Juelich JUE,
# Magdeburg MAG; CHM15k, ~60 s) and CABAM files (Kenttarova KEN, Melun MEL, Poitiers POI; CL31, 15-min means).
# Common layout: `time` in seconds since 1970-01-01 UTC; `mlh` in m above ground at the ceilometer (sensor at ground);
# `quality_level_mlh` 0 good, 1 suspect, 2 bad, 9 missing (provider: use 0 only); `qc_flag_mlh` bitmask
# (meteoswiss_qc_bits()). STRATfinder: Kotthaus et al. (2020, Kotthaus2020). CABAM: Kotthaus and Grimmond (2018a, b;
# Kotthaus2018, Kotthaus2018b). Licence: Davos CC BY-NC 4.0; E-PROFILE files "to be confirmed" (internal use only).

#' QC bits of the MeteoSwiss client files
#'
#' Bit table of `qc_flag_mlh` per method (file attributes `flag_masks`, `flag_meanings`, `qc_tests`; spec 11 C,
#' spec 13). STRATfinder: 8192 `no_temperature_overlap_model` is new in the E-PROFILE files (informational, Juelich
#' optical unit TUB120017). CABAM: 2 is `dirty_window` (value removed), 16384 `incomplete_interval` (suspect).
#'
#' @param method `"STRATfinder"` or `"CABAM"`.
#' @return Named integer vector (meaning = bit value).
meteoswiss_qc_bits <- function(method = c("STRATfinder", "CABAM")) {
  method <- match.arg(method)
  if (method == "STRATfinder") {
    b <- c(no_retrieval = 1, instrument_fault = 2, maintenance_or_degraded_instrument = 4,
           at_lowest_detectable_height = 8, at_upper_search_limit = 16, spike = 32, mlh_above_ablh = 64,
           low_daily_coverage = 128, processed_without_cloud_base_fix = 256, precipitation_flag_set = 512,
           overlap_model_mismatch = 1024, overlap_model_not_rederived = 2048, slant_viewing_geometry = 4096,
           no_temperature_overlap_model = 8192)
  } else {
    b <- c(no_retrieval = 1, dirty_window = 2, at_lowest_detectable_height = 8, at_upper_search_limit = 16,
           spike = 32, low_daily_coverage = 128, precipitation_detected = 512, incomplete_interval = 16384)
  }
  vapply(b, as.integer, 1L)
}

# Aborts if the file declares a qc bit outside meteoswiss_qc_bits() (flag_masks / flag_meanings attributes, when
# present) or sets one: a new bit must be reviewed before the file is used.
check_meteoswiss_bits <- function(nc, qc, method, path) {
  known <- meteoswiss_qc_bits(method)
  masks <- ncdf4::ncatt_get(nc, "qc_flag_mlh", "flag_masks")
  if (masks$hasatt) {
    mn <- strsplit(ncdf4::ncatt_get(nc, "qc_flag_mlh", "flag_meanings")$value, " ")[[1]]
    decl <- stats::setNames(as.integer(masks$value), mn)
    bad <- decl[is.na(known[names(decl)]) | known[names(decl)] != decl]
    if (length(bad) > 0) {
      cli::cli_abort("{.path {path}}: qc bit{?s} {.val {paste(bad, names(bad))}} not in meteoswiss_qc_bits('{method}').")
    }
  }
  extra <- bitwAnd(qc[!is.na(qc)], bitwNot(sum(known)))
  if (any(extra != 0)) cli::cli_abort("{.path {path}}: qc_flag_mlh sets unknown bits {.val {unique(extra[extra != 0])}}.")
  invisible(TRUE)
}

# Global attributes and station variables of a MeteoSwiss client file (shared by both readers).
meteoswiss_file_meta <- function(nc, path, method) {
  ga <- ncdf4::ncatt_get(nc, 0)
  get <- function(v) as.numeric(ncdf4::ncvar_get(nc, v))
  chr <- function(x) if (is.null(x)) NA_character_ else as.character(x)
  list(
    file = basename(path), method = method, title = chr(ga$title), product_version = chr(ga$product_version),
    date_created = chr(ga$date_created), licence = chr(ga$license), acknowledgement = chr(ga$acknowledgement),
    references = chr(ga$references), instrument = chr(ga$instrument), source = chr(ga$source),
    institution = chr(ga$institution), data_owner = chr(ga$data_owner), station_name = chr(ga$station_name),
    station_code = chr(ga$station_code), wigos_station_id = chr(ga$wigos_station_id),
    station_lat = get("lat"), station_lon = get("lon"), station_alt = get("station_altitude"),
    sensor_height = get("sensor_height"),
    search_upper_limit = if ("search_upper_limit" %in% names(nc$var)) get("search_upper_limit") else NA_real_,
    instrument_events = chr(ga$instrument_events), known_limitations = chr(ga$known_limitations),
    qc_tests = chr(ga$qc_tests)
  )
}

# Common checks of a MeteoSwiss client file: variables present, time in seconds since 1970-01-01 UTC, mlh above ground.
check_meteoswiss_file <- function(nc, path, need) {
  miss <- setdiff(need, names(nc$var))
  if (length(miss) > 0) cli::cli_abort("{.path {path}}: variable{?s} {.val {miss}} not found.")
  t_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  if (!grepl("^seconds since 1970-01-01( 00:00(:00)?)? UTC$", t_units)) {
    cli::cli_abort("{.path {path}}: unexpected time units {.val {t_units}} (expected seconds since 1970-01-01 UTC).")
  }
  h_ref <- ncdf4::ncatt_get(nc, "mlh", "height_reference")$value
  if (!grepl("above ground", h_ref)) cli::cli_abort("{.path {path}}: mlh height reference {.val {h_ref}} is not above ground.")
  invisible(TRUE)
}

nc_num <- function(nc, v) {
  x <- as.numeric(as.vector(ncdf4::ncvar_get(nc, v)))
  x[is.nan(x)] <- NA_real_
  x
}

#' Read a MeteoSwiss STRATfinder MLH/ABLH NetCDF file
#'
#' Reads the native profiles. Aborts unless `time` is in seconds since 1970-01-01 UTC, the heights of `mlh` are
#' above ground, the variables are present and every qc bit is known ([meteoswiss_qc_bits()]).
#'
#' @param path Path of the `.nc` file [character].
#' @return Tibble, one row per profile: `time_utc` [POSIXct UTC], `mlh_m`, `ablh_m` [m above ground at the
#'   ceilometer; `NA` = fill], `quality_mlh`, `quality_ablh` [integer: 0 good, 1 suspect, 2 bad, 9 missing],
#'   `qc_flag_mlh` [integer bitmask], `segment` [integer, homogeneous instrument segment]. Attribute `meteoswiss_meta`:
#'   list with `file`, `method` (`"STRATfinder"`), `title`, `product_version`, `date_created`, `licence`,
#'   `acknowledgement`, `references`, `instrument`, `source`, `institution`, `data_owner`, `station_name`,
#'   `station_code`, `wigos_station_id`, `station_lat` [deg N], `station_lon` [deg E], `station_alt` [m a.s.l.],
#'   `sensor_height` [m agl], `search_upper_limit` [m agl], `instrument_events`, `known_limitations`, `qc_tests`
#'   [character; `NA` if absent].
read_meteoswiss_stratfinder <- function(path) {
  if (!file.exists(path)) cli::cli_abort("MeteoSwiss file {.path {path}} not found.")
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  check_meteoswiss_file(nc, path, c("mlh", "ablh", "quality_level_mlh", "quality_level_ablh", "qc_flag_mlh",
                                    "homogeneous_segment", "station_altitude", "sensor_height", "search_upper_limit",
                                    "lat", "lon"))
  get <- function(v) as.vector(ncdf4::ncvar_get(nc, v))
  out <- dplyr::tibble(
    time_utc = as.POSIXct(as.numeric(get("time")), origin = "1970-01-01", tz = "UTC"),
    mlh_m = nc_num(nc, "mlh"), ablh_m = nc_num(nc, "ablh"),
    quality_mlh = as.integer(get("quality_level_mlh")), quality_ablh = as.integer(get("quality_level_ablh")),
    qc_flag_mlh = as.integer(get("qc_flag_mlh")), segment = as.integer(get("homogeneous_segment"))
  )
  if (is.unsorted(out$time_utc) || anyDuplicated(out$time_utc) > 0) {
    cli::cli_abort("{.path {path}}: times not strictly increasing.")
  }
  check_meteoswiss_bits(nc, out$qc_flag_mlh, "STRATfinder", path)
  attr(out, "meteoswiss_meta") <- meteoswiss_file_meta(nc, path, "STRATfinder")
  out
}

#' Read a MeteoSwiss CABAM MLH NetCDF file
#'
#' Spec 13: 15-min means of the CABAM MLH (Vaisala CL31; Kotthaus2018, Kotthaus2018b). `time` is the interval
#' centre and `time_bnds` the interval [start, end) (CABAM itself stamps the end). Heights already converted from
#' slant range with cos(`tilt_angle`), above ground. Aborts unless the intervals are 15-min intervals starting on a
#' quarter hour with the time at their centre, strictly increasing, and every qc bit is known
#' ([meteoswiss_qc_bits()]). `alh` (9 slots of further aerosol layers) is not read.
#'
#' @param path Path of the `.nc` file [character].
#' @return Tibble, one row per 15-min step: `interval_start_utc`, `interval_end_utc` [POSIXct UTC], `mlh_m`, `dlh_m`
#'   [m above ground; `NA` = fill], `quality_mlh` [integer: 0 good, 1 suspect, 2 bad, 9 missing], `qc_flag_mlh`
#'   [integer bitmask], `segment` [integer], `n_raw` [-, raw profiles in the step], `window_min` [%, minimum L1
#'   window transmission], `tilt_deg` [deg]. Attribute `meteoswiss_meta` as [read_meteoswiss_stratfinder()]
#'   (`method` = `"CABAM"`, `search_upper_limit` `NA`).
read_meteoswiss_cabam <- function(path) {
  if (!file.exists(path)) cli::cli_abort("MeteoSwiss file {.path {path}} not found.")
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  check_meteoswiss_file(nc, path, c("time_bnds", "mlh", "dlh", "quality_level_mlh", "qc_flag_mlh",
                                    "homogeneous_segment", "n_samples", "window_transmission_min", "tilt_angle",
                                    "station_altitude", "sensor_height", "lat", "lon"))
  get <- function(v) as.vector(ncdf4::ncvar_get(nc, v))
  t <- as.numeric(get("time"))
  b <- matrix(as.numeric(ncdf4::ncvar_get(nc, "time_bnds")), nrow = 2)
  bad <- b[1, ] %% 900 != 0 | b[2, ] - b[1, ] != 900 | t != b[1, ] + 450
  if (any(bad)) {
    cli::cli_abort("{.path {path}}: {sum(bad)} step{?s} not a 15-min interval on the quarter hour with the time at its centre (first {format(as.POSIXct(t[which(bad)[1]], origin = '1970-01-01', tz = 'UTC'))}).")
  }
  if (is.unsorted(t) || anyDuplicated(t) > 0) cli::cli_abort("{.path {path}}: times not strictly increasing.")
  out <- dplyr::tibble(
    interval_start_utc = as.POSIXct(b[1, ], origin = "1970-01-01", tz = "UTC"),
    interval_end_utc = as.POSIXct(b[2, ], origin = "1970-01-01", tz = "UTC"),
    mlh_m = nc_num(nc, "mlh"), dlh_m = nc_num(nc, "dlh"), quality_mlh = as.integer(get("quality_level_mlh")),
    qc_flag_mlh = as.integer(get("qc_flag_mlh")), segment = as.integer(get("homogeneous_segment")),
    n_raw = nc_num(nc, "n_samples"), window_min = nc_num(nc, "window_transmission_min"),
    tilt_deg = nc_num(nc, "tilt_angle")
  )
  check_meteoswiss_bits(nc, out$qc_flag_mlh, "CABAM", path)
  attr(out, "meteoswiss_meta") <- meteoswiss_file_meta(nc, path, "CABAM")
  out
}

#' Read a MeteoSwiss client file (STRATfinder or CABAM, from the file name)
#'
#' @param path `MeteoSwiss_STRATfinder_*.nc` or `MeteoSwiss_CABAM_*.nc` [character].
#' @return Output of [read_meteoswiss_stratfinder()] or [read_meteoswiss_cabam()].
read_meteoswiss_mlh <- function(path) {
  if (grepl("^MeteoSwiss_STRATfinder_", basename(path))) return(read_meteoswiss_stratfinder(path))
  if (grepl("^MeteoSwiss_CABAM_", basename(path))) return(read_meteoswiss_cabam(path))
  cli::cli_abort("{.path {path}}: not a MeteoSwiss_STRATfinder_* or MeteoSwiss_CABAM_* file.")
}

#' Instrument events of a MeteoSwiss client file
#'
#' Parses the global attribute `instrument_events` (`"<start> to <end> [<category>] <description> (<source>)"`,
#' separated by `" | "`; E-PROFILE files: the window-transmission days of `qc_events_stratfinder.csv`, category
#' `intervention`, and the period without a temperature overlap model).
#'
#' @param x Attribute value [character]; `NA`, `""` or `"none documented"` give no rows.
#' @return Tibble: `start_utc`, `end_utc` [POSIXct UTC], `category`, `description` [character].
meteoswiss_events <- function(x) {
  empty <- dplyr::tibble(start_utc = as.POSIXct(character(), tz = "UTC"), end_utc = as.POSIXct(character(), tz = "UTC"),
                         category = character(), description = character())
  if (length(x) != 1 || is.na(x) || !nzchar(x) || grepl("^none", x)) return(empty)
  p <- trimws(strsplit(x, " | ", fixed = TRUE)[[1]])
  # "<start>[ to <end>][ [<category>]] <description>"; dates YYYY-MM-DD[THH:MM] (a breakpoint has no end, the CABAM
  # files no category)
  d <- "(\\d{4}-\\d{2}-\\d{2}(?:T\\d{2}:\\d{2})?)"
  re <- paste0("^", d, "(?: to ", d, ")?(?: \\[([^]]+)\\])? (.*)$")
  ok <- grepl(re, p, perl = TRUE)
  if (!all(ok)) cli::cli_abort("instrument_events: cannot parse {.val {p[!ok][1]}}.")
  tm <- function(s) {
    s[!nzchar(s)] <- NA_character_
    s <- ifelse(nchar(s) == 10, paste0(s, "T00:00"), s)
    as.POSIXct(s, format = "%Y-%m-%dT%H:%M", tz = "UTC")
  }
  g <- function(i) sub(re, paste0("\\", i), p, perl = TRUE)
  cat <- g(3)
  dplyr::tibble(start_utc = tm(g(1)), end_utc = tm(g(2)), category = ifelse(nzchar(cat), cat, NA_character_),
                description = g(4))
}

#' Aggregate MeteoSwiss STRATfinder profiles to half-hours
#'
#' Spec 11 C: per half-hour [t, t + 30 min) (UTC), the median of the quality-0 MLH values (provider recommendation:
#' quality level 0 only); `n_samples` = number of distinct minutes with a quality-0 MLH (the profile times jitter by
#' about 1 s around 60 s, so two profiles can fall in one minute); `qc_flag` from the coverage `n_samples / 30` as
#' [reobs_qc()] (0 if >= 75 %, 1 if >= 50 %, otherwise `NA` and no value). Rows from the first to the last half-hour
#' with a profile. Informational bits (8 lowest detectable height, 8192 no temperature overlap model; spec 13) do not
#' change the quality level: their shares are kept as `x_` columns.
#'
#' @param r Output of [read_meteoswiss_stratfinder()].
#' @return Tibble: `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC], `pblh_m` [m agl, median quality-0 MLH
#'   where `qc_flag` is 0 or 1], `qc_flag`, `n_samples` [-, minutes], `x_coverage` [0-1], `x_ablh_m` [m agl, median
#'   quality-0 ABLH, all half-hours with one], `x_share_lowest` [0-1, share of the quality-0 MLH profiles with bit 8,
#'   at the lowest detectable height], `x_share_no_temp_overlap` [0-1, share of the quality-0 MLH profiles with bit
#'   8192], `x_segment` [integer, largest homogeneous segment index in the half-hour], `x_n_profiles` [-, all
#'   profiles].
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
  share_bit <- function(bit) {
    out <- rep(NA_real_, n)
    s <- tapply(bitwAnd(r$qc_flag_mlh[ok], bit) > 0, k[ok], mean)
    out[as.integer(names(s))] <- s
    out
  }
  seg <- rep(NA_integer_, n)
  sg <- tapply(r$segment, k, max)
  seg[as.integer(names(sg))] <- as.integer(sg)
  qc <- reobs_qc(n_min, 30L)
  qc[is.na(med)] <- NA_integer_
  st <- as.POSIXct(grid, origin = "1970-01-01", tz = "UTC")
  dplyr::tibble(timestamp_start_utc = st, timestamp_end_utc = st + 1800,
                pblh_m = ifelse(is.na(qc), NA_real_, med), qc_flag = qc, n_samples = n_min,
                x_coverage = n_min / 30, x_ablh_m = ablh, x_share_lowest = share_bit(8L),
                x_share_no_temp_overlap = share_bit(8192L), x_segment = seg, x_n_profiles = tabulate(k, nbins = n))
}

#' Aggregate MeteoSwiss CABAM 15-min steps to half-hours
#'
#' Spec 13: a half-hour [t, t + 30 min) (UTC) holds the two 15-min steps whose intervals start at t and t + 15 min
#' (checked against `time_bnds` by [read_meteoswiss_cabam()]). `pblh_m` = mean of the quality-0 steps (a mean of
#' 15-min means: the half-hourly mean when both are good); `qc_flag` 0 if both steps are quality 0, 1 if one, else
#' `NA` and no value (G 2026-11-06). Rows from the first to the last half-hour with a step.
#'
#' @param r Output of [read_meteoswiss_cabam()].
#' @return Tibble: `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC], `pblh_m` [m agl], `qc_flag`,
#'   `n_samples` [-, quality-0 steps, 0-2], `x_coverage` [0-1, `n_samples` / 2], `x_dlh_m` [m agl, mean CABAM
#'   disconnected (elevated) layer height of the steps with one; only the dirty-window removal applies to it],
#'   `x_share_lowest` [0-1, share of the quality-0 steps with bit 8], `x_segment` [integer, largest segment index],
#'   `x_n_steps` [-, steps in the file, 0-2], `x_n_profiles` [-, raw profiles of the quality-0 steps], `x_tilt_deg`
#'   [deg, mean tilt angle of the steps].
cabam_halfhour <- function(r) {
  s <- as.numeric(r$interval_start_utc)
  hh <- floor(s / 1800) * 1800
  grid <- seq(min(hh), max(hh), by = 1800)
  k <- match(hh, grid)
  n <- length(grid)
  n_steps <- tabulate(k, nbins = n)
  if (any(n_steps > 2)) cli::cli_abort("CABAM: more than two 15-min steps in a half-hour.")
  ok <- r$quality_mlh %in% 0L & !is.na(r$mlh_m)
  agg <- function(x, keep, f) {
    out <- rep(NA_real_, n)
    if (any(keep)) {
      a <- tapply(x[keep], k[keep], f)
      out[as.integer(names(a))] <- a
    }
    out
  }
  n_ok <- tabulate(k[ok], nbins = n)
  qc <- ifelse(n_ok == 2L, 0L, ifelse(n_ok == 1L, 1L, NA_integer_))
  st <- as.POSIXct(grid, origin = "1970-01-01", tz = "UTC")
  dplyr::tibble(timestamp_start_utc = st, timestamp_end_utc = st + 1800,
                pblh_m = ifelse(is.na(qc), NA_real_, agg(r$mlh_m, ok, mean)), qc_flag = qc, n_samples = n_ok,
                x_coverage = n_ok / 2, x_dlh_m = agg(r$dlh_m, !is.na(r$dlh_m), mean),
                x_share_lowest = agg(bitwAnd(r$qc_flag_mlh, 8L) > 0, ok, mean),
                x_segment = as.integer(agg(r$segment, rep(TRUE, nrow(r)), max)), x_n_steps = n_steps,
                x_n_profiles = ifelse(n_ok > 0, agg(r$n_raw, ok & !is.na(r$n_raw), sum), 0),
                x_tilt_deg = agg(r$tilt_deg, !is.na(r$tilt_deg), mean))
}

#' Licence text of a MeteoSwiss client file for the sidecar
#'
#' E-PROFILE files carry `license = "TO BE CONFIRMED ..."`: internal use only until the licence is set (provider's
#' condition, MeteoSwiss 2026-11-06; spec 13).
#'
#' @param licence File attribute `license` [character].
#' @return Character.
meteoswiss_licence <- function(licence) {
  if (!is.na(licence) && grepl("^TO BE CONFIRMED", licence)) {
    paste0("to be confirmed: internal use only, contact the provider before sharing results (MeteoSwiss condition, ",
           "2026-11-06). File attribute 'license': ", licence)
  } else {
    paste0(licence, " (file attribute 'license'); measured data are not released (project plan, phase 6)")
  }
}

#' Standardize a MeteoSwiss MLH (`meteoswiss_stratfinder_mlh` or `meteoswiss_cabam_mlh`)
#'
#' Half-hourly interval values [start, start + 30 min) from [meteoswiss_halfhour()] (STRATfinder) or
#' [cabam_halfhour()] (CABAM), `method_class = aerosol`, heights above ground at the ceilometer (used as is, G
#' 2026-10-30; the station altitude is in the sidecar). Exclusion periods of the ICOS site are applied and counted.
#' Sidecar extra fields (spec 13): `data_owner` (instrument operator; Davos: the institution), `processor`
#' (MeteoSwiss, method and version, product version) and `known_limitations` (file attribute).
#'
#' @param r Output of [read_meteoswiss_stratfinder()] or [read_meteoswiss_cabam()].
#' @param site_id ICOS station paired with the ceilometer [character].
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h ICOS site UTC offset [h] (for the exclusion periods).
#' @param measured_exclusions Output of [read_measured_exclusions()] (periods of this source only, spec 13); `NULL` =
#'   none (reported).
#' @return List: `data` (standard table; `uncertainty_m` `NA`; extra columns of the aggregation), `meta` (sidecar
#'   list), `n_excluded` [-, site exclusion periods], `n_excluded_measured` [-, measured-source periods].
standardize_meteoswiss_mlh <- function(r, site_id, exclusions, utc_offset_h, measured_exclusions = NULL) {
  mm <- attr(r, "meteoswiss_meta")
  cabam <- identical(mm$method, "CABAM")
  src <- if (cabam) "meteoswiss_cabam_mlh" else "meteoswiss_stratfinder_mlh"
  h <- if (cabam) cabam_halfhour(r) else meteoswiss_halfhour(r)
  d <- dplyr::tibble(timestamp_start_utc = h$timestamp_start_utc, timestamp_end_utc = h$timestamp_end_utc,
                     pblh_m = h$pblh_m, source_id = src, qc_flag = h$qc_flag, uncertainty_m = NA_real_,
                     n_samples = h$n_samples)
  d <- dplyr::bind_cols(d, h[grep("^x_", names(h), value = TRUE)])
  d <- flag_exclusions(d, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")
  ex <- drop_excluded_rows(d, src)
  d <- ex$data
  d$excl_flag <- NULL
  d$excl_reason <- NULL
  mx <- apply_measured_exclusions(d, site_id, src, measured_exclusions)
  d <- mx$data

  version <- regmatches(mm$source, regexpr("(STRATfinder|CABAM) [0-9]{4}[a-z]?", mm$source))
  if (length(version) == 0) version <- mm$method
  owner <- if (!is.na(mm$data_owner) && nzchar(mm$data_owner)) mm$data_owner else mm$institution
  events <- if (is.na(mm$instrument_events)) "none" else mm$instrument_events
  if (cabam) {
    alg <- paste0("CABAM (Kotthaus and Grimmond 2018a, b), ", mm$source, "; computed by MeteoSwiss (not retrieved in ",
                  "this project); heights = slant range x cos(tilt_angle), tilt = daily median of the L1 inclinometer")
    res <- "15-min means of the raw CL31 profiles (time = interval centre, interval in time_bnds; CABAM stamps the end)"
    oqc <- paste("MeteoSwiss quality_level_mlh (0 good; 1 suspect: bits 16, 32, 128, 16384 incomplete_interval;",
                 "2 bad: bit 2 dirty window, value removed; 9 missing); only quality 0 used (provider recommendation)")
    qcm <- list(`0` = "both 15-min steps of the half-hour quality 0 (mean of the two)",
                `1` = "one 15-min step quality 0 (its value)",
                `NA` = "no quality-0 step (pblh_m set to NA)")
    agg <- "mean of the quality-0 15-min steps in the half-hour (cabam_halfhour()); dlh kept as x_dlh_m, alh not read"
    cite <- "Kotthaus2018 (Part 1) and Kotthaus2018b (Part 2) (docs/references/pblh_literature_2026-10.bib)"
  } else {
    alg <- paste0("STRATfinder (Kotthaus et al. 2020), ", mm$source, "; computed by MeteoSwiss ",
                  "(not retrieved in this project); daytime search up to ", mm$search_upper_limit,
                  " m agl, night search 0-700 m")
    res <- "about 60 s profiles (time_coverage_resolution PT60S, steps 59-61 s)"
    oqc <- paste("MeteoSwiss quality_level_mlh (0 good, 1 suspect, 2 bad, 9 missing; from the bitmask",
                 "qc_flag_mlh; bits 8, 512, 2048, 8192 informational); only quality 0 used (provider recommendation)")
    qcm <- list(`0` = "quality-0 MLH in >= 75 % of the 30 minutes",
                `1` = "quality-0 MLH in >= 50 % and < 75 % of the 30 minutes",
                `NA` = "no quality-0 MLH, or in < 50 % of the minutes (pblh_m set to NA)")
    agg <- "half-hourly median of the quality-0 MLH profiles (meteoswiss_halfhour())"
    cite <- "Kotthaus2020 (docs/references/pblh_literature_2026-10.bib)"
  }
  meta <- standard_meta(
    site_id, src, "aerosol",
    instrument = paste0(mm$instrument, " Station ", mm$station_name, " (", mm$station_code,
                        if (!is.na(mm$wigos_station_id)) paste0(", WIGOS ", mm$wigos_station_id) else "", "), ",
                        mm$station_lat, " N, ", mm$station_lon, " E, ", mm$station_alt, " m a.s.l.; sensor ",
                        mm$sensor_height, " m above ground."),
    retrieval_algorithm = alg,
    height_reference = "agl",
    instrument_altitude_m = mm$station_alt,
    time_reference = "interval_mean",
    original_resolution = res,
    original_qc = oqc,
    qc_flag_meaning = qcm,
    processing_history = paste0(
      "MeteoSwiss file ", mm$file, ", product_version ", mm$product_version, ", created ", mm$date_created,
      "; read by read_meteoswiss_", tolower(mm$method), "(); ", agg, "; instrument history (file attribute ",
      "instrument_events, segments in x_segment): ", events, "; rows in exclusion periods dropped (", ex$n_excluded,
      "); ", mx$history, "; standardized by standardize_meteoswiss_mlh() (specs 11, 13)."
    ),
    licence = meteoswiss_licence(mm$licence),
    citation = paste0(cite, "; data owner ", owner, "; processing MeteoSwiss. Acknowledgement: ", mm$acknowledgement),
    data_owner = owner,
    processor = paste0("MeteoSwiss (", version, "; client product version ", mm$product_version, ")"),
    known_limitations = if (is.na(mm$known_limitations)) NULL else mm$known_limitations
  )
  list(data = d, meta = meta, n_excluded = ex$n_excluded, n_excluded_measured = mx$n_excluded)
}
