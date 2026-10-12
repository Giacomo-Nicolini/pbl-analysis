# ERA5 boundary-layer height reader (yearly CSV per station, hourly, UTC,
# instantaneous at the hour). Rules: docs/specs/01_metadata_readers.md section 5.

era5_file_regex <- function(site_id) {
  site_re <- gsub("([.\\-])", "\\\\\\1", site_id)
  paste0(
    "^", site_re,
    "__(ERA5|ERA5T)__reanalysis-era5-single-levels__boundary_layer_height__(\\d{4})__lon([+-]\\d+)__lat([+-]\\d+)\\.csv$"
  )
}

#' Locate ERA5 BLH files of a station
#'
#' Accepts only names matching
#' `<site>__(ERA5|ERA5T)__reanalysis-era5-single-levels__boundary_layer_height__<year>__lon<±n>__lat<±n>.csv`;
#' other files in the folder are ignored and listed in a message. If ERA5 and
#' ERA5T exist for the same year and offsets, ERA5 is kept (warning).
#'
#' @param site_id Station code [character].
#' @param cfg Configuration list with `era5_dir` (may contain `{site}`).
#' @return Tibble, one row per file: `path`, `file`, `product` (`ERA5`/`ERA5T`),
#'   `year` [integer], `lon_offset`, `lat_offset` [integer, grid cells from the
#'   site cell].
find_era5_files <- function(site_id, cfg) {
  dir <- resolve_path(cfg$era5_dir, site = site_id)
  all_files <- list.files(dir)
  m <- stringr::str_match(all_files, era5_file_regex(site_id))
  ok <- !is.na(m[, 1])
  if (any(!ok)) {
    cli::cli_inform(c(
      "ERA5 {.val {site_id}}: ignoring {sum(!ok)} file{?s} not matching the BLH naming pattern:",
      stats::setNames(all_files[!ok], rep("*", sum(!ok)))
    ))
  }
  out <- dplyr::tibble(
    path = file.path(dir, all_files[ok]),
    file = all_files[ok],
    product = m[ok, 2],
    year = as.integer(m[ok, 3]),
    lon_offset = as.integer(m[ok, 4]),
    lat_offset = as.integer(m[ok, 5])
  )
  if (nrow(out) == 0) {
    cli::cli_warn("ERA5 {.val {site_id}}: no BLH files found in {.path {dir}}.")
    return(out)
  }

  key <- paste(out$year, out$lon_offset, out$lat_offset)
  both <- key %in% key[out$product == "ERA5"] & out$product == "ERA5T"
  if (any(both)) {
    cli::cli_warn(c(
      "ERA5 {.val {site_id}}: ERA5 and ERA5T both present; keeping ERA5, dropping:",
      stats::setNames(out$file[both], rep("*", sum(both)))
    ))
    out <- out[!both, ]
  }
  out[order(out$lat_offset, out$lon_offset, out$year), ]
}

#' Great-circle distance (haversine)
#'
#' Haversine formula (e.g. Sinnott 1984, Sky & Telescope 68:159) on a sphere of
#' radius 6371 km (mean Earth radius).
#'
#' @param lat1,lon1,lat2,lon2 Coordinates [decimal degrees].
#' @return Distance [km].
haversine_km <- function(lat1, lon1, lat2, lon2) {
  r_earth <- 6371 # km
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  2 * r_earth * asin(pmin(1, sqrt(a)))
}

#' Read ERA5 boundary-layer height
#'
#' Reads the files listed by [find_era5_files()]. Checks: `time` vs `datetime`
#' (warn, use `time`), duplicated timestamps per offset (abort), hourly
#' continuity (report gaps), negative or missing `blh_m` (report). No values
#' are modified or removed.
#'
#' @param files Output of [find_era5_files()].
#' @param offsets `"center"` (only `lon_offset == 0 & lat_offset == 0`) or
#'   `"all"` (all grid cells found).
#' @return Tibble sorted by offsets and time: `timestamp_utc` [POSIXct UTC,
#'   instantaneous value at the hour], `blh_m` boundary-layer height [m above
#'   model ground], `product` [-], `lon_offset`, `lat_offset` [integer, cells],
#'   `grid_lat`, `grid_lon` grid-cell centre [decimal degrees], `site_lat`,
#'   `site_lon` [decimal degrees], `grid_dist_km` site to cell-centre distance [km].
read_era5_blh <- function(files, offsets = c("center", "all")) {
  offsets <- match.arg(offsets)
  if (offsets == "center") {
    files <- files[files$lon_offset == 0 & files$lat_offset == 0, ]
  }
  if (nrow(files) == 0) {
    cli::cli_abort("No ERA5 BLH files to read (offsets = {.val {offsets}}).")
  }

  read_one <- function(i) {
    d <- readr::read_csv(
      files$path[i],
      col_types = readr::cols_only(
        value = readr::col_double(),
        time = readr::col_character(),
        datetime = readr::col_character(),
        latitude = readr::col_double(),
        longitude = readr::col_double(),
        my_lat = readr::col_double(),
        my_lon = readr::col_double()
      ),
      na = c("NA", "", "-9999"),
      progress = FALSE
    )
    t <- lubridate::ymd_hms(d$time, tz = "UTC", truncated = 3, quiet = TRUE)
    if (anyNA(t)) {
      cli::cli_abort("{files$file[i]}: {sum(is.na(t))} unparseable {.field time} value{?s}.")
    }
    if ("datetime" %in% names(d)) {
      t2 <- lubridate::ymd_hms(d$datetime, tz = "UTC", truncated = 3, quiet = TRUE)
      n_diff <- sum(is.na(t2) | t2 != t)
      if (n_diff > 0) {
        cli::cli_warn("{files$file[i]}: {n_diff} row{?s} where {.field time} != {.field datetime}; using {.field time}.")
      }
    }
    dplyr::tibble(
      timestamp_utc = t,
      blh_m = d$value,
      product = files$product[i],
      lon_offset = files$lon_offset[i],
      lat_offset = files$lat_offset[i],
      grid_lat = d$my_lat,
      grid_lon = d$my_lon,
      site_lat = d$latitude,
      site_lon = d$longitude
    )
  }
  d <- dplyr::bind_rows(lapply(seq_len(nrow(files)), read_one))
  d <- d[order(d$lat_offset, d$lon_offset, d$timestamp_utc), ]
  cell <- paste0("lon", sprintf("%+d", d$lon_offset), "__lat", sprintf("%+d", d$lat_offset))

  dup <- duplicated(data.frame(cell, d$timestamp_utc))
  if (any(dup)) {
    cli::cli_abort("ERA5: {sum(dup)} duplicated timestamp{?s} within a grid cell, e.g. {.val {format(utils::head(d$timestamp_utc[dup], 3))}}")
  }

  for (k in unique(cell)) {
    t <- as.numeric(d$timestamp_utc[cell == k])
    dt <- diff(t)
    gaps <- dt != 3600 # hourly
    if (any(gaps)) {
      n_missing <- sum(dt[gaps] / 3600 - 1)
      cli::cli_warn("ERA5 {k}: {sum(gaps)} gap{?s} in the hourly series ({n_missing} missing hour{?s}).")
    }
  }
  n_na <- sum(is.na(d$blh_m))
  n_neg <- sum(d$blh_m < 0, na.rm = TRUE)
  if (n_na + n_neg > 0) {
    cli::cli_warn("ERA5: {n_na} missing and {n_neg} negative {.field blh_m} value{?s} (kept as is).")
  }

  d$grid_dist_km <- haversine_km(d$site_lat, d$site_lon, d$grid_lat, d$grid_lon)
  d
}
