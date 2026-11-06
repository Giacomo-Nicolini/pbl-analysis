# Vaisala BL-View L3 mixing-layer height (spec 13 item 5; CZ-Lnz, Lanzhot grassland CL61, ~3 km from the tower).
# Prototype on the sample files data/raw/measured/CZ-Lnz/samples/l3_DD.MM.YYYY.nc (2026-10-07, 2026-10-10): one file
# per day, ~1-min profiles; `time_aerosol` [seconds since 1970-01-01 00:00:00.000; per file comment of `time`: end of
# the averaging period]; `mean_mixing_layer_height` [m, int, fill -999] and `mean_mixing_layer_quality` [1-3, fill
# -999], both "for the period of one hour": one value per hour repeated on every profile of the hour (it changes at
# hh:00:50). Open points (site to confirm): the time zone of the stored time (the files carry `time_zone` = 2, i.e.
# possibly CEST, and the coordinates are 0); whether the hourly value of the profiles stamped hh:00-hh:59 belongs to
# [hh, hh + 1 h) or to the hour ending at hh; the meaning of the quality index (1-3). No standard file yet.

#' Read a BL-View L3 file (hourly mean mixing-layer height, prototype)
#'
#' Reads only `time_aerosol`, `mean_mixing_layer_height` and `mean_mixing_layer_quality` (spec 13), converts the
#' stored time to UTC with `time_offset_h` (no default: the site has to confirm the time zone), and collapses the
#' per-profile copies to one row per hour (the hour of the profile times, floored). Aborts if the value or the
#' quality changes within an hour.
#'
#' @param path Path of the L3 `.nc` file [character].
#' @param time_offset_h Offset of the stored time from UTC [h]: stored time minus `time_offset_h` = UTC (0 if the
#'   stored time is UTC, 2 if it is CEST) [numeric].
#' @return Tibble, one row per hour: `hour_start_utc` [POSIXct UTC, floor of the profile times; interval convention to
#'   be confirmed], `mlh_m` [m; `NA` = fill -999], `quality` [integer 1-3, `NA` = fill], `n_profiles` [-]. Attribute
#'   `blview_meta`: `file`, `sw_version`, `serial_number`, `title`, `time_offset_h`.
read_blview_l3 <- function(path, time_offset_h) {
  if (!file.exists(path)) cli::cli_abort("BL-View file {.path {path}} not found.")
  if (missing(time_offset_h) || length(time_offset_h) != 1 || !is.finite(time_offset_h)) {
    cli::cli_abort("Give {.arg time_offset_h} (stored time minus UTC, h): the time zone is to be confirmed by the site.")
  }
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  need <- c("mean_mixing_layer_height", "mean_mixing_layer_quality")
  miss <- setdiff(need, names(nc$var))
  if (length(miss) > 0) cli::cli_abort("{.path {path}}: variable{?s} {.val {miss}} not found.")
  if (is.null(nc$dim$time_aerosol)) cli::cli_abort("{.path {path}}: dimension {.val time_aerosol} not found.")
  u <- nc$dim$time_aerosol$units
  if (!grepl("^seconds since 1970-01-01 00:00:00", u)) cli::cli_abort("{.path {path}}: unexpected time units {.val {u}}.")
  t <- as.numeric(nc$dim$time_aerosol$vals) - time_offset_h * 3600
  fill <- function(x) {
    x <- as.numeric(as.vector(x))
    x[x == -999] <- NA_real_
    x
  }
  h <- fill(ncdf4::ncvar_get(nc, "mean_mixing_layer_height"))
  q <- fill(ncdf4::ncvar_get(nc, "mean_mixing_layer_quality"))
  if (length(h) != length(t) || length(q) != length(t)) cli::cli_abort("{.path {path}}: variables not on {.val time_aerosol}.")
  hr <- floor(t / 3600) * 3600
  key <- paste(h, q)
  n_val <- tapply(key, hr, function(k) length(unique(k)))
  if (any(n_val > 1)) {
    cli::cli_abort("{.path {path}}: the hourly value changes within {sum(n_val > 1)} hour{?s} (first {format(as.POSIXct(as.numeric(names(n_val)[n_val > 1][1]), origin = '1970-01-01', tz = 'UTC'))}).")
  }
  first <- !duplicated(hr)
  ga <- ncdf4::ncatt_get(nc, 0)
  out <- dplyr::tibble(hour_start_utc = as.POSIXct(hr[first], origin = "1970-01-01", tz = "UTC"), mlh_m = h[first],
                       quality = as.integer(q[first]), n_profiles = as.integer(tabulate(match(hr, hr[first]))))
  attr(out, "blview_meta") <- list(file = basename(path), sw_version = ga$sw_version,
                                   serial_number = ga$instrument_serial_number, title = ga$title,
                                   time_offset_h = time_offset_h)
  out
}
