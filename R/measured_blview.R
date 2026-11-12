# Vaisala BL-View L3 mixing-layer height (spec 13 item 5 prototype, spec 14 Part B; CZ-Lnz, Lanzhot grassland CL61,
# ~3 km from the tower). Secondary source: shown and tabulated, never pooled with the E-PROFILE sites.
#
# Versions (detected from the variables, never from the folder name; checked on the samples in
# data/raw/measured/CZ-Lnz/samples/, 2026-11-12):
# - "3.1.1": `L3_OFFLINE__<YYYYMMDDhhmm>_<settings>.nc` (site NCO extracts `mlh_L3_OFFLINE_...`): `time` [s since
#   1970-01-01, int, attribute `Units`], `date_stamp` (text "YYYY-MM-DD hh:mm:ss", equal to `time` in UTC on the
#   sample), `Mean_Layer_Height`, `Mean_Layer_QualityIndex` (1 x n). Some 3.1.1 files carry no MLH (site).
# - "3.1.2": as 3.1.1 plus `location_utc_offset` (0 on the sample) and `Mean_Layer_Calculation_Time` [s since 1970,
#   one value per hour, equal to the first profile of the hour on the sample, hh:00:55].
# - "3.2.2": `l3_<dd.mm.yyyy>.nc` (full L3, and the site extracts `mlh_l3_...`): dimension `time_aerosol` [s since
#   1970-01-01; comment "represents the end of the averaging period"], `mean_mixing_layer_height` [m] and
#   `mean_mixing_layer_quality` ("mean mixing layer quality index (1-3) for the period of one hour"), fill -999; the
#   full files also carry `time_zone` (2 on the 2026-10 samples) and global `sw_version` (instrument firmware).
# ~1-min profiles; the hourly value changes on the first profile after the full hour (hh:00:xx) in every sample.
# Vaisala BL-View User Guide M211185EN-E, "NetCDF L3 File Format" (docs.vaisala.com, read 2026-11-12):
# Mean_Layer_Height / Mean_Layer_QualityIndex "Average mixing layer height / quality index for the period of 1 hour";
# file-name time "UTC time"; `location_utc_offset` "UTC offset of each measurement". Neither document states which
# end of the hour the value refers to nor the meaning of the index values. Indicative only: Muenkel et al. (ISARS
# 2012, Vaisala + KIT slides) give the BL-View gradient-minimum index 1-3 as red / yellow / green, higher for a
# stronger gradient and a larger distance to other minima, lower with clouds, i.e. 3 = best; not confirmed for the
# L3 hourly mean (spec 14: left open).
# Time (spec 14, not proven): stamps look like UTC (3.1.x file names, offset 0), but the 3.2.2 day extracts run
# 23:54 -> 22:59 UTC and the full files carry time_zone = 2. The reader takes an explicit `time_offset_h`; the
# Part C lag test checks it.

#' Detect the BL-View L3 version from the variables of an open file
#'
#' @param nc Open `ncdf4` handle.
#' @return `"3.2.2"`, `"3.1.2"`, `"3.1.1"` or `"no_mlh"` (none of the MLH variables) [character].
blview_version <- function(nc) {
  v <- names(nc$var)
  if (all(c("mean_mixing_layer_height", "mean_mixing_layer_quality") %in% v)) return("3.2.2")
  if (all(c("Mean_Layer_Height", "Mean_Layer_QualityIndex") %in% v)) {
    return(if ("location_utc_offset" %in% v) "3.1.2" else "3.1.1")
  }
  "no_mlh"
}

#' Read a BL-View L3 file (hourly mean mixing-layer height), any version
#'
#' Reads the profile times and the MLH and quality of the detected version ([blview_version()]), converts the stored
#' time to UTC with `time_offset_h` (no default: the time zone is not confirmed) and collapses the per-profile copies
#' to one row per hour. Grouping: the profiles from hh:00:00 to hh:59:59 (they all carry the value that changes on
#' the first profile after the full hour) form the hour `hour_end_utc` = hh:00, read as the end of the averaging hour
#' [hh - 1 h, hh) (spec 14; in 3.1.2 the hourly calculation time is the first profile of the group, so the value is
#' computed at hh:00 from the hour before; an assumption for 3.1.1 and 3.2.2, checked by the Part C lag test). An
#' hour whose value or quality is not constant is kept with `mlh_m` and `quality` `NA` and `consistent = FALSE`
#' (reported, not chosen).
#'
#' @param path Path of the L3 `.nc` file [character].
#' @param time_offset_h Offset of the stored time from UTC [h]: stored time minus `time_offset_h` = UTC (0 if the
#'   stored time is UTC, 1 for CET, 2 for CEST) [numeric].
#' @return Tibble, one row per hour: `hour_end_utc` [POSIXct UTC], `mlh_m` [m; `NA` = fill -999 or inconsistent],
#'   `quality` [integer 1-3, `NA` = fill], `n_profiles` [-], `n_values` [-, distinct (height, quality) pairs in the
#'   hour], `consistent` [logical], `version`, `calc_time_utc` [POSIXct UTC, 3.1.2 only], `time_zone_h` [h, 3.2.2 full
#'   files only]. Zero rows for `"no_mlh"` files. Attribute `blview_meta`: `file`, `version`, `sw_version`,
#'   `serial_number`, `title`, `time_offset_h`, `n_profiles`, `first_profile_utc`, `last_profile_utc`,
#'   `n_inconsistent`, `date_stamp_matches_utc` (3.1.x: all `date_stamp` equal `time` formatted in UTC),
#'   `location_utc_offset` (3.1.2, distinct values).
read_blview_l3 <- function(path, time_offset_h) {
  if (!file.exists(path)) cli::cli_abort("BL-View file {.path {path}} not found.")
  if (missing(time_offset_h) || length(time_offset_h) != 1 || !is.finite(time_offset_h)) {
    cli::cli_abort("Give {.arg time_offset_h} (stored time minus UTC, h): the time zone is to be confirmed by the site.")
  }
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  ver <- blview_version(nc)
  ga <- ncdf4::ncatt_get(nc, 0)
  meta <- list(file = basename(path), version = ver, sw_version = ga$sw_version,
               serial_number = ga$instrument_serial_number, title = ga$title, time_offset_h = time_offset_h,
               n_profiles = 0L, first_profile_utc = NA, last_profile_utc = NA, n_inconsistent = 0L,
               date_stamp_matches_utc = NA, location_utc_offset = NA)
  empty <- dplyr::tibble(hour_end_utc = as.POSIXct(character(), tz = "UTC"), mlh_m = numeric(), quality = integer(),
                         n_profiles = integer(), n_values = integer(), consistent = logical(), version = character(),
                         calc_time_utc = as.POSIXct(character(), tz = "UTC"), time_zone_h = numeric())
  if (ver == "no_mlh") {
    attr(empty, "blview_meta") <- meta
    return(empty)
  }
  fill <- function(x) {
    x <- as.numeric(as.vector(x))
    x[x == -999] <- NA_real_
    x
  }
  get <- function(v) as.vector(ncdf4::ncvar_get(nc, v))
  if (ver == "3.2.2") {
    if (is.null(nc$dim$time_aerosol)) cli::cli_abort("{.path {path}}: dimension {.val time_aerosol} not found.")
    u <- nc$dim$time_aerosol$units
    if (!grepl("^seconds since 1970-01-01 00:00:00", u)) cli::cli_abort("{.path {path}}: unexpected time units {.val {u}}.")
    t_raw <- as.numeric(nc$dim$time_aerosol$vals)
    h <- fill(get("mean_mixing_layer_height"))
    q <- fill(get("mean_mixing_layer_quality"))
    tz <- if ("time_zone" %in% names(nc$var)) fill(get("time_zone")) else rep(NA_real_, length(t_raw))
    calc <- rep(NA_real_, length(t_raw))
  } else {
    if (!"time" %in% names(nc$var)) cli::cli_abort("{.path {path}}: variable {.val time} not found.")
    ua <- ncdf4::ncatt_get(nc, "time", "Units")$value
    if (!is.character(ua) || !grepl("^seconds since 1970-01-01 00:00:00", ua)) {
      cli::cli_abort("{.path {path}}: unexpected time units {.val {ua}}.")
    }
    t_raw <- as.numeric(get("time"))
    h <- fill(get("Mean_Layer_Height"))
    q <- fill(get("Mean_Layer_QualityIndex"))
    tz <- rep(NA_real_, length(t_raw))
    calc <- if ("Mean_Layer_Calculation_Time" %in% names(nc$var)) fill(get("Mean_Layer_Calculation_Time")) else
      rep(NA_real_, length(t_raw))
    if ("date_stamp" %in% names(nc$var)) {
      ds <- get("date_stamp")
      meta$date_stamp_matches_utc <- all(ds == format(as.POSIXct(t_raw, origin = "1970-01-01", tz = "UTC"),
                                                      "%Y-%m-%d %H:%M:%S"))
    }
    if (ver == "3.1.2") meta$location_utc_offset <- unique(fill(get("location_utc_offset")))
  }
  if (length(h) != length(t_raw) || length(q) != length(t_raw)) {
    cli::cli_abort("{.path {path}}: MLH variables not on the profile time ({length(h)} vs {length(t_raw)}).")
  }
  if (length(t_raw) == 0) {
    attr(empty, "blview_meta") <- meta
    return(empty)
  }
  t <- t_raw - time_offset_h * 3600
  hr <- floor(t / 3600) * 3600
  key <- paste(h, q)
  grp <- match(hr, unique(hr))
  n_val <- as.integer(tapply(key, grp, function(k) length(unique(k))))
  first <- !duplicated(grp)
  cons <- n_val == 1L
  mlh <- h[first]
  qq <- as.integer(q[first])
  mlh[!cons] <- NA_real_
  qq[!cons] <- NA_integer_
  calc_t <- calc[first] - time_offset_h * 3600
  out <- dplyr::tibble(hour_end_utc = as.POSIXct(hr[first], origin = "1970-01-01", tz = "UTC"), mlh_m = mlh,
                       quality = qq, n_profiles = tabulate(grp), n_values = n_val, consistent = cons, version = ver,
                       calc_time_utc = as.POSIXct(calc_t, origin = "1970-01-01", tz = "UTC"),
                       time_zone_h = tz[first])
  out <- out[order(out$hour_end_utc), ]
  meta$n_profiles <- length(t)
  meta$first_profile_utc <- as.POSIXct(min(t), origin = "1970-01-01", tz = "UTC")
  meta$last_profile_utc <- as.POSIXct(max(t), origin = "1970-01-01", tz = "UTC")
  meta$n_inconsistent <- sum(!cons)
  attr(out, "blview_meta") <- meta
  out
}

#' Read every BL-View L3 file under a folder
#'
#' @param dir Folder searched recursively for `.nc` files [character].
#' @param time_offset_h As in [read_blview_l3()] [h].
#' @return List: `hours` (bound output of [read_blview_l3()] with `file`), `files` (tibble per file: `file`, `path`,
#'   `version`, `sw_version`, `n_profiles`, `first_profile_utc`, `last_profile_utc`, `n_hours`, `n_hours_valid`,
#'   `n_inconsistent`, `date_stamp_matches_utc`, `location_utc_offset`, `error`). A file that cannot be read is
#'   listed with its error and skipped.
read_blview_dir <- function(dir, time_offset_h) {
  f <- sort(list.files(dir, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE))
  if (length(f) == 0) cli::cli_abort("No BL-View .nc file under {.path {dir}}.")
  res <- lapply(f, function(p) tryCatch(read_blview_l3(p, time_offset_h), error = function(e) e))
  files <- dplyr::bind_rows(lapply(seq_along(f), function(i) {
    r <- res[[i]]
    if (inherits(r, "error")) {
      return(dplyr::tibble(file = basename(f[i]), path = f[i], version = NA_character_,
                           error = gsub("\\s+", " ", conditionMessage(r))))
    }
    m <- attr(r, "blview_meta")
    dplyr::tibble(file = m$file, path = f[i], version = m$version,
                  sw_version = if (is.null(m$sw_version)) NA_character_ else as.character(m$sw_version),
                  n_profiles = m$n_profiles, first_profile_utc = m$first_profile_utc,
                  last_profile_utc = m$last_profile_utc, n_hours = nrow(r), n_hours_valid = sum(!is.na(r$mlh_m)),
                  n_inconsistent = m$n_inconsistent, date_stamp_matches_utc = m$date_stamp_matches_utc,
                  location_utc_offset = paste(m$location_utc_offset, collapse = " "), error = NA_character_)
  }))
  hours <- dplyr::bind_rows(lapply(seq_along(f), function(i) {
    r <- res[[i]]
    if (inherits(r, "error") || nrow(r) == 0) return(NULL)
    dplyr::mutate(r, file = basename(f[i]))
  }))
  list(hours = hours, files = files)
}

#' One row per hour across BL-View files and versions
#'
#' The same hour can appear in two files (overlapping day files, or a period reprocessed with a newer version). Rule
#' (spec 14; flagged for G): keep the newest version (3.2.2 > 3.1.2 > 3.1.1); within a version, the row with a valid
#' value, then the first file. The overlaps are counted, with the number where the values differ.
#'
#' @param hours `hours` of [read_blview_dir()].
#' @return List: `hours` (one row per `hour_end_utc`), `overlap` (tibble: `version_kept`, `version_dropped`, `n_hours`,
#'   `n_differ` [-, both valid and different height or quality]).
blview_dedup_hours <- function(hours) {
  rank <- match(hours$version, c("3.1.1", "3.1.2", "3.2.2"))
  o <- order(hours$hour_end_utc, -rank, is.na(hours$mlh_m), hours$file)
  h <- hours[o, ]
  dup <- duplicated(h$hour_end_utc)
  kept <- h[!dup, ]
  drop <- h[dup, ]
  ov <- if (nrow(drop) == 0) {
    dplyr::tibble(version_kept = character(), version_dropped = character(), n_hours = integer(), n_differ = integer())
  } else {
    k <- kept[match(drop$hour_end_utc, kept$hour_end_utc), ]
    differ <- !is.na(k$mlh_m) & !is.na(drop$mlh_m) & (k$mlh_m != drop$mlh_m | k$quality != drop$quality)
    x <- dplyr::tibble(version_kept = k$version, version_dropped = drop$version, differ = differ)
    dplyr::summarise(dplyr::group_by(x, version_kept, version_dropped), n_hours = dplyr::n(),
                     n_differ = sum(differ), .groups = "drop")
  }
  list(hours = kept, overlap = ov)
}

#' Report per BL-View version
#'
#' Spec 14 B3: per version, date range, days with any valid MLH, share of valid hours, quality distribution, and the
#' files without MLH.
#'
#' @param r Output of [read_blview_dir()] (hours before [blview_dedup_hours()]).
#' @return Tibble per version (`"no_mlh"`: files without the MLH variables; `"error"`: unreadable files): `version`,
#'   `n_files`, `n_files_no_valid` [-, files without any valid hour], `first_utc`, `last_utc` [POSIXct UTC, first and
#'   last valid hour end], `n_hours`, `n_hours_valid`, `share_valid` [0-1], `n_days_valid` [-], `n_q1`, `n_q2`, `n_q3`
#'   [-, valid hours by quality], `n_inconsistent` [-].
blview_version_table <- function(r) {
  f <- r$files
  f$version[is.na(f$version)] <- "error"
  if (!"n_hours_valid" %in% names(f)) f$n_hours_valid <- NA_integer_
  h <- r$hours
  vs <- intersect(c("3.1.1", "3.1.2", "3.2.2", "no_mlh", "error"), unique(f$version))
  dplyr::bind_rows(lapply(vs, function(v) {
    hv <- if (is.null(h) || nrow(h) == 0) NULL else h[h$version == v, ]
    ok <- if (is.null(hv)) logical() else !is.na(hv$mlh_m)
    dplyr::tibble(
      version = v, n_files = sum(f$version == v),
      n_files_no_valid = sum(f$version == v & !(f$n_hours_valid %in% seq_len(1e6))),
      first_utc = if (length(ok) && any(ok)) min(hv$hour_end_utc[ok]) else as.POSIXct(NA, tz = "UTC"),
      last_utc = if (length(ok) && any(ok)) max(hv$hour_end_utc[ok]) else as.POSIXct(NA, tz = "UTC"),
      n_hours = length(ok), n_hours_valid = sum(ok), share_valid = if (length(ok)) mean(ok) else NA_real_,
      n_days_valid = if (length(ok)) length(unique(as.Date(hv$hour_end_utc[ok] - 1))) else 0L,
      n_q1 = sum(ok & hv$quality %in% 1L), n_q2 = sum(ok & hv$quality %in% 2L), n_q3 = sum(ok & hv$quality %in% 3L),
      n_inconsistent = if (length(ok)) sum(!hv$consistent) else 0L
    )
  }))
}

#' Standardize the BL-View L3 MLH (`vaisala_blview_mlh`, secondary)
#'
#' Hourly interval values [hour_end - 1 h, hour_end) (spec 14 B2; the hour convention of [read_blview_l3()]),
#' `method_class = aerosol`, all quality values kept (no filtering until G decides): `pblh_m` = the hourly BL-View
#' value, `qc_flag` = `NA`, quality in `x_quality`. Exclusion periods of the ICOS site and of the source are applied
#' and counted. Internal use only.
#'
#' @param hours `hours` of [blview_dedup_hours()].
#' @param site_id ICOS station [character].
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h ICOS site UTC offset [h] (for the exclusion periods).
#' @param time_offset_h Offset applied when reading [h] (recorded in the sidecar).
#' @param measured_exclusions Output of [read_measured_exclusions()]; `NULL` = none (reported).
#' @param versions Output of [blview_version_table()] (summarised in the sidecar) [-].
#' @return List: `data` (standard table: `n_samples` = profiles in the hour; `x_quality` [1-3], `x_version`,
#'   `x_consistent`, `x_n_values`, `x_time_zone_h` [h, file variable, 3.2.2 full files]), `meta`, `n_excluded`,
#'   `n_excluded_measured` [-].
standardize_blview_mlh <- function(hours, site_id, exclusions, utc_offset_h, time_offset_h, measured_exclusions = NULL,
                                   versions = NULL) {
  src <- "vaisala_blview_mlh"
  if (anyDuplicated(hours$hour_end_utc)) cli::cli_abort("BL-View: duplicated hours (use {.fn blview_dedup_hours}).")
  d <- dplyr::tibble(timestamp_start_utc = hours$hour_end_utc - 3600, timestamp_end_utc = hours$hour_end_utc,
                     pblh_m = hours$mlh_m, source_id = src, qc_flag = NA_integer_, uncertainty_m = NA_real_,
                     n_samples = as.integer(hours$n_profiles), x_quality = hours$quality, x_version = hours$version,
                     x_consistent = hours$consistent, x_n_values = hours$n_values, x_time_zone_h = hours$time_zone_h)
  d <- d[order(d$timestamp_start_utc), ]
  d <- flag_exclusions(d, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")
  ex <- drop_excluded_rows(d, src)
  d <- ex$data
  d$excl_flag <- NULL
  d$excl_reason <- NULL
  mx <- apply_measured_exclusions(d, site_id, src, measured_exclusions)
  d <- mx$data
  vt <- if (is.null(versions)) "" else paste0(
    "versions: ", paste0(versions$version, " (", versions$n_files, " files, ", versions$n_hours_valid,
                         " valid hours)", collapse = "; "), "; ")
  meta <- standard_meta(
    site_id, src, "aerosol",
    instrument = paste("Vaisala CL61 lidar ceilometer with BL-View (Lanzhot grassland tower, about 3 km from the",
                       "CZ-Lnz tower; site, 2026-11)"),
    retrieval_algorithm = paste("Vaisala BL-View L3 hourly mean mixing-layer height and quality index 1-3 (BL-View",
                                "User Guide M211185EN-E, NetCDF L3 File Format); computed by BL-View, not retrieved in",
                                "this project; meaning of the quality values not documented (indicative: 3 = best,",
                                "Muenkel et al. ISARS 2012); all quality values kept (spec 14)"),
    height_reference = "agl (as given by BL-View)",
    instrument_altitude_m = "unknown",
    time_reference = "interval_mean",
    original_resolution = "hourly value repeated on ~1-min profiles; one row per hour",
    original_qc = "BL-View quality index 1-3 (x_quality); no filtering applied",
    qc_flag_meaning = list(`NA` = "no filtering: quality in x_quality, rule to be decided by G (spec 14)"),
    processing_history = paste0(
      vt, "read by read_blview_dir() / read_blview_l3() with time_offset_h = ", time_offset_h, " (stored time ",
      "assumed UTC; not proven, spec 14 Part C lag test); hour = profiles hh:00-hh:59 -> interval [hh - 1 h, hh) ",
      "(end of the averaging hour); hours with a changing value kept as NA (x_consistent FALSE); overlaps between ",
      "files: newest version kept (blview_dedup_hours()); rows in exclusion periods dropped (", ex$n_excluded, "); ",
      mx$history, "; standardized by standardize_blview_mlh() (spec 14)."),
    contact = "CzechGlobe (CZ-Lnz)",
    licence = "internal use only: measured data are not released (spec 14); to be confirmed with the provider",
    citation = "Vaisala BL-View User Guide M211185EN-E (docs.vaisala.com)",
    data_owner = "CzechGlobe",
    processor = "Vaisala BL-View (versions in x_version)",
    known_limitations = paste("secondary source: manufacturer algorithm, never pooled with the E-PROFILE sites;",
                              "time zone and hour convention not proven")
  )
  list(data = d, meta = meta, n_excluded = ex$n_excluded, n_excluded_measured = mx$n_excluded)
}

#' Re-read BL-View hour stamps as local clock time (diagnostic option, fix-14c)
#'
#' Hypothesis (G 2026-11-12, from the cluster lag test): the 3.2.2 stamps (files with `time_zone` = 2) are local time
#' of `tz` (e.g. Europe/Prague, CET / CEST with summer time), stored as if they were UTC. The stored hour end is
#' formatted as a clock time and re-read in `tz`. Hour ends that do not exist in `tz` (the spring-forward hour) give
#' `NA`; in the repeated autumn hour, R resolves the clock time to one of the two instants, so two stored hours can
#' map to the same UTC hour (the caller keeps the first). Diagnostics only: the standard file keeps the stored time.
#'
#' @param t Stored hour ends [POSIXct, read as UTC].
#' @param tz Olson time zone of the stamps [character], e.g. `"Europe/Prague"`.
#' @return POSIXct UTC: the hour ends in UTC (`NA` where the clock time does not exist in `tz`).
blview_local_to_utc <- function(t, tz) {
  clock <- format(t, "%Y-%m-%d %H:%M:%S", tz = "UTC")
  out <- as.POSIXct(clock, format = "%Y-%m-%d %H:%M:%S", tz = tz)
  # A clock time that does not exist in tz (spring forward) may come back shifted instead of NA on some platforms:
  # keep only the times that format back to the same clock time.
  back <- format(out, "%Y-%m-%d %H:%M:%S", tz = tz)
  out[is.na(back) | back != clock] <- NA
  attr(out, "tzone") <- "UTC"
  out
}

#' Apply [blview_local_to_utc()] to the rows of one BL-View version of a standard table
#'
#' @param d Standard rows of `vaisala_blview_mlh` (with `x_version`).
#' @param version Version whose stamps are converted [character], e.g. `"3.2.2"`.
#' @param tz Time zone [character].
#' @return List: `data` (rows converted: `timestamp_end_utc` = UTC hour end, `timestamp_start_utc` = end - 1 h; rows
#'   whose stamp does not exist in `tz` dropped; after the conversion, duplicated hours keep the first row), `n_dropped`
#'   [-, nonexistent], `n_duplicated` [-].
blview_convert_version <- function(d, version, tz) {
  i <- d$x_version %in% version
  e <- blview_local_to_utc(d$timestamp_end_utc[i], tz)
  d$timestamp_end_utc[i] <- e
  d$timestamp_start_utc[i] <- e - 3600
  bad <- i & is.na(d$timestamp_end_utc)
  d <- d[!bad, ]
  d <- d[order(d$timestamp_start_utc), ]
  dup <- duplicated(d$timestamp_start_utc)
  list(data = d[!dup, ], n_dropped = sum(bad), n_duplicated = sum(dup))
}
