# Site metadata: static table (config/sites.csv), time-varying heights
# (config/sites_timevarying.csv) and per-timestamp measurement / canopy height.
# Rules for zm and canopy height: docs/project_plan.md ("Inputs") and docs/specs/01_metadata_readers.md.

#' Read static site metadata
#'
#' Reads `config/sites.csv` and validates it: `site_id` unique; `lat`, `lon`,
#' `utc_offset_h` non-missing for every row. Other columns may be missing (`NA`).
#'
#' @param file Path to `sites.csv` [character].
#' @return Tibble, one row per site:
#'   `site_id` [-], `lat` [decimal degrees N], `lon` [decimal degrees E],
#'   `altitude_m` [m a.s.l.], `utc_offset_h` [h, local standard time minus UTC],
#'   `zm_m` measurement height [m above ground], `canopy_height_m` [m],
#'   `d_m` displacement height [m], `z0_m` roughness length [m], `igbp` [-],
#'   `slope_perc` [%], `slope_class` [-], `notes` [-].
read_sites <- function(file) {
  sites <- readr::read_csv(
    file,
    na = c("NA", ""),
    col_types = readr::cols(
      site_id = readr::col_character(),
      lat = readr::col_double(),
      lon = readr::col_double(),
      altitude_m = readr::col_double(),
      utc_offset_h = readr::col_double(),
      zm_m = readr::col_double(),
      canopy_height_m = readr::col_double(),
      d_m = readr::col_double(),
      z0_m = readr::col_double(),
      igbp = readr::col_character(),
      slope_perc = readr::col_double(),
      slope_class = readr::col_character(),
      notes = readr::col_character()
    ),
    progress = FALSE
  )
  # readr drops, without a reported problem, a last line that has fewer fields and no final newline (UK-AMo,
  # 2026-10-29): every non-empty data line must give one row.
  n_lines <- sum(nzchar(trimws(readLines(file, warn = FALSE)))) - 1L
  if (nrow(sites) != n_lines) {
    cli::cli_abort("{.file {file}}: {n_lines} data line{?s} but {nrow(sites)} row{?s} read; check the field count of each line (13 fields).")
  }

  dup <- unique(sites$site_id[duplicated(sites$site_id)])
  if (length(dup) > 0) {
    cli::cli_abort("Duplicated {.field site_id} in {.file {file}}: {.val {dup}}")
  }
  required <- c("site_id", "lat", "lon", "utc_offset_h")
  bad <- sites$site_id[!stats::complete.cases(sites[required])]
  if (length(bad) > 0) {
    cli::cli_abort(c(
      "Required metadata ({.field {required}}) missing in {.file {file}}.",
      "x" = "Sites with missing values: {.val {bad}}"
    ))
  }
  sites
}

#' Read time-varying site heights
#'
#' Reads `config/sites_timevarying.csv` (long format). `zm_m` rows are step
#' periods (`date_start` inclusive, `date_end` exclusive, `NA` = open-ended);
#' `canopy_height_m` rows are point measurements. Dates are local standard time
#' calendar dates. Aborts on unknown `variable` values and on overlapping zm periods.
#'
#' @param file Path to `sites_timevarying.csv` [character].
#' @return Tibble: `site_id` [-], `variable` [-, `zm_m` or `canopy_height_m`],
#'   `value` [m], `date_start` [Date, local standard time], `date_end` [Date,
#'   local standard time, `NA` = open-ended], `notes` [-].
read_sites_timevarying <- function(file) {
  tv <- readr::read_csv(
    file,
    na = c("NA", ""),
    col_types = readr::cols(
      site_id = readr::col_character(),
      variable = readr::col_character(),
      value = readr::col_double(),
      date_start = readr::col_date(format = "%Y-%m-%d"),
      date_end = readr::col_date(format = "%Y-%m-%d"),
      notes = readr::col_character()
    ),
    progress = FALSE
  )

  allowed <- c("zm_m", "canopy_height_m")
  bad_var <- setdiff(unique(tv$variable), allowed)
  if (length(bad_var) > 0) {
    cli::cli_abort(c(
      "Unknown {.field variable} in {.file {file}}: {.val {bad_var}}",
      "i" = "Allowed: {.val {allowed}}"
    ))
  }
  if (anyNA(tv$date_start) || anyNA(tv$site_id)) {
    cli::cli_abort("{.field site_id} and {.field date_start} must be non-missing in {.file {file}}.")
  }
  zm <- tv[tv$variable == "zm_m", ]
  for (s in unique(zm$site_id)) {
    check_zm_periods(zm[zm$site_id == s, ], s)
  }
  tv
}

# Abort if zm step periods of one site overlap. Contiguity is not required
# (gaps are allowed and give zm_source = "gap").
check_zm_periods <- function(zm, site_id) {
  zm <- zm[order(zm$date_start), ]
  n <- nrow(zm)
  if (n < 2) return(invisible(TRUE))
  end <- zm$date_end[-n]
  next_start <- zm$date_start[-1]
  overlap <- is.na(end) | end > next_start
  if (any(overlap)) {
    i <- which(overlap)
    cli::cli_abort(c(
      "Overlapping {.field zm_m} periods for {.val {site_id}}.",
      "x" = "Period starting {.val {as.character(zm$date_start[i])}} overlaps the next one."
    ))
  }
  invisible(TRUE)
}

# Local standard time calendar date (00:00 local) -> POSIXct UTC.
local_date_to_utc <- function(date, utc_offset_h) {
  as.POSIXct(format(date, "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC") -
    utc_offset_h * 3600
}

#' Measurement height and canopy height per timestamp
#'
#' Rules (docs/project_plan.md, "Inputs"; spec 01 section 2). Dates in
#' `sites_tv` are local standard time dates starting at 00:00 local, converted
#' to UTC with the site `utc_offset_h`.
#' - zm: step function, valid from `date_start` (inclusive) to `date_end`
#'   (exclusive, `NA` = open-ended) -> `"timevarying"`; before the first period
#'   the first value -> `"extrapolated_before_first"`; between periods or after a
#'   closed last period `NA` -> `"gap"`; no time-varying record -> `sites` value
#'   -> `"static"` (`NA` -> `"missing"`).
#' - canopy height: linear interpolation between point measurements ->
#'   `"interpolated"`; constant before the first / after the last measurement,
#'   and everywhere for a single measurement -> `"held_constant"`; no record ->
#'   `sites` value -> `"static"`; all missing -> `NA`, `"missing"`.
#'
#' @param site_id Station code [character, length 1].
#' @param timestamp_utc Timestamps [POSIXct, UTC].
#' @param sites Output of [read_sites()].
#' @param sites_tv Output of [read_sites_timevarying()].
#' @return Tibble, one row per input timestamp: `timestamp_utc` [POSIXct UTC],
#'   `zm_m` measurement height [m above ground], `zm_source` [-],
#'   `canopy_height_m` [m], `hc_source` [-].
get_site_heights <- function(site_id, timestamp_utc, sites, sites_tv) {
  meta <- sites[sites$site_id == site_id, ]
  if (nrow(meta) != 1) {
    cli::cli_abort("Site {.val {site_id}} not found (or not unique) in site metadata.")
  }
  off <- meta$utc_offset_h
  t <- as.numeric(timestamp_utc)
  n <- length(t)
  tv <- sites_tv[sites_tv$site_id == site_id, ]

  # --- zm: step function ---
  zm_tv <- tv[tv$variable == "zm_m", ]
  zm_tv <- zm_tv[order(zm_tv$date_start), ]
  zm <- rep(NA_real_, n)
  zm_source <- rep(NA_character_, n)
  if (nrow(zm_tv) == 0) {
    zm[] <- meta$zm_m
    zm_source[] <- if (is.na(meta$zm_m)) "missing" else "static"
  } else {
    check_zm_periods(zm_tv, site_id)
    start <- as.numeric(local_date_to_utc(zm_tv$date_start, off))
    end <- as.numeric(local_date_to_utc(zm_tv$date_end, off))
    end[is.na(end)] <- Inf
    idx <- findInterval(t, start)
    before <- !is.na(idx) & idx == 0
    zm[before] <- zm_tv$value[1]
    zm_source[before] <- "extrapolated_before_first"
    after <- !is.na(idx) & idx > 0
    inside <- after
    inside[after] <- t[after] < end[idx[after]]
    zm[inside] <- zm_tv$value[idx[inside]]
    zm_source[inside] <- "timevarying"
    zm_source[after & !inside] <- "gap"
  }

  # --- canopy height: linear interpolation of point measurements ---
  hc_tv <- tv[tv$variable == "canopy_height_m" & !is.na(tv$value), ]
  hc <- rep(NA_real_, n)
  hc_source <- rep(NA_character_, n)
  if (nrow(hc_tv) == 0) {
    hc[] <- meta$canopy_height_m
    hc_source[] <- if (is.na(meta$canopy_height_m)) "missing" else "static"
  } else {
    tm <- as.numeric(local_date_to_utc(hc_tv$date_start, off))
    ok <- !is.na(t)
    if (length(unique(tm)) == 1) {
      hc[ok] <- mean(hc_tv$value)
      hc_source[ok] <- "held_constant"
    } else {
      # rule = 2: constant outside the measured range; ties (same date) averaged
      hc[ok] <- stats::approx(tm, hc_tv$value, xout = t[ok], rule = 2, ties = mean)$y
      within <- t >= min(tm) & t <= max(tm)
      hc_source[ok] <- ifelse(within[ok], "interpolated", "held_constant")
    }
  }

  dplyr::tibble(
    timestamp_utc = timestamp_utc,
    zm_m = zm,
    zm_source = zm_source,
    canopy_height_m = hc,
    hc_source = hc_source
  )
}
