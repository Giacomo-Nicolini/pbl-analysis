# Hard, date-based exclusion periods (config/exclusion_periods.csv).
# Rows are flagged, never dropped silently; every step reports the counts.

#' Read exclusion periods
#'
#' Reads `config/exclusion_periods.csv`. Dates are local standard time calendar
#' dates; `date_end` is inclusive (the whole day is excluded).
#'
#' @param file Path to `exclusion_periods.csv` [character].
#' @return Tibble: `site_id` [-], `date_start` [Date, local standard time],
#'   `date_end` [Date, local standard time, inclusive], `scope` [-], `reason` [-].
read_exclusions <- function(file) {
  ex <- readr::read_csv(
    file,
    na = c("NA", ""),
    col_types = readr::cols(
      site_id = readr::col_character(),
      date_start = readr::col_date(format = "%Y-%m-%d"),
      date_end = readr::col_date(format = "%Y-%m-%d"),
      scope = readr::col_character(),
      reason = readr::col_character()
    ),
    progress = FALSE
  )
  if (anyNA(ex[c("site_id", "date_start", "date_end", "scope")])) {
    cli::cli_abort("{.field site_id}, {.field date_start}, {.field date_end} and {.field scope} must be non-missing in {.file {file}}.")
  }
  bad <- ex$date_end < ex$date_start
  if (any(bad)) {
    cli::cli_abort("{.field date_end} before {.field date_start} in {.file {file}} for {.val {ex$site_id[bad]}}.")
  }
  ex
}

#' Flag records inside exclusion periods
#'
#' Interval per period: from `date_start` 00:00 local standard time (inclusive)
#' to `date_end + 1 day` 00:00 local standard time (exclusive), both converted to
#' UTC with `utc_offset_h`. Rows are flagged, not dropped; the number of flagged
#' records is always reported (also when zero).
#'
#' @param data Data frame with a UTC time column.
#' @param site_id Station code [character, length 1].
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h Site UTC offset [h, local standard time minus UTC].
#' @param time_col Name of the POSIXct UTC column compared with the periods
#'   [character]; default `"timestamp_start_utc"`.
#' @return `data` with two added columns: `excl_flag` [logical] and
#'   `excl_reason` [character, `NA` if not excluded; reasons of overlapping
#'   periods joined with `"; "`].
flag_exclusions <- function(data, site_id, exclusions, utc_offset_h,
                            time_col = "timestamp_start_utc") {
  if (!time_col %in% names(data)) {
    cli::cli_abort("Column {.field {time_col}} not found in {.arg data}.")
  }
  if (length(utc_offset_h) != 1 || is.na(utc_offset_h)) {
    cli::cli_abort("{.arg utc_offset_h} must be a single non-missing value [h].")
  }
  ex <- exclusions[exclusions$site_id == site_id, ]
  known_scopes <- "all"
  bad_scope <- setdiff(unique(ex$scope), known_scopes)
  if (length(bad_scope) > 0) {
    cli::cli_abort(c(
      "Unknown exclusion {.field scope} for {.val {site_id}}: {.val {bad_scope}}",
      "i" = "Implemented: {.val {known_scopes}}"
    ))
  }

  t <- data[[time_col]]
  flag <- rep(FALSE, length(t))
  reason <- rep(NA_character_, length(t))
  for (i in seq_len(nrow(ex))) {
    start <- local_date_to_utc(ex$date_start[i], utc_offset_h)
    end <- local_date_to_utc(ex$date_end[i] + 1, utc_offset_h)
    inside <- !is.na(t) & t >= start & t < end
    reason[inside] <- ifelse(is.na(reason[inside]), ex$reason[i],
                             paste(reason[inside], ex$reason[i], sep = "; "))
    flag <- flag | inside
  }
  cli::cli_inform(
    "Exclusions {.val {site_id}}: {sum(flag)} of {length(flag)} record{?s} flagged ({nrow(ex)} period{?s})."
  )
  data$excl_flag <- flag
  data$excl_reason <- reason
  data
}

#' Drop flagged exclusion records
#'
#' @param data Data frame with logical column `excl_flag` (from [flag_exclusions()]).
#' @return `data` without the rows where `excl_flag` is `TRUE`; the number of
#'   removed rows is reported.
drop_excluded <- function(data) {
  if (!"excl_flag" %in% names(data)) {
    cli::cli_abort("Column {.field excl_flag} not found; run {.fn flag_exclusions} first.")
  }
  drop <- !is.na(data$excl_flag) & data$excl_flag
  cli::cli_inform("Dropped {sum(drop)} excluded record{?s} of {nrow(data)}.")
  data[!drop, , drop = FALSE]
}
