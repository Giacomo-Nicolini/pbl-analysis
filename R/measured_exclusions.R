# Exclusion periods of a measured PBLH source only (config/measured_exclusions.csv; spec 13, G 2026-11-06).
# Unlike config/exclusion_periods.csv (which removes the ICOS site's own records), these periods remove only the values
# of one measured source, e.g. months with a dirty ceilometer window. Applied by the standard-file builders
# (standardize_meteoswiss_mlh(), standardize_reobs_mlh(), standardize_cloudnet_mlh()); counts always reported.

#' Read the measured-source exclusion periods
#'
#' Reads `config/measured_exclusions.csv`: `site_id` (ICOS station of the pairing), `source_id` (standard source id,
#' e.g. `meteoswiss_cabam_mlh`), `date_start`, `date_end` (UTC calendar dates, `YYYY-MM-DD`; `date_end` inclusive:
#' the whole UTC day is excluded; measured sources are in UTC), `reason`.
#'
#' @param file Path to `measured_exclusions.csv` [character].
#' @return Tibble: `site_id`, `source_id` [-], `date_start`, `date_end` [Date, UTC, inclusive], `reason` [-].
read_measured_exclusions <- function(file) {
  ex <- readr::read_csv(
    file, na = c("NA", ""),
    col_types = readr::cols(site_id = readr::col_character(), source_id = readr::col_character(),
                            date_start = readr::col_date(format = "%Y-%m-%d"),
                            date_end = readr::col_date(format = "%Y-%m-%d"), reason = readr::col_character()),
    progress = FALSE
  )
  need <- c("site_id", "source_id", "date_start", "date_end", "reason")
  miss <- setdiff(need, names(ex))
  if (length(miss)) cli::cli_abort("{.path {file}}: missing column{?s} {.field {miss}}.")
  bad <- which(!stats::complete.cases(ex[need]))
  if (length(bad)) cli::cli_abort("{.path {file}}: every field is required (rows {bad}).")
  bad <- which(ex$date_end < ex$date_start)
  if (length(bad)) cli::cli_abort("{.path {file}}: {.field date_end} before {.field date_start} (rows {bad}).")
  bad <- which(!ex$source_id %in% style_ids("source"))
  if (length(bad)) cli::cli_abort("{.path {file}}: unknown {.field source_id} {.val {unique(ex$source_id[bad])}} (rows {bad}).")
  ex[need]
}

#' Drop the values of a measured source inside its exclusion periods
#'
#' Interval per period: `date_start` 00:00 UTC (inclusive) to `date_end + 1 day` 00:00 UTC (exclusive), compared with
#' `timestamp_start_utc`. Rows are dropped (as for the site exclusion periods); the count is always reported, also
#' when zero.
#'
#' @param d Standard table of one source (`timestamp_start_utc` [POSIXct UTC]).
#' @param site_id ICOS station [character]; @param source_id Standard source id [character].
#' @param mex Output of [read_measured_exclusions()]; `NULL` = no table (reported).
#' @return List: `data` (`d` without the excluded rows), `n_excluded` [-], `periods` [character, the applied periods
#'   with reasons, `"none"` if none], `history` [character, sentence for the sidecar `processing_history`].
apply_measured_exclusions <- function(d, site_id, source_id, mex) {
  if (is.null(mex)) {
    cli::cli_inform("{source_id}: no measured-source exclusion table given.")
    return(list(data = d, n_excluded = 0L, periods = "none",
                history = "no measured-source exclusion table given"))
  }
  ex <- mex[mex$site_id == site_id & mex$source_id == source_id, ]
  t <- d$timestamp_start_utc
  drop <- rep(FALSE, length(t))
  for (i in seq_len(nrow(ex))) {
    t0 <- as.POSIXct(ex$date_start[i], tz = "UTC")
    t1 <- as.POSIXct(ex$date_end[i] + 1, tz = "UTC")
    drop <- drop | (!is.na(t) & t >= t0 & t < t1)
  }
  periods <- if (nrow(ex)) paste0(ex$date_start, " to ", ex$date_end, " (", ex$reason, ")", collapse = "; ") else "none"
  cli::cli_inform("{source_id}: {sum(drop)} of {length(drop)} row{?s} in {nrow(ex)} measured-source exclusion period{?s} dropped.")
  list(data = d[!drop, , drop = FALSE], n_excluded = sum(drop), periods = periods,
       history = paste0("rows in measured-source exclusion periods (config/measured_exclusions.csv) dropped (",
                        sum(drop), "; periods: ", periods, ")"))
}
