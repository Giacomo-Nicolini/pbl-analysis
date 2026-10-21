# Standard PBLH format: one data table (csv.gz) + YAML sidecar per site x source.
# Definition: docs/project_plan.md ("Standard PBLH format");
# implementation rules: docs/specs/04_era5_standard_compare.md section 4 and
# docs/specs/standard_format.md.

standard_required_cols <- c("timestamp_start_utc", "timestamp_end_utc", "pblh_m", "source_id")
standard_recommended_cols <- c("qc_flag", "uncertainty_m", "n_samples")

#' Method classes of the standard PBLH format
#'
#' @return Character vector: `turbulence` (Doppler lidar TKE / sigma_w, sodar),
#'   `aerosol` (ceilometer, backscatter lidar), `thermodynamic` (radiosonde,
#'   Ri_b, microwave radiometer), `model` (ERA5), `parameterisation` (PBLHmod),
#'   `unknown` (`docs/project_plan.md`).
standard_method_classes <- function() {
  c("turbulence", "aerosol", "thermodynamic", "model", "parameterisation", "unknown")
}

#' Optional sidecar fields counted in `metadata_completeness`
#'
#' The automatic fields `created_utc` and `code_version` are excluded (spec 04,
#' amendment 4).
#'
#' @return Character vector of 12 field names.
standard_optional_fields <- function() {
  c(
    "instrument", "retrieval_algorithm", "height_reference", "instrument_altitude_m",
    "time_reference", "original_resolution", "original_qc", "qc_flag_meaning",
    "processing_history", "contact", "licence", "citation"
  )
}

#' Build a standard sidecar
#'
#' @param site_id Station code [character].
#' @param source_id Source identifier, e.g. `"era5_blh_center"` [character].
#' @param method_class One of [standard_method_classes()].
#' @param ... Optional fields of [standard_optional_fields()]; missing ones are
#'   set to `"unknown"`. Unknown names abort.
#' @return Named list (sidecar without the automatic fields).
standard_meta <- function(site_id, source_id, method_class, ...) {
  opt <- list(...)
  bad <- setdiff(names(opt), standard_optional_fields())
  if (length(bad) > 0) {
    cli::cli_abort("Unknown sidecar field{?s}: {.field {bad}}.")
  }
  meta <- list(site_id = site_id, source_id = source_id, method_class = method_class)
  for (f in standard_optional_fields()) {
    meta[[f]] <- if (is.null(opt[[f]])) "unknown" else opt[[f]]
  }
  meta
}

# TRUE if a sidecar value is unknown: NULL, NA, empty or "unknown".
is_unknown_field <- function(x) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  if (is.list(x)) return(FALSE)
  all(is.na(x)) || (is.character(x) && length(x) == 1 && tolower(trimws(x)) %in% c("unknown", ""))
}

#' Metadata completeness of a sidecar
#'
#' Fraction of the optional fields ([standard_optional_fields()], 12 fields;
#' `created_utc` and `code_version` excluded) that are not `unknown`.
#'
#' @param meta Sidecar list.
#' @return Fraction [0-1].
metadata_completeness <- function(meta) {
  f <- standard_optional_fields()
  mean(!vapply(f, function(k) is_unknown_field(meta[[k]]), logical(1)))
}

#' Validate a standard PBLH table and its sidecar
#'
#' Checks: required columns and types (`timestamp_*_utc` POSIXct, `pblh_m`
#' numeric, `source_id` character equal to the sidecar); recommended columns
#' numeric; other columns prefixed `x_`; no missing timestamps; start <= end;
#' `time_reference: instantaneous` implies start == end, `interval_mean`
#' start < end; no duplicated start times; `pblh_m` finite or `NA`; sidecar
#' `site_id`, `source_id` known and `method_class` valid. All problems are
#' collected and reported in one error.
#'
#' @param df Standard table.
#' @param meta Sidecar list.
#' @return `TRUE` invisibly; aborts listing all problems otherwise.
validate_standard_pblh <- function(df, meta) {
  p <- character(0)
  miss <- setdiff(standard_required_cols, names(df))
  if (length(miss) > 0) p <- c(p, paste0("missing required column(s): ", paste(miss, collapse = ", ")))
  extra <- setdiff(names(df), c(standard_required_cols, standard_recommended_cols))
  extra <- extra[!startsWith(extra, "x_")]
  if (length(extra) > 0) p <- c(p, paste0("extra column(s) without 'x_' prefix: ", paste(extra, collapse = ", ")))

  for (k in c("site_id", "source_id", "method_class")) {
    if (is_unknown_field(meta[[k]])) p <- c(p, paste0("sidecar field '", k, "' missing or unknown"))
  }
  if (!is.null(meta$method_class) && !meta$method_class %in% standard_method_classes()) {
    p <- c(p, paste0("invalid method_class '", meta$method_class, "'"))
  }

  ts_ok <- TRUE
  for (k in intersect(c("timestamp_start_utc", "timestamp_end_utc"), names(df))) {
    if (!inherits(df[[k]], "POSIXct")) {
      p <- c(p, paste0(k, " is not POSIXct"))
      ts_ok <- FALSE
    } else if (anyNA(df[[k]])) {
      p <- c(p, paste0(k, ": ", sum(is.na(df[[k]])), " missing value(s)"))
      ts_ok <- FALSE
    }
  }
  if (ts_ok && all(c("timestamp_start_utc", "timestamp_end_utc") %in% names(df))) {
    s <- df$timestamp_start_utc
    e <- df$timestamp_end_utc
    if (any(s > e)) p <- c(p, paste0(sum(s > e), " row(s) with start > end"))
    tr <- meta$time_reference
    if (identical(tr, "instantaneous") && any(s != e)) {
      p <- c(p, paste0("time_reference instantaneous but ", sum(s != e), " row(s) with start != end"))
    }
    if (identical(tr, "interval_mean") && any(s == e)) {
      p <- c(p, paste0("time_reference interval_mean but ", sum(s == e), " row(s) with start == end"))
    }
    if (anyDuplicated(s)) p <- c(p, paste0(sum(duplicated(s)), " duplicated start time(s)"))
  }
  if ("pblh_m" %in% names(df)) {
    if (!is.numeric(df$pblh_m)) {
      p <- c(p, "pblh_m is not numeric")
    } else if (any(is.infinite(df$pblh_m))) {
      p <- c(p, paste0("pblh_m: ", sum(is.infinite(df$pblh_m)), " infinite value(s)"))
    }
  }
  if ("source_id" %in% names(df)) {
    if (!is.character(df$source_id)) {
      p <- c(p, "source_id is not character")
    } else if (!is.null(meta$source_id) && any(df$source_id != meta$source_id)) {
      p <- c(p, paste0("source_id in the table differs from the sidecar ('", meta$source_id, "')"))
    }
  }
  for (k in intersect(standard_recommended_cols, names(df))) {
    if (!is.numeric(df[[k]]) && !all(is.na(df[[k]]))) p <- c(p, paste0(k, " is not numeric"))
  }

  if (length(p) > 0) {
    cli::cli_abort(c(
      "Standard PBLH table {.val {meta$site_id}} / {.val {meta$source_id}} is not valid:",
      stats::setNames(p, rep("x", length(p)))
    ))
  }
  invisible(TRUE)
}

# Order columns: required, recommended (added as NA if absent), x_ columns.
standard_order_cols <- function(df) {
  if (!"qc_flag" %in% names(df)) df$qc_flag <- NA_integer_
  if (!"uncertainty_m" %in% names(df)) df$uncertainty_m <- NA_real_
  if (!"n_samples" %in% names(df)) df$n_samples <- NA_integer_
  x <- sort(names(df)[startsWith(names(df), "x_")])
  df[c(standard_required_cols, standard_recommended_cols, x)]
}

# Git version of the working tree ("<short hash>[-dirty]"), or "unknown".
code_version <- function() {
  v <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(here::here()), "describe", "--always", "--dirty", "--abbrev=7"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0)
  )
  if (length(v) == 1 && nzchar(v) && is.null(attr(v, "status"))) v else "unknown"
}

fmt_iso_utc <- function(t) format(t, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

#' Paths of a standard file pair
#'
#' @param dir Directory of the station, `<standard_dir>/<site>` [character].
#' @param site_id,source_id Station code and source identifier [character].
#' @return Named character: `data` (`<site>__<source>.csv.gz`), `meta` (`.yml`).
standard_paths <- function(dir, site_id, source_id) {
  base <- file.path(dir, paste0(site_id, "__", source_id))
  c(data = paste0(base, ".csv.gz"), meta = paste0(base, ".yml"))
}

#' Write a standard PBLH table and its sidecar
#'
#' Validates with [validate_standard_pblh()], then writes
#' `<dir>/<site>__<source_id>.csv.gz` (times as ISO 8601 UTC
#' `YYYY-MM-DDTHH:MM:SSZ`, missing values `NA`) and the YAML sidecar with
#' `metadata_completeness`, `created_utc` and `code_version` added.
#'
#' @param df Standard table.
#' @param meta Sidecar list (e.g. from [standard_meta()]).
#' @param dir Output directory, created if needed [character].
#' @return Named character of the two paths, invisibly.
write_standard_pblh <- function(df, meta, dir) {
  validate_standard_pblh(df, meta)
  df <- standard_order_cols(df)
  meta$metadata_completeness <- round(metadata_completeness(meta), 4)
  meta$created_utc <- fmt_iso_utc(Sys.time())
  meta$code_version <- code_version()
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- standard_paths(dir, meta$site_id, meta$source_id)
  out <- df
  out$timestamp_start_utc <- fmt_iso_utc(out$timestamp_start_utc)
  out$timestamp_end_utc <- fmt_iso_utc(out$timestamp_end_utc)
  readr::write_csv(out, paths[["data"]], na = "NA")
  yaml::write_yaml(meta, paths[["meta"]])
  invisible(paths)
}

#' Read a standard PBLH table and its sidecar
#'
#' @param path Path of the `.csv.gz` data file; the sidecar is the `.yml` with
#'   the same base name.
#' @return Tibble (timestamps POSIXct UTC, `pblh_m` [m], `qc_flag` integer,
#'   `uncertainty_m` [m], `n_samples` integer, `x_` columns as guessed) with
#'   the sidecar list as attribute `meta`. Validated.
read_standard_pblh <- function(path) {
  meta_path <- sub("\\.csv(\\.gz)?$", ".yml", path)
  if (!file.exists(meta_path)) {
    cli::cli_abort("Sidecar {.file {meta_path}} not found.")
  }
  meta <- yaml::read_yaml(meta_path)
  df <- readr::read_csv(
    path,
    col_types = readr::cols(
      timestamp_start_utc = readr::col_character(),
      timestamp_end_utc = readr::col_character(),
      pblh_m = readr::col_double(),
      source_id = readr::col_character(),
      qc_flag = readr::col_integer(),
      uncertainty_m = readr::col_double(),
      n_samples = readr::col_integer(),
      .default = readr::col_guess()
    ),
    na = "NA",
    progress = FALSE
  )
  for (k in c("timestamp_start_utc", "timestamp_end_utc")) {
    df[[k]] <- as.POSIXct(df[[k]], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  }
  validate_standard_pblh(df, meta)
  attr(df, "meta") <- meta
  df
}
