# Convert PBLH sources to the standard format (R/standard_format.R).
# Pure functions: each returns list(data = standard table, meta = sidecar, ...);
# writing is done by the caller. Rules: docs/specs/04_era5_standard_compare.md
# section 5 and amendments 2-5.

# Drop rows with excl_flag == TRUE; report the count (also zero).
drop_excluded_rows <- function(d, label) {
  ex <- !is.na(d$excl_flag) & d$excl_flag
  cli::cli_inform("{label}: {sum(ex)} of {nrow(d)} row{?s} in exclusion periods dropped.")
  list(data = d[!ex, ], n_excluded = sum(ex))
}

# Coerce ISO 8601 character timestamps (e.g. from a csv) to POSIXct UTC.
as_utc_time <- function(x) {
  if (inherits(x, "POSIXct")) return(lubridate::with_tz(x, "UTC"))
  lubridate::ymd_hms(x, tz = "UTC")
}

#' Standardize PBLHmod output
#'
#' `qc_flag`: 0 no flag; 1 forcing interpolated (`flag_ustar_interp`,
#' `flag_L_interp` or `flag_H_interp`); 2 h held by the model
#' (`flag_sign_inconsistent`); 3 below the measurement height
#' (`flag_below_zm`); the highest applicable value wins. Rows in exclusion
#' periods (`excl_flag`) are dropped and counted.
#'
#' @param out Output of [run_pblh_model()] or its derived csv (same columns).
#' @param label Run label, e.g. `"default"` [character]; `source_id =
#'   "pblhmod_<label>"`.
#' @param opts Optional [pblh_options()] of the run, written to the sidecar
#'   (`NULL`: taken from attribute `pblh_options` of `out` if present).
#' @return List: `data` (standard table: `pblh_m` [m above ground], `qc_flag`,
#'   `x_regime`, `x_method`, `x_h_source`), `meta` (sidecar), `n_excluded`.
standardize_pblhmod <- function(out, label, opts = NULL) {
  source_id <- paste0("pblhmod_", label)
  site_id <- unique(out$site_id)
  if (length(site_id) != 1) cli::cli_abort("{.arg out} must hold exactly one site.")
  if (is.null(opts)) opts <- attr(out, "pblh_options")
  ex <- drop_excluded_rows(out, source_id)
  d <- ex$data
  fl <- function(k) if (k %in% names(d)) !is.na(d[[k]]) & as.logical(d[[k]]) else rep(FALSE, nrow(d))
  qc <- integer(nrow(d))
  qc[fl("flag_ustar_interp") | fl("flag_L_interp") | fl("flag_H_interp")] <- 1L
  qc[fl("flag_sign_inconsistent")] <- 2L
  qc[fl("flag_below_zm")] <- 3L

  data <- dplyr::tibble(
    timestamp_start_utc = as_utc_time(d$timestamp_start_utc),
    timestamp_end_utc = as_utc_time(d$timestamp_end_utc),
    pblh_m = as.numeric(d$pblh_m),
    source_id = source_id,
    qc_flag = qc,
    x_regime = as.character(d$regime),
    x_method = as.character(d$method),
    x_h_source = as.character(d$H_source)
  )
  opt_txt <- if (is.null(opts)) {
    "options: unknown"
  } else {
    o <- opts[setdiff(names(opts), "label")]
    paste0("options: ", paste(names(o), vapply(o, function(v) paste(format(v), collapse = ","), ""),
                              sep = " = ", collapse = "; "))
  }
  meta <- standard_meta(
    site_id, source_id, "parameterisation",
    instrument = "eddy-covariance system (ICOS L2: FLUXES, FLUXNET, METEO)",
    retrieval_algorithm = paste(
      "PBLHmod v2 (spec 02): Kljun et al. 2015 (doi:10.5194/gmd-8-3695-2015) App. B diagnostics",
      "(B1 stable, B2/B3 neutral) and Batchvarova & Gryning 1991 slab model (B5) in the convective window"
    ),
    height_reference = "agl",
    time_reference = "interval_mean",
    original_resolution = "30 min",
    original_qc = "logical flags of run_pblh_model() (R/pblh_model.R)",
    qc_flag_meaning = list(
      `0` = "no flag",
      `1` = "forcing interpolated (flag_ustar_interp, flag_L_interp or flag_H_interp)",
      `2` = "h held by the model (flag_sign_inconsistent)",
      `3` = "below the measurement height (flag_below_zm)"
    ),
    processing_history = paste0(
      "run_pblh_model() (scripts/02_run_pblh_model.R, run '", label, "'; ", opt_txt,
      "); raw output, no smoothing, clipping or filling; rows in exclusion periods dropped (",
      ex$n_excluded, "); standardized by standardize_pblhmod() (spec 04)."
    )
  )
  list(data = data, meta = meta, n_excluded = ex$n_excluded)
}

# Short provenance text from a parsed ERA5 header (parse_era5_header()).
era5_provenance <- function(header) {
  if (is.null(header)) return("ERA5 extraction details unknown")
  g <- function(k) if (is.null(header[[k]])) "unknown" else format(header[[k]])
  paste0(
    "ERA5 ", g("dataset"), " / ", g("variable"), " extracted by C. Trotta (Oct 2026) with ", g("script"),
    " (method ", g("method"), ", ARCO zarr ", g("geochunked_url"), ", downloaded ", g("download_utc"),
    "): nearest 0.25 deg grid point to the tower (", g("latitude_requested"), " N, ", g("longitude_requested"),
    " E; cell centre ", g("latitude_grid"), " N, ", g("longitude_grid"), " E) and its 8 neighbours (3x3 cells)"
  )
}

era5_licence <- "Licence to use Copernicus Products"
era5_citation <- paste(
  "Hersbach et al. 2020 (doi:10.1002/qj.3803); ERA5 hourly data on single levels",
  "(doi:10.24381/cds.adbb2d47) and on pressure levels (doi:10.24381/cds.bd0915c6)"
)

#' Standardize ERA5 boundary-layer height
#'
#' Instantaneous values (`timestamp_start_utc == timestamp_end_utc` =
#' `timestamp_utc`), `qc_flag = NA` (no QC applied; amendment 3), 3x3 spread
#' as extra columns. Rows in the site exclusion periods (converted with the
#' site UTC offset, [flag_exclusions()]) are dropped and counted.
#'
#' @param era5_center BLH of the tower cell: output of
#'   [read_era5_single_level()] (`cells = "center"`, or `"all"`: the tower
#'   cell is selected).
#' @param era5_bilin Output of [era5_bilinear()] on `blh_m`.
#' @param spread Output of [era5_spread()] on `blh_m`.
#' @param variant `"center"` (`source_id = "era5_blh_center"`) or
#'   `"bilinear"` (`"era5_blh_bilinear"`).
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h Site UTC offset [h].
#' @param header Optional parsed header of the tower-cell BLH file (column
#'   `header` of [find_era5_files()]) for `processing_history`.
#' @return List: `data` (standard table: `pblh_m` [m above model surface],
#'   `x_blh_sd9_m`, `x_blh_range9_m` [m], `x_blh_cv9` [-], `x_n9` [-]),
#'   `meta`, `n_excluded`.
standardize_era5_blh <- function(era5_center, era5_bilin, spread, variant = c("center", "bilinear"),
                                 exclusions, utc_offset_h, header = NULL) {
  variant <- match.arg(variant)
  if ("cell_row" %in% names(era5_center)) {
    era5_center <- era5_center[era5_center$cell_row == 0 & era5_center$cell_col == 0, ]
  }
  site_id <- unique(era5_center$site_id)
  if (length(site_id) != 1) cli::cli_abort("{.arg era5_center} must hold exactly one site.")
  source_id <- paste0("era5_blh_", variant)

  base <- if (variant == "center") {
    dplyr::tibble(timestamp_utc = era5_center$timestamp_utc, pblh_m = era5_center$blh_m)
  } else {
    dplyr::tibble(timestamp_utc = era5_bilin$timestamp_utc, pblh_m = era5_bilin$blh_m_bilin)
  }
  sp <- dplyr::tibble(
    timestamp_utc = spread$timestamp_utc, x_blh_sd9_m = spread$sd9, x_blh_range9_m = spread$range9,
    x_blh_cv9 = spread$cv9, x_n9 = as.integer(spread$n9)
  )
  d <- dplyr::left_join(base, sp, by = "timestamp_utc")
  d <- d[order(d$timestamp_utc), ]
  d$timestamp_start_utc <- d$timestamp_utc
  d <- flag_exclusions(d, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")
  ex <- drop_excluded_rows(d, source_id)
  d <- ex$data

  data <- dplyr::tibble(
    timestamp_start_utc = d$timestamp_utc,
    timestamp_end_utc = d$timestamp_utc,
    pblh_m = d$pblh_m,
    source_id = source_id,
    qc_flag = NA_integer_,
    x_blh_sd9_m = d$x_blh_sd9_m,
    x_blh_range9_m = d$x_blh_range9_m,
    x_blh_cv9 = d$x_blh_cv9,
    x_n9 = d$x_n9
  )
  meta <- standard_meta(
    site_id, source_id, "model",
    instrument = "none (ECMWF ERA5 reanalysis)",
    retrieval_algorithm = paste(
      "IFS Cy41r2 diagnostic: bulk Richardson number = 0.25 from the lowest model level, no parcel excess",
      "(ECMWF 2016, Part IV Sect. 3.10.1, Eq. 3.90)"
    ),
    height_reference = "above model surface",
    time_reference = "instantaneous",
    original_resolution = if (variant == "center") {
      "hourly, 0.25 deg grid cell containing the tower"
    } else {
      "hourly, 0.25 deg grid, bilinear interpolation of the 4 surrounding cell centres to the tower position"
    },
    qc_flag_meaning = list(`NA` = "no QC applied (ERA5 reanalysis)"),
    processing_history = paste0(
      era5_provenance(header), ". Read with read_era5_single_level()",
      if (variant == "bilinear") "; interpolated with era5_bilinear()" else "",
      "; x_ columns = spread over the 3x3 cells (era5_spread()); rows in exclusion periods dropped (",
      ex$n_excluded, "); standardized by standardize_era5_blh() (spec 04)."
    ),
    licence = era5_licence,
    citation = era5_citation
  )
  list(data = data, meta = meta, n_excluded = ex$n_excluded)
}

#' Detect run-mean fill values in a PBLH series
#'
#' ETC `proc_PBL_height()` fills gaps with the mean of the day or night values
#' of a processing run; the run boundaries are not known in general. Rule
#' (spec 04 section 5, amendment 7): a value, rounded to 0.01 m, other than the
#' floor value, is a fill value if it occurs in at least `min_n` rows in total
#' **and** in at least `min_frac` of the non-missing rows of at least one
#' calendar year (the year stands in for the processing run). Only the rows
#' of the qualifying years are fill rows; rows with the same value in other
#' years are counted in `n_outside` (coincidences).
#'
#' @param x PBLH [m].
#' @param year Calendar year of each row (local standard time, from
#'   `TIMESTAMP_START`) [integer, same length as `x`].
#' @param min_n Minimum number of occurrences in the whole record [-]; default 50.
#' @param min_frac Minimum fraction of the non-missing rows of a year [-];
#'   default 0.005.
#' @param floor_value Floor value excluded from the search [m]; default 50.
#' @return Tibble sorted by `n_total` (decreasing), one row per fill value:
#'   `value` [m, rounded to 0.01], `n_total` [-], `years` [list of integer,
#'   qualifying years], `share` [list of numeric, fraction of the non-missing
#'   rows of each qualifying year], `n_flagged` [-, rows in qualifying years],
#'   `n_outside` [-, rows in other years].
detect_fill_values <- function(x, year, min_n = 50, min_frac = 0.005, floor_value = 50) {
  if (length(year) != length(x)) cli::cli_abort("{.arg year} must have the length of {.arg x}.")
  empty <- dplyr::tibble(value = numeric(0), n_total = integer(0), years = list(), share = list(),
                         n_flagged = integer(0), n_outside = integer(0))
  ok <- !is.na(x)
  r <- sprintf("%.2f", x[ok]) # rounded to 0.01 m, compared as text
  y <- year[ok]
  if (length(r) == 0) return(empty)
  tab <- table(r)
  cand <- names(tab)[tab >= min_n & abs(as.numeric(names(tab)) - floor_value) > 1e-9]
  if (length(cand) == 0) return(empty)
  n_year <- table(y) # non-missing rows per year
  rows <- lapply(cand, function(v) {
    ny <- table(y[r == v])
    share <- as.numeric(ny) / as.numeric(n_year[names(ny)])
    q <- share >= min_frac
    if (!any(q)) return(NULL)
    dplyr::tibble(
      value = as.numeric(v), n_total = as.integer(tab[[v]]),
      years = list(as.integer(names(ny)[q])), share = list(share[q]),
      n_flagged = as.integer(sum(ny[q])), n_outside = as.integer(sum(ny[!q]))
    )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) return(empty)
  out[order(-out$n_total), ]
}

# TRUE for rows holding a fill value in one of its qualifying years.
is_fill_row <- function(x, year, fills) {
  out <- rep(FALSE, length(x))
  if (nrow(fills) == 0) return(out)
  r <- sprintf("%.2f", x)
  for (i in seq_len(nrow(fills))) {
    out <- out | (!is.na(x) & r == sprintf("%.2f", fills$value[i]) & year %in% fills$years[[i]])
  }
  out
}

# One line per fill value: total rows, qualifying years with share, rows outside.
format_fill_values <- function(fills) {
  if (nrow(fills) == 0) return(character(0))
  vapply(seq_len(nrow(fills)), function(i) {
    yrs <- paste0(fills$years[[i]], " (", sprintf("%.2f", 100 * fills$share[[i]]), " %)", collapse = ", ")
    sprintf("%.2f m: %d rows in total; years %s: %d rows flagged; %d rows outside those years",
            fills$value[i], fills$n_total[i], yrs, fills$n_flagged[i], fills$n_outside[i])
  }, character(1))
}

#' Standardize the operational ICOS PBLH
#'
#' The `PBLH` column of FLUXES L2 (ETC processing: legacy `calc_PBL_height()`
#' + `proc_PBL_height()`, `docs/legacy_review.md`). `qc_flag`: 0 as computed;
#' 1 floored (`pblh_m == 50`); 2 run-mean fill (value and year from
#' [detect_fill_values()], amendment 7); `NA` where `pblh_m` is missing. Rows
#' outside the FLUXES file (if `src_fluxes` is present) and rows in exclusion
#' periods are dropped (the latter counted).
#'
#' @param d Output of [read_icos_site()] (`PBLH` [m], `TIMESTAMP_START`
#'   [character, `YYYYMMDDHHMM` local standard time], `excl_flag`).
#' @param ... Passed to [detect_fill_values()].
#' @return List: `data` (standard table, `pblh_m` [m]), `meta`, `n_excluded`,
#'   `fill_values` (output of [detect_fill_values()]), `n_floor` [-].
standardize_icos_oper <- function(d, ...) {
  site_id <- unique(d$site_id)
  if (length(site_id) != 1) cli::cli_abort("{.arg d} must hold exactly one site.")
  for (k in c("PBLH", "TIMESTAMP_START")) {
    if (!k %in% names(d)) cli::cli_abort("Column {.field {k}} not found.")
  }
  source_id <- "icos_oper"
  if ("src_fluxes" %in% names(d)) d <- d[!is.na(d$src_fluxes) & d$src_fluxes, ]
  ex <- drop_excluded_rows(d, source_id)
  d <- ex$data
  x <- as.numeric(d$PBLH)
  year <- as.integer(substr(d$TIMESTAMP_START, 1, 4)) # local standard time
  fills <- detect_fill_values(x, year, ...)
  floor_value <- list(...)$floor_value
  if (is.null(floor_value)) floor_value <- 50
  is_floor <- !is.na(x) & abs(x - floor_value) < 1e-9
  is_fill <- is_fill_row(x, year, fills)
  qc <- ifelse(is.na(x), NA_integer_, 0L)
  qc[is_floor] <- 1L
  qc[is_fill] <- 2L
  fill_txt <- format_fill_values(fills)
  cli::cli_inform(c(
    "{source_id}: {sum(is_floor)} floored row{?s} ({floor_value} m); {nrow(fills)} fill value{?s} detected:",
    stats::setNames(fill_txt, rep("*", length(fill_txt)))
  ))

  data <- dplyr::tibble(
    timestamp_start_utc = d$timestamp_start_utc,
    timestamp_end_utc = d$timestamp_end_utc,
    pblh_m = x,
    source_id = source_id,
    qc_flag = as.integer(qc)
  )
  meta <- standard_meta(
    site_id, source_id, "parameterisation",
    instrument = "eddy-covariance system (ICOS ETC operational processing)",
    retrieval_algorithm = "legacy calc_PBL_height() (ETC implementation of Kljun et al. 2015 App. B; docs/legacy_review.md)",
    height_reference = "agl",
    time_reference = "interval_mean",
    original_resolution = "30 min",
    original_qc = "none in the file; qc_flag reconstructed from the floor value and the detected fill values",
    qc_flag_meaning = list(
      `0` = "as computed",
      `1` = paste0("floored at ", floor_value, " m (proc_PBL_height)"),
      `2` = "run-mean fill (proc_PBL_height; value and years detected by detect_fill_values(), spec 04 amendment 7)"
    ),
    processing_history = paste0(
      "ICOS ETC operational processing: calc_PBL_height() then proc_PBL_height() (negative values -> NA; ",
      "floor at 50 m; remaining NA filled with the mean of all day (hour of TIMESTAMP_END 06-18, local ",
      "standard time) or night values of the processing run). Read from ICOSETC_", site_id,
      "_FLUXES_INTERIM_L2 (PBLH) with read_icos_site(); fill values detected: ",
      if (nrow(fills) == 0) "none" else paste(fill_txt, collapse = "; "),
      "; rows in exclusion periods dropped (", ex$n_excluded, "); standardized by standardize_icos_oper() (spec 04)."
    ),
    licence = "CC BY 4.0"
  )
  list(data = data, meta = meta, n_excluded = ex$n_excluded, fill_values = fills, n_floor = sum(is_floor))
}
