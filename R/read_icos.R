# ICOS L2 readers (FLUXES, FLUXNET FLUXMET HH, METEO).
# Input facts: CLAUDE.md "Input data". Rules: docs/specs/01_metadata_readers.md section 4.
# Timestamps in the files are local standard time (no DST); converted to UTC here.
# Units are converted at this boundary only (degC -> K, kPa -> Pa); columns in
# K get a unit suffix (_K, _PA). Declared units are sanity-checked on the column
# median before conversion (abort if implausible). No filtering, no gap-filling.

icos_na <- c("-9999", "-9999.0", "-9999.00", "NA", "")

icos_cols <- list(
  fluxes = c(
    "TIMESTAMP_START", "TIMESTAMP_END", "H", "H_UNCLEANED", "H_SSITC_TEST", "USTAR",
    "MO_LENGTH", "ZL", "T_SONIC", "WS", "WD", "U_SIGMA", "V_SIGMA", "W_SIGMA", "TAU",
    "TAU_SSITC_TEST", "PBLH", "FETCH_MAX", "FETCH_50", "FETCH_70", "FETCH_80",
    "FETCH_90", "FOOTPRINT_TA_CONTR", "FOOTPRINT_FLAG"
  ),
  fluxnet = c(
    "TIMESTAMP_START", "TIMESTAMP_END", "TA_F", "TA_F_QC", "PA_F", "PA_F_QC",
    "H_F_MDS", "H_F_MDS_QC", "H_CORR", "NIGHT", "SW_IN_POT", "WS_F", "P_F"
  ),
  meteo = c(
    "TIMESTAMP_START", "TIMESTAMP_END", "TA", "PA", "P", "SW_IN", "WS",
    paste0("TA_", 1:10)
  )
)

# Declared units of the raw columns handled at the reader boundary.
# T_SONIC is already in K in the L2 FLUXES files (checked on BE-Bra, spec 01):
# not converted, only renamed T_SONIC_K.
icos_degc_cols <- c("TA_F", "TA", paste0("TA_", 1:10))
icos_kelvin_cols <- c("T_SONIC")
icos_kpa_cols <- c("PA_F", "PA")

# Plausible range of the column MEDIAN for each declared unit (unit sanity check,
# not a data filter): degC [-40, 50], K [220, 330], kPa [50, 110].
icos_unit_median_ranges <- list(
  degC = c(-40, 50),
  K = c(220, 330),
  kPa = c(50, 110)
)

# Time columns never prefixed / never unit-converted.
icos_time_cols <- c("TIMESTAMP_START", "TIMESTAMP_END", "timestamp_start_utc", "timestamp_end_utc")

#' Locate the ICOS L2 zip files of a station
#'
#' Looks in `<icos_l2_dir>/<site_id>_output_L2I_2026/` for
#' `ICOSETC_<site>_FLUXES_INTERIM_L2.zip` (required),
#' `ICOSETC_<site>_FLUXNET_FLUXMET_HH_*.zip` and `ICOSETC_<site>_METEO_INTERIM_L2.zip`
#' (optional, warning if missing). Aborts if a pattern matches more than one file.
#'
#' @param site_id Station code, e.g. `"BE-Bra"` [character].
#' @param cfg Configuration list (from [get_cfg()]) with element `icos_l2_dir`.
#' @param dir_suffix Station folder suffix [character]; default `"_output_L2I_2026"`.
#' @return Named list of absolute paths: `fluxes`, `fluxnet`, `meteo` (`NULL` if missing).
find_icos_files <- function(site_id, cfg, dir_suffix = "_output_L2I_2026") {
  dir <- file.path(resolve_path(cfg$icos_l2_dir, site = site_id), paste0(site_id, dir_suffix))
  if (!dir.exists(dir)) {
    cli::cli_abort("ICOS L2 folder for {.val {site_id}} not found: {.path {dir}}")
  }
  site_re <- gsub("([.\\-])", "\\\\\\1", site_id)
  patterns <- c(
    fluxes = paste0("^ICOSETC_", site_re, "_FLUXES_INTERIM_L2\\.zip$"),
    fluxnet = paste0("^ICOSETC_", site_re, "_FLUXNET_FLUXMET_HH_.*\\.zip$"),
    meteo = paste0("^ICOSETC_", site_re, "_METEO_INTERIM_L2\\.zip$")
  )
  out <- list(fluxes = NULL, fluxnet = NULL, meteo = NULL)
  for (k in names(patterns)) {
    hits <- list.files(dir, pattern = patterns[[k]], full.names = TRUE)
    if (length(hits) > 1) {
      cli::cli_abort(c(
        "More than one {.val {k}} file for {.val {site_id}} in {.path {dir}}:",
        stats::setNames(basename(hits), rep("*", length(hits)))
      ))
    }
    if (length(hits) == 0) {
      if (k == "fluxes") {
        cli::cli_abort("Required FLUXES file missing for {.val {site_id}} in {.path {dir}}.")
      }
      cli::cli_warn("Optional {toupper(k)} file missing for {.val {site_id}} in {.path {dir}}.")
      next
    }
    out[k] <- list(hits)
  }
  out
}

# Parse "YYYYMMDDHHMM" local standard time -> POSIXct UTC (local clock - offset).
icos_time_to_utc <- function(x, utc_offset_h) {
  t <- as.POSIXct(x, format = "%Y%m%d%H%M", tz = "UTC")
  t - utc_offset_h * 3600
}

#' Read one ICOS L2 CSV (zipped or plain) with checks
#'
#' Internal workhorse for [read_icos_fluxes()], [read_icos_fluxnet()] and
#' [read_icos_meteo()]: reads only `cols` (absent columns -> warning), converts
#' timestamps to UTC, checks duplicates (abort), order (sort + warn), step and
#' duration != 30 min (warn), completes the regular 30-min grid, converts units.
#'
#' @param path Path to a zip holding one CSV, or to a CSV [character].
#' @param cols Columns to read [character].
#' @param utc_offset_h Site UTC offset [h, local standard time minus UTC].
#' @param label File label for messages [character].
#' @return Tibble on a complete 30-min grid with `timestamp_start_utc`,
#'   `timestamp_end_utc` [POSIXct UTC], original `TIMESTAMP_START`,
#'   `TIMESTAMP_END` [character, local standard time; `NA` on added rows],
#'   `.in_file` [logical, row present in the file], and the data columns
#'   (unit-converted columns renamed with `_K` / `_PA`). Attribute `read_info`:
#'   one-row tibble with file, period, row counts and missing columns.
read_icos_csv <- function(path, cols, utc_offset_h, label = basename(path)) {
  if (length(utc_offset_h) != 1 || is.na(utc_offset_h)) {
    cli::cli_abort("{.arg utc_offset_h} must be a single non-missing value [h].")
  }
  if (!file.exists(path)) cli::cli_abort("File not found: {.path {path}}")

  d <- readr::read_csv(
    path,
    col_select = dplyr::any_of(cols),
    col_types = readr::cols(
      .default = readr::col_double(),
      TIMESTAMP_START = readr::col_character(),
      TIMESTAMP_END = readr::col_character()
    ),
    na = icos_na,
    progress = FALSE,
    show_col_types = FALSE
  )
  if (nrow(d) == 0) cli::cli_abort("{label}: no data rows in {.path {path}}.")
  missing_cols <- setdiff(cols, names(d))
  if (any(c("TIMESTAMP_START", "TIMESTAMP_END") %in% missing_cols)) {
    cli::cli_abort("{label}: {.field TIMESTAMP_START}/{.field TIMESTAMP_END} not found.")
  }
  if (length(missing_cols) > 0) {
    cli::cli_warn("{label}: {length(missing_cols)} requested column{?s} absent: {.field {missing_cols}}")
  }

  # Local standard time -> UTC
  d$timestamp_start_utc <- icos_time_to_utc(d$TIMESTAMP_START, utc_offset_h)
  d$timestamp_end_utc <- icos_time_to_utc(d$TIMESTAMP_END, utc_offset_h)
  bad <- is.na(d$timestamp_start_utc)
  if (any(bad)) {
    cli::cli_abort("{label}: {sum(bad)} unparseable {.field TIMESTAMP_START} value{?s}, e.g. {.val {utils::head(d$TIMESTAMP_START[bad], 3)}}")
  }

  dup <- duplicated(d$timestamp_start_utc)
  if (any(dup)) {
    cli::cli_abort("{label}: {sum(dup)} duplicated timestamp{?s}, e.g. {.val {utils::head(d$TIMESTAMP_START[dup], 3)}}")
  }
  if (is.unsorted(d$timestamp_start_utc)) {
    cli::cli_warn("{label}: timestamps not in chronological order; sorted.")
    d <- d[order(d$timestamp_start_utc), ]
  }

  step <- 1800 # half-hourly, s
  dt <- diff(as.numeric(d$timestamp_start_utc))
  n_step <- sum(dt != step)
  if (n_step > 0) {
    cli::cli_warn("{label}: {n_step} step{?s} between consecutive records != 30 min.")
  }
  off_grid <- (as.numeric(d$timestamp_start_utc) - as.numeric(d$timestamp_start_utc[1])) %% step != 0
  if (any(off_grid)) {
    cli::cli_abort("{label}: {sum(off_grid)} timestamp{?s} not on the 30-min grid, e.g. {.val {utils::head(d$TIMESTAMP_START[off_grid], 3)}}")
  }
  dur <- as.numeric(d$timestamp_end_utc) - as.numeric(d$timestamp_start_utc)
  n_dur <- sum(is.na(dur) | dur != step)
  if (n_dur > 0) {
    cli::cli_warn("{label}: {n_dur} record{?s} with TIMESTAMP_END - TIMESTAMP_START != 30 min.")
  }

  check_icos_units(d, label)

  d$.in_file <- TRUE
  n_file <- nrow(d)
  d <- complete_halfhour_grid(d, label)

  d <- convert_icos_units(d)
  d <- dplyr::relocate(d, dplyr::all_of(c("timestamp_start_utc", "timestamp_end_utc")))

  attr(d, "read_info") <- dplyr::tibble(
    file = label,
    path = path,
    first_start_utc = min(d$timestamp_start_utc),
    last_start_utc = max(d$timestamp_start_utc),
    n_rows_file = n_file,
    n_rows_added = nrow(d) - n_file,
    missing_cols = paste(missing_cols, collapse = ", ")
  )
  d
}

# Put a table on a complete 30-min grid (first to last timestamp_start_utc).
# Added rows: NA data, timestamp_end_utc = start + 30 min, .in_file = FALSE.
complete_halfhour_grid <- function(d, label) {
  grid <- dplyr::tibble(timestamp_start_utc = seq(
    min(d$timestamp_start_utc), max(d$timestamp_start_utc), by = 1800
  ))
  out <- dplyr::left_join(grid, d, by = "timestamp_start_utc")
  added <- is.na(out$.in_file)
  out$.in_file[added] <- FALSE
  out$timestamp_end_utc[added] <- out$timestamp_start_utc[added] + 1800
  cli::cli_inform("{label}: {sum(added)} row{?s} added to complete the 30-min grid ({nrow(out)} rows).")
  out
}

# Unit sanity check on raw columns: abort if the median of a column is outside
# the plausible range of its declared unit (icos_unit_median_ranges).
check_icos_units <- function(d, label) {
  declared <- list(degC = icos_degc_cols, K = icos_kelvin_cols, kPa = icos_kpa_cols)
  for (u in names(declared)) {
    rng <- icos_unit_median_ranges[[u]]
    for (v in intersect(declared[[u]], names(d))) {
      if (all(is.na(d[[v]]))) next
      med <- stats::median(d[[v]], na.rm = TRUE)
      if (med < rng[1] || med > rng[2]) {
        cli::cli_abort(c(
          "{label}: unit check failed for {.field {v}} (declared {u}).",
          "x" = "Median = {signif(med, 6)}, expected within [{rng[1]}, {rng[2]}] {u}.",
          "i" = "The file may use a different unit; check before converting."
        ))
      }
    }
  }
  invisible(TRUE)
}

# Unit conversion at the reader boundary; columns renamed with unit suffix.
convert_icos_units <- function(d) {
  for (v in intersect(icos_degc_cols, names(d))) {
    d[[v]] <- d[[v]] + 273.15 # degC -> K (SI definition, 0 degC = 273.15 K)
    names(d)[names(d) == v] <- paste0(v, "_K")
  }
  for (v in intersect(icos_kelvin_cols, names(d))) {
    names(d)[names(d) == v] <- paste0(v, "_K") # already K: rename only
  }
  for (v in intersect(icos_kpa_cols, names(d))) {
    d[[v]] <- d[[v]] * 1000 # kPa -> Pa
    names(d)[names(d) == v] <- paste0(v, "_PA")
  }
  d
}

# Prefix data columns (not time columns, not .in_file).
prefix_cols <- function(d, prefix) {
  keep <- names(d) %in% c(icos_time_cols, ".in_file")
  names(d)[!keep] <- paste0(prefix, names(d)[!keep])
  d
}

#' Read an ICOS L2 FLUXES file
#'
#' @param path Path to `ICOSETC_<site>_FLUXES_INTERIM_L2.zip` [character].
#' @param utc_offset_h Site UTC offset [h, local standard time minus UTC].
#' @return Tibble on a complete 30-min grid (see [read_icos_csv()]). Columns
#'   (ICOS names and units unless suffixed): `timestamp_start_utc`,
#'   `timestamp_end_utc` [POSIXct UTC]; `TIMESTAMP_START`, `TIMESTAMP_END`
#'   [character, YYYYMMDDHHMM local standard time]; `H`, `H_UNCLEANED` sensible
#'   heat flux [W m-2]; `H_SSITC_TEST`, `TAU_SSITC_TEST` QC flags [-, 0/1/2];
#'   `USTAR` friction velocity [m s-1]; `MO_LENGTH` Obukhov length [m]; `ZL`
#'   stability parameter [-]; `T_SONIC_K` sonic temperature [K]; `WS` wind speed
#'   [m s-1]; `WD` wind direction [deg from N]; `U_SIGMA`, `V_SIGMA`, `W_SIGMA`
#'   standard deviations of wind components [m s-1]; `TAU` momentum flux
#'   [kg m-1 s-2]; `PBLH` operational PBL height [m]; `FETCH_MAX`, `FETCH_50`,
#'   `FETCH_70`, `FETCH_80`, `FETCH_90` footprint distances [m];
#'   `FOOTPRINT_TA_CONTR` target-area contribution [%]; `FOOTPRINT_FLAG` [-];
#'   `src_fluxes` [logical, row present in the file].
read_icos_fluxes <- function(path, utc_offset_h) {
  d <- read_icos_csv(path, icos_cols$fluxes, utc_offset_h, label = "FLUXES")
  info <- attr(d, "read_info")
  names(d)[names(d) == ".in_file"] <- "src_fluxes"
  attr(d, "read_info") <- info
  d
}

#' Read an ICOS L2 FLUXNET FLUXMET HH file
#'
#' @param path Path to `ICOSETC_<site>_FLUXNET_FLUXMET_HH_<years>_v*.zip` [character].
#' @param utc_offset_h Site UTC offset [h, local standard time minus UTC].
#' @return Tibble on a complete 30-min grid (see [read_icos_csv()]), data columns
#'   prefixed `fn_`: `fn_TA_F_K` air temperature, gap-filled [K]; `fn_TA_F_QC`
#'   [-]; `fn_PA_F_PA` air pressure, gap-filled [Pa]; `fn_PA_F_QC` [-];
#'   `fn_H_F_MDS` sensible heat flux, MDS gap-filled [W m-2]; `fn_H_F_MDS_QC`
#'   [-]; `fn_H_CORR` energy-balance corrected H [W m-2]; `fn_NIGHT` night flag
#'   [-, 1 = night]; `fn_SW_IN_POT` potential shortwave radiation [W m-2];
#'   `fn_WS_F` wind speed, gap-filled [m s-1]; `fn_P_F` precipitation, gap-filled
#'   [mm]; plus the time columns and `src_fluxnet` [logical].
read_icos_fluxnet <- function(path, utc_offset_h) {
  d <- read_icos_csv(path, icos_cols$fluxnet, utc_offset_h, label = "FLUXNET")
  info <- attr(d, "read_info")
  d <- prefix_cols(d, "fn_")
  names(d)[names(d) == ".in_file"] <- "src_fluxnet"
  attr(d, "read_info") <- info
  d
}

#' Read an ICOS L2 METEO file
#'
#' @param path Path to `ICOSETC_<site>_METEO_INTERIM_L2.zip` [character].
#' @param utc_offset_h Site UTC offset [h, local standard time minus UTC].
#' @return Tibble on a complete 30-min grid (see [read_icos_csv()]), data columns
#'   prefixed `met_`: `met_TA_K` air temperature [K]; `met_PA_PA` air pressure
#'   [Pa]; `met_P` precipitation [mm]; `met_SW_IN` incoming shortwave radiation
#'   [W m-2]; `met_WS` wind speed [m s-1]; `met_TA_1_K` ... `met_TA_10_K`
#'   temperature profile [K]; plus the time columns and `src_meteo` [logical].
read_icos_meteo <- function(path, utc_offset_h) {
  d <- read_icos_csv(path, icos_cols$meteo, utc_offset_h, label = "METEO")
  info <- attr(d, "read_info")
  d <- prefix_cols(d, "met_")
  names(d)[names(d) == ".in_file"] <- "src_meteo"
  attr(d, "read_info") <- info
  d
}

#' Read and merge all ICOS L2 data of a station
#'
#' Full outer join of FLUXES, FLUXNET and METEO on `timestamp_start_utc`, put on
#' a complete 30-min grid, with site heights ([get_site_heights()], evaluated at
#' `timestamp_start_utc`) and exclusion flags ([flag_exclusions()]).
#' `TIMESTAMP_START`/`TIMESTAMP_END` are coalesced from the files (priority
#' FLUXES, FLUXNET, METEO) and are `NA` on rows absent from all files.
#'
#' @param site_id Station code [character].
#' @param cfg Configuration list with `icos_l2_dir`.
#' @param sites Output of [read_sites()].
#' @param sites_tv Output of [read_sites_timevarying()].
#' @param exclusions Output of [read_exclusions()].
#' @return Tibble, one row per half-hour, sorted by time: `site_id`,
#'   `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC], `TIMESTAMP_START`,
#'   `TIMESTAMP_END` [character, local standard time], `src_fluxes`,
#'   `src_fluxnet`, `src_meteo` [logical, row present in that file], the data
#'   columns of the three readers, `zm_m` [m], `zm_source`, `canopy_height_m`
#'   [m], `hc_source`, `excl_flag` [logical], `excl_reason`. Attribute
#'   `read_info`: one row per file read (period, rows added, missing columns).
read_icos_site <- function(site_id, cfg, sites, sites_tv, exclusions) {
  meta <- sites[sites$site_id == site_id, ]
  if (nrow(meta) != 1) cli::cli_abort("Site {.val {site_id}} not found in site metadata.")
  off <- meta$utc_offset_h
  files <- find_icos_files(site_id, cfg)

  tabs <- list(fluxes = read_icos_fluxes(files$fluxes, off))
  if (!is.null(files$fluxnet)) tabs$fluxnet <- read_icos_fluxnet(files$fluxnet, off)
  if (!is.null(files$meteo)) tabs$meteo <- read_icos_meteo(files$meteo, off)
  info <- dplyr::bind_rows(lapply(tabs, attr, "read_info"))

  # Join: suffix time columns per table, coalesce afterwards.
  sfx <- c(fluxes = "fx", fluxnet = "fn", meteo = "met")
  tcols <- c("timestamp_end_utc", "TIMESTAMP_START", "TIMESTAMP_END")
  for (k in names(tabs)) {
    i <- names(tabs[[k]]) %in% tcols
    names(tabs[[k]])[i] <- paste0(names(tabs[[k]])[i], ".", sfx[[k]])
  }
  d <- Reduce(function(a, b) dplyr::full_join(a, b, by = "timestamp_start_utc"), tabs)
  for (v in tcols) {
    cand <- intersect(paste0(v, ".", sfx), names(d))
    d[[v]] <- do.call(dplyr::coalesce, unname(as.list(d[cand])))
    d[cand] <- NULL
  }
  for (s in c("src_fluxes", "src_fluxnet", "src_meteo")) {
    if (!s %in% names(d)) d[[s]] <- FALSE
  }
  d$.in_file <- TRUE
  d <- d[order(d$timestamp_start_utc), ]
  d <- complete_halfhour_grid(d, paste0(site_id, " merged"))
  d$.in_file <- NULL
  for (s in c("src_fluxes", "src_fluxnet", "src_meteo")) d[[s]][is.na(d[[s]])] <- FALSE

  d$site_id <- site_id
  d <- dplyr::relocate(d, dplyr::all_of(c(
    "site_id", "timestamp_start_utc", "timestamp_end_utc", "TIMESTAMP_START",
    "TIMESTAMP_END", "src_fluxes", "src_fluxnet", "src_meteo"
  )))

  h <- get_site_heights(site_id, d$timestamp_start_utc, sites, sites_tv)
  d <- dplyr::bind_cols(d, h[c("zm_m", "zm_source", "canopy_height_m", "hc_source")])
  d <- flag_exclusions(d, site_id, exclusions, off, time_col = "timestamp_start_utc")
  attr(d, "read_info") <- info
  d
}

#' Plausibility / coverage report for ICOS data (report only, never filters)
#'
#' Default ranges (spec 01, for reporting only; not physical filters):
#' `H` [-200, 800] W m-2; `USTAR` [0, 3] m s-1; `T_SONIC_K` [223, 323] K;
#' `fn_PA_F_PA` [50000, 110000] Pa; `PBLH` [10, 5000] m; `|MO_LENGTH|` > 0 m.
#' Rows of `ranges` with `NA` bounds give coverage only (`n_out_of_range` = `NA`).
#' By default, rows flagged by [flag_exclusions()] (`excl_flag` TRUE) are not
#' counted; their number is reported.
#'
#' @param data Output of [read_icos_site()] (or a reader).
#' @param include_excluded If `FALSE` (default), skip rows with `excl_flag`
#'   `TRUE` (no effect if `data` has no `excl_flag` column) [logical].
#' @param ranges Tibble with `variable`, `lower`, `upper` (units of the
#'   variable), `abs` [logical, compare `|x|`] and optional `strict_lower`
#'   [logical, require `x > lower`]; see [icos_default_ranges()]. Bounds are
#'   otherwise inclusive.
#' @return Tibble: `variable`, `n_total` [-, rows counted], `n_excluded` [-,
#'   rows flagged as excluded; part of `n_total` only if `include_excluded`], `n_valid` [-, non-missing],
#'   `pct_valid` [%], `n_out_of_range` [-, among valid], `min`, `max` [units of
#'   the variable, valid values].
icos_range_report <- function(data, ranges = icos_default_ranges(),
                              include_excluded = FALSE) {
  excl <- if ("excl_flag" %in% names(data)) !is.na(data$excl_flag) & data$excl_flag else rep(FALSE, nrow(data))
  n_excluded <- sum(excl)
  if (!include_excluded) {
    data <- data[!excl, , drop = FALSE]
    cli::cli_inform("Range report: {n_excluded} excluded row{?s} skipped.")
  } else {
    cli::cli_inform("Range report: {n_excluded} excluded row{?s} included.")
  }
  rows <- lapply(seq_len(nrow(ranges)), function(i) {
    v <- ranges$variable[i]
    n_total <- nrow(data)
    if (!v %in% names(data)) {
      return(dplyr::tibble(variable = v, n_total = n_total, n_excluded = n_excluded, n_valid = 0L,
                           pct_valid = 0, n_out_of_range = NA_integer_,
                           min = NA_real_, max = NA_real_))
    }
    x <- data[[v]]
    x <- x[!is.na(x)]
    y <- if (isTRUE(ranges$abs[i])) abs(x) else x
    lo <- ranges$lower[i]
    hi <- ranges$upper[i]
    n_oor <- if (is.na(lo) && is.na(hi)) {
      NA_integer_
    } else if (isTRUE(ranges$strict_lower[i])) {
      sum(y <= lo | (!is.na(hi) & y > hi))
    } else {
      sum((!is.na(lo) & y < lo) | (!is.na(hi) & y > hi))
    }
    dplyr::tibble(
      variable = v, n_total = n_total, n_excluded = n_excluded, n_valid = length(x),
      pct_valid = if (n_total > 0) 100 * length(x) / n_total else NA_real_,
      n_out_of_range = as.integer(n_oor),
      min = if (length(x)) min(x) else NA_real_,
      max = if (length(x)) max(x) else NA_real_
    )
  })
  dplyr::bind_rows(rows)
}

#' Default plausibility ranges for [icos_range_report()]
#'
#' @return Tibble `variable`, `lower`, `upper` (units of the variable), `abs`
#'   [logical], `strict_lower` [logical, `x > lower` required instead of `>=`].
icos_default_ranges <- function() {
  dplyr::tibble(
    variable = c("H", "USTAR", "T_SONIC_K", "fn_PA_F_PA", "PBLH", "MO_LENGTH"),
    lower = c(-200, 0, 223, 50000, 10, 0),
    upper = c(800, 3, 323, 110000, 5000, NA),
    abs = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    strict_lower = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  )
}
