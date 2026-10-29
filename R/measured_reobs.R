# SIRTA ReOBS hourly mixing-layer height at Palaiseau (spec 10 Part A).
# Dataset: SIRTA 1D ReOBS (AERIS), Chiriaco et al. (2018, Chiriaco2018); MLH from the Lufft CHM15k with
# STRATfinder and from the Vaisala CL31 with CABAM v2019.3 (Kotthaus et al. 2018, Kotthaus2018; Kotthaus et al.
# 2020, Kotthaus2020). Hourly means of the native values: `time` = start of the hour [UTC], `time_bnds` = [t, t + 1 h).
# Licence: "Inspired by CC-BY" (file attribute `license`); acknowledgement text in the global attribute.

#' ReOBS MLH sources: variables and nominal sample counts
#'
#' Nominal number of native values per hour: CHM15k STRATfinder 1-min values (60), CL31 CABAM 15-min values (4)
#' (maximum of `_nb` in product 2026.08.1, spec 10 A2).
#'
#' @return Tibble: `instrument` (`chm15k`, `cl31`), `source_id`, `var` (NetCDF MLH variable), `nominal_nb` [-],
#'   `instrument_name`, `algorithm`, `native_resolution`, `citation_keys` [-].
reobs_mlh_sources <- function() {
  dplyr::tibble(
    instrument = c("chm15k", "cl31"),
    source_id = c("sirta_mlh_chm15k", "sirta_mlh_cl31"),
    var = c("mlh_chm15k", "mlh_cl31"),
    nominal_nb = c(60L, 4L),
    instrument_name = c("Lufft CHM15k automatic lidar ceilometer (SIRTA, Palaiseau)",
                        "Vaisala CL31 automatic lidar ceilometer (SIRTA, Palaiseau)"),
    algorithm = c("STRATfinder (Kotthaus et al. 2020), for high-SNR ceilometers",
                  "CABAM v2019.3 (Kotthaus and Grimmond 2018; Kotthaus et al. 2020), for low-SNR ceilometers"),
    native_resolution = c("1 min", "15 min"),
    citation_keys = c("Chiriaco2018; Kotthaus2020", "Chiriaco2018; Kotthaus2018; Kotthaus2020")
  )
}

#' Read a SIRTA ReOBS 1D NetCDF file
#'
#' Reads the hourly ReOBS file directly or from a zip archive: `ncdf4` cannot read inside a zip, so the single
#' `reobs_*.nc` member is extracted to a temporary directory that is removed on exit (never unzipped into the
#' repository). For every variable in `vars`, the value and, where present, the ancillary `<var>_std` and
#' `<var>_nb` are returned. Fill values (`NaN` for values, `-999` for counts) become `NA`. Aborts if `time` is
#' not "minutes since 1970-01-01", if `time_bnds` is not [time, time + 60 min), or if a variable is missing.
#'
#' @param path Path of the `.nc` file or of a `.zip` containing exactly one `reobs_*.nc` [character].
#' @param vars ReOBS variables to read [character]; heights in m above ground, `tcf_cl31` in oktas.
#' @return Tibble: `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC, from `time_bnds`], and per variable
#'   `<var>` [unit of the file], `<var>_std` [same unit], `<var>_nb` [integer, number of native values]. Attribute
#'   `reobs_meta`: list with `file`, `title`, `product_version`, `reobs_version`, `uuid`, `licence`,
#'   `acknowledgement`, `references`, `metadata_link`, `station_name`, `station_lat` [deg N], `station_lon`
#'   [deg E], `station_alt` [m a.s.l.].
read_reobs <- function(path, vars = c("mlh_chm15k", "mlh_cl31", "cbh_chm15k", "cbh_cl31", "tcf_cl31",
                                      "ablh_chm15k")) {
  if (!file.exists(path)) cli::cli_abort("ReOBS file {.path {path}} not found.")
  nc_path <- path
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    members <- utils::unzip(path, list = TRUE)$Name
    nc_member <- members[grepl("(^|/)reobs_[^/]*\\.nc$", members)]
    if (length(nc_member) != 1) {
      cli::cli_abort("{.path {path}}: expected exactly one reobs_*.nc member, found {length(nc_member)}.")
    }
    tmp <- tempfile("reobs_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    utils::unzip(path, files = nc_member, exdir = tmp, junkpaths = TRUE)
    nc_path <- file.path(tmp, basename(nc_member))
  }
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE, after = FALSE)

  t_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  if (!identical(t_units, "minutes since 1970-01-01")) {
    cli::cli_abort("{.path {path}}: unexpected time units {.val {t_units}} (expected minutes since 1970-01-01, UTC).")
  }
  tm <- as.numeric(ncdf4::ncvar_get(nc, "time"))
  tb <- ncdf4::ncvar_get(nc, "time_bnds")
  tb <- matrix(as.numeric(tb), nrow = 2)
  if (ncol(tb) != length(tm) || any(tb[1, ] != tm) || any(tb[2, ] - tb[1, ] != 60)) {
    cli::cli_abort("{.path {path}}: time_bnds is not [time, time + 60 min) on every row.")
  }
  to_utc <- function(min) as.POSIXct(min * 60, origin = "1970-01-01", tz = "UTC")
  out <- dplyr::tibble(timestamp_start_utc = to_utc(tb[1, ]), timestamp_end_utc = to_utc(tb[2, ]))

  miss <- setdiff(vars, names(nc$var))
  if (length(miss) > 0) cli::cli_abort("{.path {path}}: variable{?s} {.val {miss}} not found.")
  get <- function(v, count = FALSE) {
    x <- as.vector(ncdf4::ncvar_get(nc, v))
    if (count) {
      x <- as.integer(x)
      x[!is.na(x) & x == -999L] <- NA_integer_
    } else {
      x <- as.numeric(x)
      x[is.nan(x)] <- NA_real_
    }
    x
  }
  for (v in vars) {
    out[[v]] <- get(v)
    if (paste0(v, "_std") %in% names(nc$var)) out[[paste0(v, "_std")]] <- get(paste0(v, "_std"))
    if (paste0(v, "_nb") %in% names(nc$var)) out[[paste0(v, "_nb")]] <- get(paste0(v, "_nb"), count = TRUE)
  }

  ga <- ncdf4::ncatt_get(nc, 0)
  st <- function(v) if (v %in% names(nc$var)) ncdf4::ncvar_get(nc, v) else NA
  attr(out, "reobs_meta") <- list(
    file = basename(nc_path), title = ga$title, product_version = ga$product_version,
    reobs_version = ga$reobs_version, uuid = ga$uuid, licence = ga$license, acknowledgement = ga$acknowledgement,
    references = ga$references, metadata_link = ga$metadata_link,
    station_name = as.character(st("station_name")), station_lat = as.numeric(st("station_lat")),
    station_lon = as.numeric(st("station_lon")), station_alt = as.numeric(st("station_alt"))
  )
  out
}

#' Quality flag of an hourly ReOBS MLH from its sample count
#'
#' Spec 10 A3: `0` if `nb` >= 75 % of the nominal count, `1` if >= 50 %, otherwise `NA` (value not used).
#'
#' @param nb Number of native values in the hour [-].
#' @param nominal Nominal number per hour [-] (60 for CHM15k, 4 for CL31).
#' @return Integer `0`, `1` or `NA`.
reobs_qc <- function(nb, nominal) {
  frac <- nb / nominal
  out <- rep(NA_integer_, length(nb))
  out[!is.na(frac) & frac >= 0.5] <- 1L
  out[!is.na(frac) & frac >= 0.75] <- 0L
  out
}

#' Standardize a ReOBS MLH source (`sirta_mlh_chm15k` or `sirta_mlh_cl31`)
#'
#' Hourly interval means [start, start + 1 h), `method_class = aerosol`. `pblh_m` is the ReOBS hourly mean MLH
#' where [reobs_qc()] is 0 or 1, otherwise `NA`; `n_samples` = `_nb`; `uncertainty_m` = `_std` (standard
#' deviation of the native values within the hour, not a retrieval uncertainty). Rows from the first to the last
#' hour with a value or a count for this instrument. Exclusion periods of the ICOS site are applied and counted.
#'
#' @param r Output of [read_reobs()] (with `mlh_<instrument>`, `_std`, `_nb`, `cbh_<instrument>`, `tcf_cl31`,
#'   and `ablh_chm15k` for the CHM15k).
#' @param instrument `"chm15k"` or `"cl31"`.
#' @param site_id ICOS station paired with SIRTA [character].
#' @param exclusions Output of [read_exclusions()].
#' @param utc_offset_h ICOS site UTC offset [h] (for the exclusion periods).
#' @return List: `data` (standard table with `x_cbh_m` [m agl, cloud base of the same instrument],
#'   `x_tcf_oktas` [oktas, CL31 total cloud fraction], `x_coverage` [`_nb` / nominal, 0-1], and for the CHM15k
#'   `x_ablh_m` [m agl, aerosol-layer top, not a mixing height]), `meta` (sidecar list), `n_excluded` [-].
standardize_reobs_mlh <- function(r, instrument, site_id, exclusions, utc_offset_h) {
  src <- reobs_mlh_sources()
  src <- src[src$instrument == instrument, ]
  if (nrow(src) != 1) cli::cli_abort("Unknown ReOBS instrument {.val {instrument}} (chm15k or cl31).")
  v <- src$var
  need <- c(v, paste0(v, c("_std", "_nb")), paste0("cbh_", instrument), "tcf_cl31",
            if (instrument == "chm15k") "ablh_chm15k")
  miss <- setdiff(need, names(r))
  if (length(miss) > 0) cli::cli_abort("ReOBS table: column{?s} {.field {miss}} missing.")
  rm <- attr(r, "reobs_meta")

  has <- which(!is.na(r[[v]]) | !is.na(r[[paste0(v, "_nb")]]))
  if (length(has) == 0) cli::cli_abort("ReOBS {.val {v}}: no data.")
  r <- r[seq(min(has), max(has)), ]
  qc <- reobs_qc(r[[paste0(v, "_nb")]], src$nominal_nb)
  qc[is.na(r[[v]])] <- NA_integer_
  d <- dplyr::tibble(
    timestamp_start_utc = r$timestamp_start_utc, timestamp_end_utc = r$timestamp_end_utc,
    pblh_m = ifelse(is.na(qc), NA_real_, r[[v]]), source_id = src$source_id, qc_flag = qc,
    uncertainty_m = ifelse(is.na(qc), NA_real_, r[[paste0(v, "_std")]]),
    n_samples = r[[paste0(v, "_nb")]],
    x_cbh_m = r[[paste0("cbh_", instrument)]], x_tcf_oktas = r$tcf_cl31,
    x_coverage = r[[paste0(v, "_nb")]] / src$nominal_nb
  )
  if (instrument == "chm15k") d$x_ablh_m <- r$ablh_chm15k
  d <- flag_exclusions(d, site_id, exclusions, utc_offset_h, time_col = "timestamp_start_utc")
  ex <- drop_excluded_rows(d, src$source_id)
  d <- ex$data
  d$excl_flag <- NULL
  d$excl_reason <- NULL

  meta <- standard_meta(
    site_id, src$source_id, "aerosol",
    instrument = src$instrument_name,
    retrieval_algorithm = paste0(src$algorithm, "; mixing-layer height from attenuated backscatter, computed by ",
                                 "SIRTA/AERIS (not retrieved in this project)"),
    height_reference = "agl",
    instrument_altitude_m = rm$station_alt,
    time_reference = "interval_mean",
    original_resolution = paste0(src$native_resolution, " MLH values, hourly mean in ReOBS (cell_methods time: mean)"),
    original_qc = "ReOBS quality control at the native resolution (Chiriaco et al. 2018); hourly mean, std and count",
    qc_flag_meaning = list(`0` = paste0("n_samples >= 75 % of nominal (", src$nominal_nb, " per hour)"),
                           `1` = "n_samples >= 50 % and < 75 % of nominal",
                           `NA` = "no value, or n_samples < 50 % of nominal (pblh_m set to NA)"),
    processing_history = paste0(
      "SIRTA 1D ReOBS file ", rm$file, ", product_version ", rm$product_version, ", ReOBS ", rm$reobs_version,
      ", uuid ", rm$uuid, " (", rm$metadata_link, "); read by read_reobs() from the zip without permanent ",
      "extraction; uncertainty_m = within-hour std of the native values; rows in exclusion periods dropped (",
      ex$n_excluded, "); standardized by standardize_reobs_mlh() (spec 10)."
    ),
    licence = paste0(rm$licence, " (file attribute 'license'; ", rm$metadata_link, ")"),
    citation = paste0(src$citation_keys, " (docs/references/pblh_literature_2026-10.bib). Acknowledgement: ",
                      rm$acknowledgement)
  )
  list(data = d, meta = meta, n_excluded = ex$n_excluded)
}
