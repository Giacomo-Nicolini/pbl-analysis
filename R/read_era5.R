# ERA5 readers for the C. Trotta delivery (Oct 2026): hourly single levels,
# 6-hourly pressure levels and static fields, 0.25 deg, tower cell + 8
# neighbours. Rules: docs/specs/04_era5_standard_compare.md section 2;
# file facts: CLAUDE.md ("Input data", ERA5).

era5_grid_res_deg <- 0.25 # ERA5 regular lat-lon grid of this delivery [deg]

#' ERA5 variable registry
#'
#' Single source of truth for the ERA5 variables used in this project: file
#' name, output column, units and time convention. Units and conventions from
#' the ECMWF parameter database (https://codes.ecmwf.int/grib/param-db/:
#' blh 159, cbh 228023, sshf 146, t 130, z 129, lsm 172, sdor 160) and the ERA5
#' documentation (Hersbach et al. 2020, doi:10.1002/qj.3803; ERA5 data
#' documentation, ECMWF Confluence, "parameter listings" and "accumulations").
#'
#' @return Tibble, one row per variable: `variable` (name in the file),
#'   `family` (`"single_level"`, `"pressure_level"`, `"static"`), `column`
#'   (output column), `unit_in`, `unit_out`, `conversion`, `time_reference`,
#'   `header_units` (value of the header key `units`; static files only, `NA`
#'   where the header carries no units).
era5_variables <- function() {
  dplyr::tribble(
    ~variable, ~family, ~column, ~unit_in, ~unit_out, ~conversion, ~time_reference, ~header_units,
    # blh: IFS diagnostic, bulk Richardson number 0.25 (ECMWF2016, Part IV, Sect. 3.10.1, Eq. 3.90)
    "boundary_layer_height", "single_level", "blh_m", "m", "m", "none",
    "instantaneous at timestamp_utc", NA,
    # cbh: missing when there is no cloud (ECMWF parameter 228023)
    "cloud_base_height", "single_level", "cbh_m", "m", "m", "none; NA = no cloud",
    "instantaneous at timestamp_utc", NA,
    # sshf: accumulated over the hour ending at the timestamp, ECMWF flux sign
    # (positive downward); ERA5 documentation, "accumulations"
    "surface_sensible_heat_flux", "single_level", "H_era5_Wm2",
    "J m-2 (hourly accumulation, positive downward)", "W m-2 (positive upward)", "-x / 3600",
    "interval (timestamp_utc - 1 h, timestamp_utc]", NA,
    "temperature", "pressure_level", "T_K", "K", "K", "none", "instantaneous at timestamp_utc", NA,
    # geopotential height zg = z / g0 (WMO); differs from geometric height by < 0.5 % below 10 km
    "geopotential", "pressure_level", "z_m2s2", "m2 s-2", "m2 s-2; zg_m [m]", "zg_m = z / g",
    "instantaneous at timestamp_utc", NA,
    "z", "static", "orog_m", "m2 s-2", "m", "orog_m = z / g", "static", "m**2 s**-2",
    "lsm", "static", "lsm", "0-1", "0-1", "none", "static", "(0 - 1)",
    "sdor", "static", "sdor_m", "m", "m", "none", "static", "m"
  )
}

era5_end_marker <- "# ==== FINE METADATA (dati sotto) ===="

#' Parse the metadata header of an ERA5 file
#'
#' Reads the `# key: value` lines between the first line and the end marker
#' `# ==== FINE METADATA (dati sotto) ====` (split at the first `": "`,
#' trimmed). Keys starting with `latitude_`, `longitude_`, `grid_`, `pixel_`,
#' and `n_rows`, `distance_km` are converted to numeric when the value starts
#' with a number (the leading signed number is taken, e.g.
#' `"+0 (+ = north, - = south)"` -> 0); otherwise they stay character (e.g.
#' `longitude_conversion`, `grid_scanning`). `levels` becomes an integer
#' vector (`NA` for `"n/a"`).
#'
#' @param path ERA5 csv file [character].
#' @param max_lines Lines searched for the end marker [integer].
#' @return Named list of header values, plus `n_header_lines` [integer]:
#'   lines to skip before the csv header.
parse_era5_header <- function(path, max_lines = 200L) {
  lines <- readLines(path, n = max_lines, warn = FALSE, encoding = "UTF-8")
  end <- which(trimws(lines) == era5_end_marker)
  if (length(end) == 0) {
    cli::cli_abort("{.file {basename(path)}}: ERA5 header end marker {.val {era5_end_marker}} not found in the first {max_lines} lines.")
  }
  end <- end[1]
  body <- lines[seq_len(end - 1)][-1]
  body <- body[grepl("^#\\s*[^:]+: ", body)]
  body <- sub("^#\\s*", "", body)
  pos <- regexpr(": ", body, fixed = TRUE)
  keys <- trimws(substr(body, 1, pos - 1))
  vals <- trimws(substr(body, pos + 2, nchar(body)))
  out <- stats::setNames(as.list(vals), keys)

  num_key <- grepl("^(latitude_|longitude_|grid_|pixel_)", keys) | keys %in% c("n_rows", "distance_km")
  num_re <- "^[+-]?[0-9]+(\\.[0-9]*)?([eE][+-]?[0-9]+)?"
  for (k in keys[num_key]) {
    v <- out[[k]]
    if (grepl(num_re, v)) out[[k]] <- as.numeric(regmatches(v, regexpr(num_re, v)))
  }
  if (!is.null(out$levels)) {
    lv <- out$levels
    out$levels <- if (tolower(lv) %in% c("n/a", "na", "")) {
      NA_integer_
    } else {
      as.integer(trimws(strsplit(lv, ",", fixed = TRUE)[[1]]))
    }
  }
  out$n_header_lines <- end
  out
}

# Regex of the file names of one family; groups: variable, first time, last
# time, north, south, east, west.
era5_file_regex <- function(site_id, family) {
  site_re <- gsub("([.\\-])", "\\\\\\1", site_id)
  offs <- "__north_([+-]\\d+)__south_([+-]\\d+)__east_([+-]\\d+)__west_([+-]\\d+)\\.csv$"
  switch(family,
    single_level = paste0("^", site_re, "__reanalysis_era5_single_level__([a-z0-9_]+?)__(\\d{8}T\\d{4})__(\\d{8}T\\d{4})", offs),
    pressure_level = paste0("^", site_re, "__reanalysis_era5_pressure_levels-pl__([a-z0-9_]+?)__(\\d{8}T\\d{4})__(\\d{8}T\\d{4})", offs),
    static = paste0("^", site_re, "__single-levels__(z|lsm|sdor)__(\\d{12})__(\\d{12})", offs),
    cli::cli_abort("Unknown ERA5 family {.val {family}}.")
  )
}

# Parse ERA5 time strings ("YYYY-MM-DD HH:MM[:SS]") to POSIXct UTC; abort on failure.
parse_era5_time <- function(x, label) {
  t <- lubridate::ymd_hms(x, tz = "UTC", truncated = 1, quiet = TRUE)
  if (anyNA(t)) {
    cli::cli_abort("{label}: {sum(is.na(t))} unparseable time value{?s}, e.g. {.val {utils::head(x[is.na(t)], 3)}}.")
  }
  t
}

#' Locate ERA5 files of a station and assign their grid cell
#'
#' Lists the files of one family, keeps those matching the naming pattern of
#' the family and a variable of [era5_variables()] (other files are ignored and
#' listed in a message) and reads their headers.
#'
#' The cell position comes from the header, never from the file name: the
#' tower cell is the file whose (`latitude_grid`, `longitude_grid`) is closest
#' to (`latitude_requested`, `longitude_requested`) and must contain the tower
#' (within half a grid step); `cell_row = round((grid_lat - tower_grid_lat) /
#' 0.25)` (+1 = north), `cell_col = round((grid_lon - tower_grid_lon) / 0.25)`
#' (+1 = east). The name is a cross-check only: single/pressure-level files use
#' direction flags (`north_+1` -> row +1, `south_+1` -> row -1, `east_+1` ->
#' col +1, `west_+1` -> col -1), static files signed offsets in `north_` (row)
#' and `east_` (col) with `south_+0__west_+0`; a disagreement aborts. Per
#' variable the 9 cells of a 3x3 block are expected: missing cells warn,
#' duplicated cells or cells outside the block abort.
#'
#' @param site_id Station code [character].
#' @param cfg Configuration list with `era5_sl_dir`, `era5_pl_dir`,
#'   `era5_static_dir` (may contain `{site}`, `{var}`).
#' @param family `"single_level"`, `"pressure_level"` or `"static"`.
#' @param variable Optional variable names [character]; `NULL` = all variables
#'   of the family in [era5_variables()].
#' @param sites Optional output of [read_sites()]; if given, the requested
#'   coordinates in the headers are checked with [check_era5_site_coords()].
#' @return Tibble, one row per file: `site_id`, `family`, `variable`, `path`,
#'   `file`, `grid_lat`, `grid_lon` cell centre [deg], `lat_requested`,
#'   `lon_requested` tower [deg], `grid_dist_km` tower to cell centre [km],
#'   `cell_row`, `cell_col` [integer, cells; +1 = north / east],
#'   `time_first`, `time_last` [POSIXct UTC, from the name],
#'   `header_n_rows` [integer], `header_first_time`, `header_last_time`
#'   [POSIXct UTC], `n_header_lines` [integer], `header` [list of the parsed
#'   header].
find_era5_files <- function(site_id, cfg, family = c("single_level", "pressure_level", "static"),
                            variable = NULL, sites = NULL) {
  family <- match.arg(family)
  reg <- era5_variables()
  known <- reg$variable[reg$family == family]
  if (is.null(variable)) variable <- known
  bad <- setdiff(variable, known)
  if (length(bad) > 0) {
    cli::cli_abort("ERA5 {family}: unknown variable{?s} {.val {bad}}; known: {.val {known}}.")
  }

  dirs <- switch(family,
    single_level = resolve_path(cfg$era5_sl_dir, site = site_id),
    pressure_level = resolve_path(cfg$era5_pl_dir, site = site_id),
    static = unique(vapply(variable, function(v) resolve_path(cfg$era5_static_dir, site = site_id, var = v),
                           character(1)))
  )
  all_files <- unlist(lapply(dirs, list.files, full.names = TRUE), use.names = FALSE)
  all_files <- unique(all_files[!dir.exists(all_files)])
  m <- stringr::str_match(basename(all_files), era5_file_regex(site_id, family))
  ok <- !is.na(m[, 1]) & m[, 2] %in% known
  # Static folder holds all stations locally: list only this station's ignored files.
  own <- startsWith(basename(all_files), paste0(site_id, "__"))
  ign <- !ok & (own | family != "static")
  if (any(ign)) {
    cli::cli_inform(c(
      "ERA5 {.val {site_id}} {family}: ignoring {sum(ign)} file{?s} not matching the naming pattern or the registry:",
      stats::setNames(basename(all_files[ign]), rep("*", sum(ign)))
    ))
  }
  keep <- ok & m[, 2] %in% variable
  m <- m[keep, , drop = FALSE]
  paths <- all_files[keep]
  if (length(paths) == 0) {
    cli::cli_abort("ERA5 {.val {site_id}} {family}: no files for {.val {variable}} in {.path {dirs}}.")
  }

  hdr <- lapply(paths, parse_era5_header)
  hnum <- function(key) {
    vapply(seq_along(hdr), function(i) {
      v <- hdr[[i]][[key]]
      if (is.null(v) || !is.numeric(v)) {
        cli::cli_abort("{.file {basename(paths[i])}}: header key {.field {key}} missing or not numeric.")
      }
      v
    }, numeric(1))
  }
  hchr <- function(key) vapply(hdr, function(h) if (is.null(h[[key]])) NA_character_ else h[[key]], character(1))
  name_fmt <- if (family == "static") "%Y%m%d%H%M" else "%Y%m%dT%H%M"
  ftime <- function(x) as.POSIXct(x, format = name_fmt, tz = "UTC")

  f <- dplyr::tibble(
    site_id = site_id,
    family = family,
    variable = m[, 2],
    path = paths,
    file = basename(paths),
    grid_lat = hnum("latitude_grid"),
    grid_lon = hnum("longitude_grid"),
    lat_requested = hnum("latitude_requested"),
    lon_requested = hnum("longitude_requested"),
    grid_res = hnum("grid_resolution_deg"),
    time_first = ftime(m[, 3]),
    time_last = ftime(m[, 4]),
    name_n = as.integer(m[, 5]), name_s = as.integer(m[, 6]),
    name_e = as.integer(m[, 7]), name_w = as.integer(m[, 8]),
    header_n_rows = as.integer(hnum("n_rows")),
    header_first_time = parse_era5_time(hchr("first_time"), "ERA5 header first_time"),
    header_last_time = parse_era5_time(hchr("last_time"), "ERA5 header last_time"),
    n_header_lines = vapply(hdr, function(h) h$n_header_lines, integer(1)),
    header = hdr
  )

  # Header consistency: station and variable.
  pid <- hchr("point_id")
  if (any(pid != site_id, na.rm = TRUE) || anyNA(pid)) {
    cli::cli_abort("ERA5 {.val {site_id}}: header {.field point_id} differs from the station in {.file {f$file[is.na(pid) | pid != site_id]}}.")
  }
  hvar <- hchr("variable")
  if (any(is.na(hvar) | hvar != f$variable)) {
    i <- which(is.na(hvar) | hvar != f$variable)[1]
    cli::cli_abort("{.file {f$file[i]}}: header variable {.val {hvar[i]}} differs from the name ({.val {f$variable[i]}}).")
  }
  if (any(abs(f$grid_res - era5_grid_res_deg) > 1e-9)) {
    cli::cli_abort("ERA5 {.val {site_id}}: {.field grid_resolution_deg} must be {era5_grid_res_deg}; found {.val {unique(f$grid_res)}}.")
  }
  if (length(unique(round(f$lat_requested, 6))) > 1 || length(unique(round(f$lon_requested, 6))) > 1) {
    cli::cli_abort("ERA5 {.val {site_id}}: requested coordinates differ between files.")
  }

  # Cell position from the header, per variable.
  f$grid_dist_km <- haversine_km(f$lat_requested, f$lon_requested, f$grid_lat, f$grid_lon)
  f$cell_row <- NA_integer_
  f$cell_col <- NA_integer_
  for (v in unique(f$variable)) {
    i <- which(f$variable == v)
    t0 <- i[which.min(f$grid_dist_km[i])]
    half <- era5_grid_res_deg / 2 + 1e-6
    if (abs(f$grid_lat[t0] - f$lat_requested[t0]) > half || abs(f$grid_lon[t0] - f$lon_requested[t0]) > half) {
      cli::cli_abort("ERA5 {.val {site_id}} {v}: no file for the cell containing the tower (closest cell centre {f$grid_lat[t0]}, {f$grid_lon[t0]}).")
    }
    f$cell_row[i] <- as.integer(round((f$grid_lat[i] - f$grid_lat[t0]) / era5_grid_res_deg))
    f$cell_col[i] <- as.integer(round((f$grid_lon[i] - f$grid_lon[t0]) / era5_grid_res_deg))
  }

  # Cross-check with the name.
  if (family == "static") {
    bad_conv <- f$name_s != 0 | f$name_w != 0
    if (any(bad_conv)) {
      cli::cli_abort("ERA5 static: expected {.code south_+0__west_+0} in {.file {f$file[bad_conv]}}.")
    }
    name_row <- f$name_n
    name_col <- f$name_e
  } else {
    flags <- c(f$name_n, f$name_s, f$name_e, f$name_w)
    if (any(!flags %in% c(0L, 1L))) {
      cli::cli_abort("ERA5 {family}: direction flags must be +0 or +1.")
    }
    name_row <- f$name_n - f$name_s
    name_col <- f$name_e - f$name_w
  }
  dis <- name_row != f$cell_row | name_col != f$cell_col
  if (any(dis)) {
    i <- which(dis)[1]
    cli::cli_abort(c(
      "ERA5 {.file {f$file[i]}}: cell position in the name and in the header disagree.",
      "x" = "Name: row {name_row[i]}, col {name_col[i]}; header (grid {f$grid_lat[i]}, {f$grid_lon[i]}): row {f$cell_row[i]}, col {f$cell_col[i]}."
    ))
  }

  # 3x3 block per variable.
  out_blk <- abs(f$cell_row) > 1 | abs(f$cell_col) > 1
  if (any(out_blk)) {
    cli::cli_abort("ERA5 {.val {site_id}}: cells outside the 3x3 block in {.file {f$file[out_blk]}}.")
  }
  dup <- duplicated(f[c("variable", "cell_row", "cell_col")])
  if (any(dup)) {
    cli::cli_abort("ERA5 {.val {site_id}} {family}: duplicated cell in {.file {f$file[dup]}}.")
  }
  for (v in unique(f$variable)) {
    have <- paste(f$cell_row[f$variable == v], f$cell_col[f$variable == v])
    all9 <- expand.grid(r = -1:1, c = -1:1)
    miss <- all9[!paste(all9$r, all9$c) %in% have, ]
    if (nrow(miss) > 0) {
      cli::cli_warn("ERA5 {.val {site_id}} {v}: {nrow(miss)} of 9 cells missing (row, col): {.val {paste0('(', miss$r, ', ', miss$c, ')')}}.")
    }
  }

  if (!is.null(sites)) check_era5_site_coords(f, sites)

  f <- f[order(f$variable, f$cell_row, f$cell_col), ]
  f[c(
    "site_id", "family", "variable", "path", "file", "grid_lat", "grid_lon", "lat_requested",
    "lon_requested", "grid_dist_km", "cell_row", "cell_col", "time_first", "time_last",
    "header_n_rows", "header_first_time", "header_last_time", "n_header_lines", "header"
  )]
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

#' Check the tower coordinates used for the ERA5 extraction
#'
#' Warns when `lat_requested` / `lon_requested` in the ERA5 headers differ from
#' the site metadata by more than `tol_deg` (sign of an extraction with
#' outdated coordinates). Values are never modified.
#'
#' @param files Output of [find_era5_files()] (one station).
#' @param sites Output of [read_sites()].
#' @param tol_deg Tolerance [decimal degrees]; default 0.001.
#' @return `TRUE` invisibly if within tolerance, `FALSE` (with a warning) otherwise.
check_era5_site_coords <- function(files, sites, tol_deg = 0.001) {
  site_id <- unique(files$site_id)
  if (length(site_id) != 1) {
    cli::cli_abort("ERA5 coordinate check: {.arg files} must hold exactly one site.")
  }
  meta <- sites[sites$site_id == site_id, ]
  if (nrow(meta) != 1) {
    cli::cli_abort("ERA5 coordinate check: site {.val {site_id}} not found in site metadata.")
  }
  era_lat <- unique(files$lat_requested)
  era_lon <- unique(files$lon_requested)
  dlat <- max(abs(era_lat - meta$lat))
  dlon <- max(abs(era_lon - meta$lon))
  if (dlat > tol_deg || dlon > tol_deg) {
    cli::cli_warn(c(
      "ERA5 {.val {site_id}}: requested coordinates differ from site metadata by more than {tol_deg} deg.",
      "x" = "ERA5 headers: lat {.val {era_lat}}, lon {.val {era_lon}}; site metadata: lat {meta$lat}, lon {meta$lon}.",
      "i" = "ERA5 may have been extracted for outdated coordinates."
    ))
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# Read the data block of one ERA5 file (row of find_era5_files()). Checks the
# number of rows and first/last time against the header. Returns a tibble with
# `timestamp_utc` [POSIXct UTC] and the numeric data columns.
read_era5_data <- function(f, time_col) {
  d <- readr::read_csv(
    f$path,
    skip = f$n_header_lines,
    col_types = do.call(readr::cols, c(
      stats::setNames(list(readr::col_character()), time_col),
      list(.default = readr::col_double())
    )),
    na = c("", "NA", "NaN", "nan"),
    progress = FALSE
  )
  if (!time_col %in% names(d)) {
    cli::cli_abort("{.file {f$file}}: time column {.field {time_col}} not found.")
  }
  t <- parse_era5_time(d[[time_col]], f$file)
  if (nrow(d) != f$header_n_rows) {
    cli::cli_abort("{.file {f$file}}: {nrow(d)} rows read, header {.field n_rows} = {f$header_n_rows}.")
  }
  if (nrow(d) > 0 && (min(t) != f$header_first_time || max(t) != f$header_last_time)) {
    cli::cli_abort(c(
      "{.file {f$file}}: time range differs from the header.",
      "x" = "Data {format(min(t))} to {format(max(t))}; header {format(f$header_first_time)} to {format(f$header_last_time)}."
    ))
  }
  d[[time_col]] <- NULL
  dplyr::bind_cols(dplyr::tibble(timestamp_utc = t), d)
}

# Duplicated times per cell abort; gaps in the regular series warn.
check_era5_times <- function(d, step_s, label) {
  key <- paste(d$cell_row, d$cell_col)
  dup <- duplicated(data.frame(key, d$timestamp_utc))
  if (any(dup)) {
    cli::cli_abort("{label}: {sum(dup)} duplicated timestamp{?s} within a cell, e.g. {.val {format(utils::head(d$timestamp_utc[dup], 3))}}.")
  }
  for (k in unique(key)) {
    dt <- diff(sort(as.numeric(d$timestamp_utc[key == k])))
    gaps <- dt != step_s
    if (any(gaps)) {
      cli::cli_warn("{label}, cell (row col) {k}: {sum(gaps)} gap{?s} ({sum(dt[gaps] / step_s - 1)} missing step{?s}).")
    }
  }
  invisible(TRUE)
}

#' Read ERA5 hourly single-level files
#'
#' Reads one variable (all cells or the tower cell) and converts it at the
#' reader boundary as in [era5_variables()]: `surface_sensible_heat_flux`
#' (hourly accumulation [J m-2], positive downward) becomes `H_era5_Wm2 =
#' -sshf / 3600` [W m-2, positive upward] over the interval
#' (`timestamp_utc - 1 h`, `timestamp_utc`] (ERA5 documentation,
#' "accumulations"). Checks (values are never modified): rows and first/last
#' time vs header (abort), duplicated times (abort), hourly gaps (warn),
#' median `blh_m` in [20, 3000] m (abort), diurnal sign of `H_era5_Wm2` with
#' [check_era5_h_sign()] (abort).
#'
#' @param files Output of [find_era5_files()] (family `"single_level"`).
#' @param variable `"boundary_layer_height"`, `"cloud_base_height"` or
#'   `"surface_sensible_heat_flux"`.
#' @param cells `"all"` or `"center"` (tower cell only).
#' @return Long tibble sorted by cell and time: `site_id`, `cell_row`,
#'   `cell_col` [integer], `grid_lat`, `grid_lon` [deg], `grid_dist_km` [km],
#'   `timestamp_utc` [POSIXct UTC], and `blh_m` [m above model surface,
#'   instantaneous] or `cbh_m` [m above surface, instantaneous; `NA` = no
#'   cloud] or `H_era5_Wm2` [W m-2, positive upward, mean over the hour] with
#'   `timestamp_start_utc`, `timestamp_end_utc` [POSIXct UTC].
read_era5_single_level <- function(files, variable, cells = c("all", "center")) {
  cells <- match.arg(cells)
  reg <- era5_variables()
  reg <- reg[reg$family == "single_level" & reg$variable == variable, ]
  if (nrow(reg) != 1) {
    cli::cli_abort("Unknown ERA5 single-level variable {.val {variable}}.")
  }
  f <- files[files$family == "single_level" & files$variable == variable, ]
  if (cells == "center") f <- f[f$cell_row == 0 & f$cell_col == 0, ]
  if (nrow(f) == 0) {
    cli::cli_abort("No ERA5 files for {.val {variable}} (cells = {.val {cells}}).")
  }

  d <- dplyr::bind_rows(lapply(seq_len(nrow(f)), function(i) {
    x <- read_era5_data(f[i, ], "time")
    if (!variable %in% names(x)) {
      cli::cli_abort("{.file {f$file[i]}}: column {.field {variable}} not found.")
    }
    dplyr::tibble(
      site_id = f$site_id[i], cell_row = f$cell_row[i], cell_col = f$cell_col[i],
      grid_lat = f$grid_lat[i], grid_lon = f$grid_lon[i], grid_dist_km = f$grid_dist_km[i],
      timestamp_utc = x$timestamp_utc, value = x[[variable]]
    )
  }))
  d <- d[order(d$cell_row, d$cell_col, d$timestamp_utc), ]
  label <- paste("ERA5", f$site_id[1], variable)
  check_era5_times(d, 3600, label)

  n_na <- sum(is.na(d$value))
  if (variable == "cloud_base_height") {
    cli::cli_inform("{label}: {n_na} of {nrow(d)} values missing ({round(100 * n_na / max(1, nrow(d)), 1)} %, no cloud).")
  } else if (n_na > 0) {
    cli::cli_warn("{label}: {n_na} missing value{?s} (kept as is).")
  }

  if (variable == "boundary_layer_height") {
    med <- stats::median(d$value, na.rm = TRUE)
    if (!is.finite(med) || med < 20 || med > 3000) {
      cli::cli_abort("{label}: median {round(med, 1)} outside [20, 3000] m; check the unit.")
    }
  }
  if (variable == "surface_sensible_heat_flux") {
    d$value <- -d$value / 3600 # J m-2 per hour, positive down -> W m-2, positive up
    d$timestamp_start_utc <- d$timestamp_utc - 3600
    d$timestamp_end_utc <- d$timestamp_utc
  }
  names(d)[names(d) == "value"] <- reg$column
  if (variable == "surface_sensible_heat_flux") check_era5_h_sign(d)
  d
}

#' Diurnal sign check of the ERA5 sensible heat flux
#'
#' The mean diurnal cycle of `H_era5_Wm2` (positive upward), by hour of local
#' solar time of the interval midpoint (UTC + longitude / 15 h; equation of
#' time neglected, < 17 min), must peak within `tol_h` of solar noon in every
#' cell. A failure indicates a missing or double sign flip.
#'
#' @param d Output of [read_era5_single_level()] for
#'   `surface_sensible_heat_flux`.
#' @param tol_h Tolerance around solar noon [h]; default 3.
#' @return Tibble invisibly: `cell_row`, `cell_col`, `peak_solar_hour` [h,
#'   bin centre]. Aborts if any cell fails.
check_era5_h_sign <- function(d, tol_h = 3) {
  mid <- d$timestamp_end_utc - 1800
  hod <- as.numeric(mid - lubridate::floor_date(mid, "day"), units = "hours")
  solar_h <- (hod + d$grid_lon / 15) %% 24
  d$bin <- floor(solar_h)
  cyc <- dplyr::summarise(d, .by = c(cell_row, cell_col, bin), H = mean(H_era5_Wm2, na.rm = TRUE))
  peak <- dplyr::summarise(
    cyc[is.finite(cyc$H), ], .by = c(cell_row, cell_col),
    peak_solar_hour = bin[which.max(H)] + 0.5
  )
  off <- abs(((peak$peak_solar_hour - 12) + 12) %% 24 - 12)
  if (any(off > tol_h)) {
    bad <- peak[off > tol_h, ]
    cli::cli_abort(c(
      "ERA5 sensible heat flux: mean diurnal maximum not within {tol_h} h of solar noon.",
      "x" = "Cells (row col / peak solar hour): {.val {paste0(bad$cell_row, ' ', bad$cell_col, ' / ', bad$peak_solar_hour)}}.",
      "i" = "Expected H_era5_Wm2 = -sshf / 3600 (sshf positive downward); check the sign convention."
    ))
  }
  invisible(peak)
}

#' Read ERA5 pressure-level temperature and geopotential
#'
#' Reads the wide files (`temperature_<p>hPa` [K], `geopotential_<p>hPa`
#' [m2 s-2]) of the selected cells, reshapes them to long format, joins
#' temperature and geopotential on cell x time x level (abort on mismatch) and
#' derives geopotential height `zg_m = z / g` (g = 9.80665 m s-2; WMO
#' geopotential height, within 0.5 % of geometric height below 10 km) and
#' potential temperature `theta_K = T_K * (1000 / p_hPa)^(r_d / cp)`
#' (Poisson's equation, e.g. Wallace & Hobbs 2006, Eq. 3.54; r_d and cp from
#' [pblh_constants()]). Checks: levels in the columns vs header `levels`,
#' rows and times vs header (abort), duplicated times (abort), 6-hourly gaps
#' (warn), median `T_K` at 1000 hPa in [230, 320] K (abort), `zg_m` increasing
#' with decreasing pressure (abort).
#'
#' @param files Output of [find_era5_files()] (family `"pressure_level"`).
#' @param cells `"all"` or `"center"`.
#' @param const Constants from [pblh_constants()].
#' @return Long tibble: `site_id`, `cell_row`, `cell_col` [integer],
#'   `timestamp_utc` [POSIXct UTC, instantaneous], `p_hPa` [hPa], `T_K` [K],
#'   `z_m2s2` [m2 s-2], `zg_m` [m], `theta_K` [K].
read_era5_pressure_levels <- function(files, cells = c("all", "center"), const = pblh_constants()) {
  cells <- match.arg(cells)
  f <- files[files$family == "pressure_level", ]
  if (cells == "center") f <- f[f$cell_row == 0 & f$cell_col == 0, ]
  for (v in c("temperature", "geopotential")) {
    if (!any(f$variable == v)) cli::cli_abort("No ERA5 pressure-level files for {.val {v}}.")
  }

  read_one <- function(i) {
    x <- read_era5_data(f[i, ], "time")
    v <- f$variable[i]
    cols <- setdiff(names(x), "timestamp_utc")
    mm <- stringr::str_match(cols, paste0("^", v, "_(\\d+)hPa$"))
    if (anyNA(mm[, 1])) {
      cli::cli_abort("{.file {f$file[i]}}: unexpected column{?s} {.val {cols[is.na(mm[, 1])]}}.")
    }
    lev <- as.integer(mm[, 2])
    hlev <- f$header[[i]]$levels
    if (!setequal(lev, hlev)) {
      cli::cli_abort("{.file {f$file[i]}}: levels in the columns ({.val {lev}}) differ from header {.field levels} ({.val {hlev}}).")
    }
    n <- nrow(x)
    dplyr::tibble(
      site_id = f$site_id[i], cell_row = f$cell_row[i], cell_col = f$cell_col[i],
      timestamp_utc = rep(x$timestamp_utc, times = length(lev)),
      p_hPa = rep(lev, each = n),
      value = unlist(x[cols], use.names = FALSE)
    )
  }
  parts <- lapply(seq_len(nrow(f)), read_one)
  tk <- dplyr::bind_rows(parts[f$variable == "temperature"])
  zz <- dplyr::bind_rows(parts[f$variable == "geopotential"])
  label <- paste("ERA5", f$site_id[1], "pressure levels")
  check_era5_times(tk[tk$p_hPa == tk$p_hPa[1], ], 6 * 3600, paste(label, "temperature"))
  check_era5_times(zz[zz$p_hPa == zz$p_hPa[1], ], 6 * 3600, paste(label, "geopotential"))

  by <- c("site_id", "cell_row", "cell_col", "timestamp_utc", "p_hPa")
  n_t_only <- nrow(dplyr::anti_join(tk, zz, by = by))
  n_z_only <- nrow(dplyr::anti_join(zz, tk, by = by))
  if (n_t_only + n_z_only > 0) {
    cli::cli_abort("{label}: temperature and geopotential do not match ({n_t_only} cell x time x level only in temperature, {n_z_only} only in geopotential).")
  }
  d <- dplyr::inner_join(dplyr::rename(tk, T_K = value), dplyr::rename(zz, z_m2s2 = value), by = by)
  d$zg_m <- d$z_m2s2 / const$g                                   # geopotential height
  d$theta_K <- d$T_K * (1000 / d$p_hPa)^(const$r_d / const$cp) # Poisson's equation
  d <- d[order(d$cell_row, d$cell_col, d$timestamp_utc, -d$p_hPa), ]

  if (1000 %in% d$p_hPa) {
    med <- stats::median(d$T_K[d$p_hPa == 1000], na.rm = TRUE)
    if (!is.finite(med) || med < 230 || med > 320) {
      cli::cli_abort("{label}: median T at 1000 hPa = {round(med, 1)} outside [230, 320] K; check the unit.")
    }
  }
  grp <- paste(d$cell_row, d$cell_col, as.numeric(d$timestamp_utc))
  dz <- stats::ave(d$zg_m, grp, FUN = function(z) c(Inf, diff(z)))
  n_bad <- sum(dz <= 0, na.rm = TRUE)
  if (n_bad > 0) {
    cli::cli_abort("{label}: geopotential height not increasing with decreasing pressure in {n_bad} level pair{?s}.")
  }
  d
}

#' Read ERA5 static fields of the 3x3 cells
#'
#' Surface geopotential `z` (orography height `orog_m = z / g`, g = 9.80665
#' m s-2), land-sea mask `lsm` and standard deviation of sub-grid orography
#' `sdor` (ECMWF parameters 129, 172, 160). The header `units` must match the
#' registry ([era5_variables()]); otherwise abort.
#'
#' @param site_id Station code [character].
#' @param cfg Configuration list with `era5_static_dir`.
#' @param sites Optional output of [read_sites()] for the coordinate check.
#' @param const Constants from [pblh_constants()].
#' @return Tibble, one row per cell: `site_id`, `cell_row`, `cell_col`
#'   [integer], `grid_lat`, `grid_lon` [deg], `grid_dist_km` [km], `orog_m`
#'   [m a.s.l., model orography], `lsm` [0-1], `sdor_m` [m].
read_era5_static <- function(site_id, cfg, sites = NULL, const = pblh_constants()) {
  f <- find_era5_files(site_id, cfg, family = "static", sites = sites)
  reg <- era5_variables()
  reg <- reg[reg$family == "static", ]
  vals <- lapply(seq_len(nrow(f)), function(i) {
    r <- reg[reg$variable == f$variable[i], ]
    u <- f$header[[i]]$units
    if (is.null(u) || u != r$header_units) {
      cli::cli_abort("{.file {f$file[i]}}: header units {.val {u}}, expected {.val {r$header_units}}.")
    }
    x <- read_era5_data(f[i, ], "valid_time_utc")
    if (nrow(x) != 1 || !f$variable[i] %in% names(x)) {
      cli::cli_abort("{.file {f$file[i]}}: expected one row with column {.field {f$variable[i]}}.")
    }
    v <- x[[f$variable[i]]]
    if (f$variable[i] == "z") v <- v / const$g # surface geopotential -> orography height
    dplyr::tibble(
      variable = r$column, cell_row = f$cell_row[i], cell_col = f$cell_col[i],
      grid_lat = f$grid_lat[i], grid_lon = f$grid_lon[i], grid_dist_km = f$grid_dist_km[i], value = v
    )
  })
  v <- dplyr::bind_rows(vals)
  pos <- dplyr::distinct(v, cell_row, cell_col, grid_lat, grid_lon, grid_dist_km)
  if (anyDuplicated(pos[c("cell_row", "cell_col")])) {
    cli::cli_abort("ERA5 static {.val {site_id}}: cell centres differ between variables.")
  }
  out <- pos
  for (col in c("orog_m", "lsm", "sdor_m")) {
    x <- v[v$variable == col, c("cell_row", "cell_col", "value")]
    names(x)[3] <- col
    out <- dplyr::left_join(out, x, by = c("cell_row", "cell_col"))
  }
  out <- dplyr::bind_cols(dplyr::tibble(site_id = rep(site_id, nrow(out))), out)
  out[order(out$cell_row, out$cell_col), ]
}
