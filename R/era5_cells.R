# ERA5 grid cells around a station: cell metadata, bilinear interpolation to
# the tower position, 3x3 spread as representativeness metric.
# Rules: docs/specs/04_era5_standard_compare.md section 3.

#' ERA5 cell metadata of a station
#'
#' Describes the ERA5 tower cell (position, distance, model orography vs site
#' altitude, land-sea mask, sub-grid orography) and the range of the static
#' fields over the 3x3 block.
#'
#' @param site_id Station code [character].
#' @param cfg Configuration list with `era5_static_dir`.
#' @param sites Output of [read_sites()] (`lat`, `lon` [deg], `altitude_m` [m a.s.l.]).
#' @param static Optional output of [read_era5_static()]; read if `NULL`.
#' @return Tibble, one row: `site_id`, `site_lat`, `site_lon` [deg],
#'   `site_altitude_m` [m a.s.l.], `cell_lat`, `cell_lon` tower-cell centre
#'   [deg], `cell_dist_km` tower to cell centre [km], `orog_m` tower-cell model
#'   orography [m a.s.l.], `dalt_m = site_altitude_m - orog_m` [m], `lsm`
#'   [0-1], `sdor_m` [m]; over the 3x3 block: `orog_min_m`, `orog_max_m` [m],
#'   `lsm_min`, `lsm_max` [0-1], `sdor_max_m` [m], `n_cells` [-].
era5_cell_metadata <- function(site_id, cfg, sites, static = NULL) {
  meta <- sites[sites$site_id == site_id, ]
  if (nrow(meta) != 1) {
    cli::cli_abort("Site {.val {site_id}} not found in site metadata.")
  }
  if (is.null(static)) static <- read_era5_static(site_id, cfg, sites)
  tw <- static[static$cell_row == 0 & static$cell_col == 0, ]
  if (nrow(tw) != 1) {
    cli::cli_abort("ERA5 static {.val {site_id}}: tower cell missing.")
  }
  rng <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
  dplyr::tibble(
    site_id = site_id,
    site_lat = meta$lat,
    site_lon = meta$lon,
    site_altitude_m = meta$altitude_m,
    cell_lat = tw$grid_lat,
    cell_lon = tw$grid_lon,
    cell_dist_km = haversine_km(meta$lat, meta$lon, tw$grid_lat, tw$grid_lon),
    orog_m = tw$orog_m,
    dalt_m = meta$altitude_m - tw$orog_m,
    lsm = tw$lsm,
    sdor_m = tw$sdor_m,
    orog_min_m = rng(static$orog_m, min),
    orog_max_m = rng(static$orog_m, max),
    lsm_min = rng(static$lsm, min),
    lsm_max = rng(static$lsm, max),
    sdor_max_m = rng(static$sdor_m, max),
    n_cells = nrow(static)
  )
}

#' Bilinear interpolation of an ERA5 field to the tower position
#'
#' Standard bilinear interpolation on the regular lat-lon grid (e.g. Press et
#' al. 2007, Numerical Recipes, 3rd ed., Sect. 3.6) from the 4 cell centres
#' surrounding the tower: the tower cell (row 0, col 0) and its neighbours in
#' the direction of the sign of `site - cell centre` in latitude (row) and
#' longitude (col). With fractional distances `fy = (site_lat - lat0) / (lat1 -
#' lat0)` and `fx = (site_lon - lon0) / (lon1 - lon0)` the weights are
#' `(1 - fy)(1 - fx)`, `fy (1 - fx)`, `(1 - fy) fx`, `fy fx` (sum 1). A time
#' with a missing value in a cell of non-zero weight gives `NA`.
#'
#' @param era5_long Long ERA5 tibble with `cell_row`, `cell_col`, `grid_lat`,
#'   `grid_lon` [deg], `timestamp_utc` and `value_col` (e.g. output of
#'   [read_era5_single_level()] with `cells = "all"`).
#' @param site_lat,site_lon Tower position [deg].
#' @param value_col Name of the column to interpolate [character].
#' @return Tibble `timestamp_utc`, `<value_col>_bilin` [unit of `value_col`],
#'   with attribute `weights` (tibble `cell_row`, `cell_col`, `weight` [-]).
era5_bilinear <- function(era5_long, site_lat, site_lon, value_col) {
  if (!value_col %in% names(era5_long)) {
    cli::cli_abort("Column {.field {value_col}} not found.")
  }
  pos <- dplyr::distinct(era5_long, cell_row, cell_col, grid_lat, grid_lon)
  p0 <- pos[pos$cell_row == 0 & pos$cell_col == 0, ]
  if (nrow(p0) != 1) cli::cli_abort("ERA5 bilinear: tower cell (row 0, col 0) missing.")
  dr <- if (site_lat >= p0$grid_lat) 1L else -1L
  dc <- if (site_lon >= p0$grid_lon) 1L else -1L
  lat1 <- p0$grid_lat + dr * era5_grid_res_deg
  lon1 <- p0$grid_lon + dc * era5_grid_res_deg
  fy <- (site_lat - p0$grid_lat) / (lat1 - p0$grid_lat)
  fx <- (site_lon - p0$grid_lon) / (lon1 - p0$grid_lon)
  if (fy < 0 || fy > 1 || fx < 0 || fx > 1) {
    cli::cli_abort("ERA5 bilinear: tower ({site_lat}, {site_lon}) is not inside the tower cell's neighbourhood.")
  }
  w <- dplyr::tibble(
    cell_row = c(0L, dr, 0L, dr),
    cell_col = c(0L, 0L, dc, dc),
    weight = c((1 - fy) * (1 - fx), fy * (1 - fx), (1 - fy) * fx, fy * fx)
  )

  times <- sort(unique(era5_long$timestamp_utc))
  out <- numeric(length(times))
  for (k in seq_len(nrow(w))) {
    sel <- era5_long$cell_row == w$cell_row[k] & era5_long$cell_col == w$cell_col[k]
    if (w$weight[k] == 0) next
    if (!any(sel)) {
      cli::cli_warn("ERA5 bilinear: cell (row {w$cell_row[k]}, col {w$cell_col[k]}) missing; result NA.")
      out[] <- NA_real_
      next
    }
    v <- era5_long[[value_col]][sel][match(times, era5_long$timestamp_utc[sel])]
    out <- out + w$weight[k] * v
  }
  res <- dplyr::tibble(timestamp_utc = times, value = out)
  names(res)[2] <- paste0(value_col, "_bilin")
  attr(res, "weights") <- w
  res
}

#' Spread of an ERA5 field over the 3x3 cells
#'
#' Per time, statistics over the available cells (`n9 < 9` when cells or
#' values are missing). The coefficient of variation is a representativeness
#' metric of the tower cell (spec 04; `docs/project_plan.md`, ERA5
#' representativeness).
#'
#' @param era5_long Long ERA5 tibble with `timestamp_utc` and `value_col`.
#' @param value_col Name of the column [character].
#' @return Tibble per `timestamp_utc`: `mean9`, `sd9` (sample sd, `NA` if
#'   `n9 < 2`), `min9`, `max9`, `range9 = max9 - min9` [unit of `value_col`],
#'   `cv9 = sd9 / mean9` [-], `n9` [number of non-missing cells].
era5_spread <- function(era5_long, value_col) {
  if (!value_col %in% names(era5_long)) {
    cli::cli_abort("Column {.field {value_col}} not found.")
  }
  d <- dplyr::tibble(timestamp_utc = era5_long$timestamp_utc, v = era5_long[[value_col]])
  d <- d[!is.na(d$v), ]
  if (nrow(d) == 0) {
    return(dplyr::tibble(
      timestamp_utc = sort(unique(era5_long$timestamp_utc)), mean9 = NA_real_, sd9 = NA_real_,
      min9 = NA_real_, max9 = NA_real_, range9 = NA_real_, cv9 = NA_real_, n9 = 0L
    ))
  }
  out <- dplyr::summarise(
    d, .by = timestamp_utc,
    mean9 = mean(v),
    sd9 = if (dplyr::n() > 1) stats::sd(v) else NA_real_,
    min9 = min(v),
    max9 = max(v),
    n9 = dplyr::n()
  )
  out$range9 <- out$max9 - out$min9
  out$cv9 <- out$sd9 / out$mean9
  out <- out[order(out$timestamp_utc), c("timestamp_utc", "mean9", "sd9", "min9", "max9", "range9", "cv9", "n9")]
  # times where all cells are missing: keep them with n9 = 0
  all_t <- sort(unique(era5_long$timestamp_utc))
  if (length(all_t) > nrow(out)) {
    out <- dplyr::left_join(dplyr::tibble(timestamp_utc = all_t), out, by = "timestamp_utc")
    out$n9[is.na(out$n9)] <- 0L
  }
  out
}
