# Footprint computations with the FFP model of Kljun et al. (2015), `Kljun2015` (spec 08).
# Kljun's R code (calc_footprint_FFP.R v1.3, calc_footprint_FFP_climatology.R v1.31, (c) N. Kljun) is local only
# (cfg$ffp_code_dir, gitignored): it is sourced into a dedicated environment, never copied or modified.
# All calls pass umean only (z0 = NaN): deviation from ETC, which passes both (the single-footprint function then
# scales x with z0 and Psi_M, the climatology with umean; spec 08, G 2026-10-26).

#' Load Kljun's FFP code
#'
#' @param dir Folder with `calc_footprint_FFP.R` and `calc_footprint_FFP_climatology.R`.
#' @return Environment with `calc_footprint_FFP` and `calc_footprint_FFP_climatology`,
#'   or `NULL` (with a message) if the files are missing.
load_ffp_code <- function(dir) {
  f <- file.path(dir, c("calc_footprint_FFP.R", "calc_footprint_FFP_climatology.R"))
  if (!all(file.exists(f))) {
    cli::cli_inform("Kljun FFP code not found in {.path {dir}} (cfg$ffp_code_dir; local only, not in the repository).")
    return(NULL)
  }
  env <- new.env(parent = baseenv())
  # the climatology file calls require(EBImage) at the top: not needed with smooth_data = 0 (smoothing done by
  # ffp_smooth()), so its warning is suppressed.
  for (p in f) suppressWarnings(suppressMessages(sys.source(p, envir = env)))
  if (!all(c("calc_footprint_FFP", "calc_footprint_FFP_climatology") %in% ls(env))) {
    cli::cli_abort("FFP functions not found after sourcing {.path {dir}}.")
  }
  env
}

# Constants of the crosswind-integrated parameterisation (Kljun et al. 2015, Eq. 17, as in the code).
ffp_constants <- function() list(a = 1.4524, b = -1.9914, c = 1.4622, d = 0.1359, k = 0.4, xstar_end = 30)

#' Crosswind-integrated footprint distances from the paper's equations
#'
#' Kljun et al. (2015): scaled footprint F* = a (X* - d)^b exp(-c / (X* - d)) (Eq. 14,
#' constants Eq. 17) on the code's X* grid (`nx` + 1 points from d to 30, first point
#' dropped), real scale with umean (Eq. 6): x = X* zm / (1 - zm/h) * umean / u* * k,
#' f = F* / (zm / (1 - zm/h) * umean / u* * k). Distances as ETC
#' (`legacy/FFP_calc_parallel_dev.R`): peak = argmax of f, x_r = x where cumsum(f) is
#' closest to r sum(f). Because the X* grid does not depend on the inputs, x_r = X*_r times
#' the scale factor; verified against the code (spec 08).
#'
#' @param zm Measurement height above d [m].
#' @param h Boundary layer height [m].
#' @param umean Mean wind speed [m s-1].
#' @param ustar Friction velocity [m s-1].
#' @param r Fractions [0-1].
#' @param nx Number of X* grid points (code default 1000).
#' @return Tibble: `x_peak`, `x_<100 r>` [m, upwind distance]; `NA` where the inputs
#'   are missing or zm >= h.
ffp_ci_paper <- function(zm, h, umean, ustar, r = c(0.5, 0.7, 0.8, 0.9), nx = 1000) {
  k <- ffp_constants()
  xs <- seq(k$d, k$xstar_end, length.out = nx + 2)[-1]
  fs <- k$a * (xs - k$d)^k$b * exp(-k$c / (xs - k$d))
  ok <- !is.na(fs)
  xs <- xs[ok]
  fs <- fs[ok]
  cs <- cumsum(fs)
  i_r <- vapply(r, function(p) which.min(abs(cs - p * sum(fs))), integer(1))
  scale <- zm / (1 - zm / h) * umean / ustar * k$k
  scale[!(zm < h)] <- NA_real_
  out <- dplyr::tibble(x_peak = xs[which.max(fs)] * scale)
  for (j in seq_along(r)) out[[paste0("x_", round(100 * r[j]))]] <- xs[i_r[j]] * scale
  out
}

#' Crosswind-integrated footprint distances from Kljun's code (one half-hour)
#'
#' Used to verify [ffp_ci_paper()]. `calc_footprint_FFP()` with umean only (z0 = NaN),
#' `rslayer = 1`, no contours; peak = argmax of `f_ci`, x_r by cumulative `f_ci` (ETC).
#'
#' @param env Output of [load_ffp_code()].
#' @inheritParams ffp_ci_paper
#' @param L Obukhov length [m].
#' @param sigmav Lateral velocity standard deviation [m s-1].
#' @param wind_dir Wind direction [deg].
#' @return One-row tibble as [ffp_ci_paper()], plus `flag_err`.
ffp_ci_code <- function(env, zm, h, umean, L, sigmav, ustar, wind_dir, r = c(0.5, 0.7, 0.8, 0.9)) {
  junk <- utils::capture.output(f <- env$calc_footprint_FFP(zm = zm, z0 = NaN, umean = umean, h = h, ol = L,
                                                          sigmav = sigmav, ustar = ustar, wind_dir = wind_dir,
                                                          r = NaN, rslayer = 1))
  out <- dplyr::tibble(x_peak = NA_real_)
  for (p in r) out[[paste0("x_", round(100 * p))]] <- NA_real_
  out$flag_err <- f$flag_err
  if (f$flag_err == 0) {
    x <- as.numeric(f$x_ci)
    fc <- as.numeric(f$f_ci)
    cs <- cumsum(fc)
    out$x_peak <- x[which.max(fc)]
    for (p in r) out[[paste0("x_", round(100 * p))]] <- x[which.min(abs(cs - p * sum(fc)))]
  }
  out
}

#' Fixed footprint grid
#'
#' As `calc_footprint_FFP_climatology()`: x = seq(xmin, xmax, by = dx), y likewise;
#' matrices have one row per y and one column per x (x east, y north, tower at 0, 0).
#'
#' @param domain c(xmin, xmax, ymin, ymax) [m].
#' @param dx Cell size [m] (dy = dx).
#' @return List: `x`, `y` [m], `rho` [m, distance from the tower, matrix], `dx`, `domain`.
ffp_grid <- function(domain = c(-1000, 1000, -1000, 1000), dx = 1) {
  x <- seq(domain[1], domain[2], by = dx)
  y <- seq(domain[3], domain[4], by = dx)
  rho <- sqrt(outer(y^2, x^2, "+"))
  list(x = x, y = y, rho = rho, dx = dx, domain = domain)
}

#' Footprint field of one half-hour on the fixed grid (Kljun's code)
#'
#' `calc_footprint_FFP_climatology()` called for one half-hour, umean only (z0 = NaN),
#' `rslayer = 1`, `smooth_data = 0` (smoothing by [ffp_smooth()]), no contours. The
#' climatology of several half-hours is the mean of these fields (the code sums the
#' fields and divides by the number of valid half-hours), smoothed afterwards.
#'
#' @param env Output of [load_ffp_code()].
#' @param grid Output of [ffp_grid()].
#' @inheritParams ffp_ci_code
#' @return Matrix [m-2] (rows y, columns x), or `NULL` if `flag_err != 0`.
ffp_field_halfhour <- function(env, grid, zm, h, umean, L, sigmav, ustar, wind_dir) {
  junk <- utils::capture.output(f <- without_ebimage_notice(env$calc_footprint_FFP_climatology(
    zm = zm, z0 = NaN, umean = umean, h = h, ol = L, sigmav = sigmav, ustar = ustar, wind_dir = wind_dir,
    domain = grid$domain, dx = grid$dx, dy = grid$dx, r = NaN, rslayer = 1, smooth_data = 0, pulse = 0
  )))
  if (!identical(f$flag_err, 0) && !identical(f$flag_err, 0L)) return(NULL)
  f$fclim_2d
}

# calc_footprint_FFP_climatology() calls require(EBImage) at every call; EBImage is not needed with smooth_data = 0.
# Mute only that message and warning.
without_ebimage_notice <- function(expr) {
  withCallingHandlers(expr,
    warning = function(w) if (grepl("EBImage", conditionMessage(w))) invokeRestart("muffleWarning"),
    message = function(m) if (grepl("EBImage", conditionMessage(m))) invokeRestart("muffleMessage")
  )
}

#' Smoothing of a footprint field as in calc_footprint_FFP_climatology()
#'
#' The code applies the 3x3 kernel (0.05 0.1 0.05 / 0.1 0.4 0.1 / 0.05 0.1 0.05) twice,
#' each time padding the field with one row/column of zeros, filtering with
#' `EBImage::filter2()` and cropping back. `filter2()`'s circular boundary then only
#' wraps the zero padding, so this equals a zero-padded 3x3 convolution (the kernel
#' is symmetric); verified in the tests against an FFT circular convolution.
#'
#' @param f Footprint field (matrix).
#' @param times Number of passes (code: 2).
#' @return Smoothed field, same size.
ffp_smooth <- function(f, times = 2) {
  k <- matrix(c(0.05, 0.1, 0.05, 0.1, 0.4, 0.1, 0.05, 0.1, 0.05), 3)
  nr <- nrow(f)
  nc <- ncol(f)
  for (t in seq_len(times)) {
    p <- matrix(0, nr + 2, nc + 2)
    p[2:(nr + 1), 2:(nc + 1)] <- f
    o <- matrix(0, nr, nc)
    for (a in 1:3) for (b in 1:3) o <- o + k[a, b] * p[a:(a + nr - 1), b:(b + nc - 1)]
    f <- o
  }
  f
}

#' Metrics of a (smoothed) footprint field
#'
#' Source regions as the code's contours: for fraction r, the threshold fr is the
#' value at which the cumulative sum of the sorted field (times dx dy) is closest to r;
#' the region is the set of cells with f >= fr (Kljun et al. 2015, Sect. 5.3).
#'
#' @param f Footprint field [m-2] on `grid`.
#' @param grid Output of [ffp_grid()].
#' @param mask Target-area mask (logical matrix on `grid`), or `NULL`.
#' @param r Fractions [0-1].
#' @param keep_region Fraction whose region is returned as cell indices (e.g. 0.8), or `NULL`.
#' @return List: `metrics` (one-row tibble: `f_domain` [-, footprint fraction inside the
#'   domain], `f_target_area` [-], `area_<r>` [m2], `xmax_<r>` [m, largest distance from the
#'   tower], `edge_<r>` [logical, region touches the domain edge]) and `region` (integer
#'   cell indices, or `NULL`).
ffp_field_metrics <- function(f, grid, mask = NULL, r = c(0.5, 0.7, 0.8, 0.9), keep_region = 0.8) {
  da <- grid$dx^2
  v <- c(f)
  o <- order(v, decreasing = TRUE)
  cs <- cumsum(v[o]) * da
  m <- dplyr::tibble(f_domain = sum(v) * da, f_target_area = if (is.null(mask)) NA_real_ else sum(v[c(mask)]) * da)
  nr <- nrow(f)
  nc <- ncol(f)
  edge <- rep(FALSE, length(v))
  edge[c(row(f) == 1 | row(f) == nr | col(f) == 1 | col(f) == nc)] <- TRUE
  region <- NULL
  for (p in r) {
    n_in <- which.min(abs(cs - p))
    idx <- o[seq_len(n_in)]
    lab <- round(100 * p)
    m[[paste0("area_", lab)]] <- n_in * da
    m[[paste0("xmax_", lab)]] <- max(grid$rho[idx])
    m[[paste0("edge_", lab)]] <- any(edge[idx])
    if (!is.null(keep_region) && isTRUE(all.equal(p, keep_region))) region <- sort(idx)
  }
  list(metrics = m, region = region)
}

#' Jaccard index of two regions
#'
#' @param a,b Integer cell indices.
#' @return |a and b| / |a or b| [-].
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

#' Point in polygon (even-odd rule)
#'
#' @param px,py Point coordinates [m].
#' @param rings List of rings, each a two-column matrix of vertices [m] (outer rings and
#'   holes; the even-odd rule handles holes and multipolygons).
#' @return Logical, one per point.
point_in_rings <- function(px, py, rings) {
  inside <- rep(FALSE, length(px))
  for (rg in rings) {
    vx <- rg[, 1]
    vy <- rg[, 2]
    n <- length(vx)
    j <- n
    for (i in seq_len(n)) {
      cross <- ((vy[i] > py) != (vy[j] > py)) & (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])
      inside <- xor(inside, cross)
      j <- i
    }
  }
  inside
}

#' Target-area mask on the footprint grid
#'
#' Reads the ICOS target-area KML, projects it to a local azimuthal equidistant
#' projection centred on the tower (metres, x east, y north, matching the FFP grid) and
#' marks the grid cells whose centre lies inside.
#'
#' @param kml Path of `<site>_TA.kml`.
#' @param lat,lon Tower coordinates [decimal degrees].
#' @param grid Output of [ffp_grid()].
#' @return Logical matrix on `grid`; attribute `area_m2` (polygon area).
target_area_mask <- function(kml, lat, lon, grid) {
  p <- sf::st_zm(sf::st_read(kml, quiet = TRUE))
  crs <- sprintf("+proj=aeqd +lat_0=%.8f +lon_0=%.8f +x_0=0 +y_0=0 +datum=WGS84 +units=m", lat, lon)
  p <- sf::st_transform(p, crs)
  co <- sf::st_coordinates(p)
  grp <- interaction(as.data.frame(co[, grepl("^L[0-9]", colnames(co)), drop = FALSE]), drop = TRUE)
  rings <- lapply(split(seq_len(nrow(co)), grp), function(i) co[i, c("X", "Y"), drop = FALSE])
  px <- rep(grid$x, each = length(grid$y))
  py <- rep(grid$y, times = length(grid$x))
  m <- matrix(point_in_rings(px, py, rings), nrow = length(grid$y), ncol = length(grid$x))
  attr(m, "area_m2") <- as.numeric(sum(sf::st_area(p)))
  m
}
