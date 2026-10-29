# Range compression of a PBLH source against a measured MLH (spec 10 A5): Theil-Sen slope of ln(source) on
# ln(measured) with a day-block bootstrap interval, and the median log ratio by tercile of the measured MLH.
# Slope < 1: the source compresses the range (too deep when the layer is shallow, too shallow when it is deep;
# candidate finding of spec 07).

#' Options of the comparison with a measured MLH (spec 10)
#'
#' @param aerosol_window Apparent solar time window [h, `c(start, end)`, start included, end excluded] for
#'   aerosol-based references (late morning to mid-afternoon: Schween2014, Kotthaus2018); applied to the interval
#'   midpoint (G 2026-10-29).
#' @param n_boot Bootstrap replicates of the Theil-Sen interval [-].
#' @param conf Confidence level of the bootstrap interval [0-1].
#' @param seed Random seed of the bootstrap [integer].
#' @param boot_tol Relative tolerance of the Theil-Sen slope in each bootstrap replicate ([theil_sen()]; the point
#'   estimate uses 1e-7) [-].
#' @return Named list.
measured_compare_options <- function(aerosol_window = c(10, 16), n_boot = 1000L, conf = 0.95, seed = 20261029L,
                                     boot_tol = 1e-4) {
  list(aerosol_window = aerosol_window, n_boot = as.integer(n_boot), conf = conf, seed = as.integer(seed),
       boot_tol = boot_tol)
}

#' Number of strict inversions of a sequence
#'
#' Pairs `i < j` with `q[i] > q[j]`, counted level by level as in a bottom-up merge sort: at block size `s`, each
#' element of the right half of a `2s` block is compared with the left half of the same block by `findInterval()`
#' on block-offset ranks (vectorised, O(n log^2 n)).
#'
#' @param q Numeric vector without `NA`.
#' @return Number of inversions [-, double].
count_inversions <- function(q) {
  n <- length(q)
  if (n < 2) return(0)
  q <- rank(q, ties.method = "min") # equal values are not inversions
  pos <- seq_len(n) - 1
  total <- 0
  s <- 1
  while (s < n) {
    blk <- pos %/% (2 * s)
    right <- (pos %% (2 * s)) >= s
    off <- blk * (n + 1)
    lk <- sort(off[!right] + q[!right])
    hi <- findInterval(off[right] + n, lk)          # left elements of the block
    lo <- findInterval(off[right] + q[right], lk)   # ... with rank <= q
    total <- total + sum(as.numeric(hi - lo))
    s <- 2 * s
  }
  total
}

#' Theil-Sen slope
#'
#' Median of the pairwise slopes `(y_j - y_i) / (x_j - x_i)` over the pairs with `x_i != x_j` (Theil 1950; Sen
#' 1968), without forming all pairs: with the points ordered by `(x, y)`, the number of slopes below `b` is the
#' number of inversions of `y - b x` ([count_inversions()]; within a group of equal `x` the `y` are sorted, so such
#' pairs never count), and each order statistic is found by bisection on `b` (even number of slopes: mean of the
#' two middle ones). The bisection starts from a bracket taken from a deterministic sample of pairwise slopes,
#' checked by the exact counts and widened if needed (speed only; no random numbers). Intercept = median of
#' `y - slope x`.
#'
#' @param x,y Numeric vectors (pairs with a non-finite value dropped).
#' @param tol Tolerance of the slope, relative to `max(1, |slope|)` [-].
#' @return List: `slope` [unit of y / unit of x], `intercept` [unit of y], `n` (points), `n_pairs` (pairs with
#'   distinct x) [-]; `NA` slope with fewer than 2 distinct x.
theil_sen <- function(x, y, tol = 1e-7) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  o <- order(x, y)
  x <- x[o]
  y <- y[o]
  g <- as.numeric(table(x))
  n_pairs <- n * (n - 1) / 2 - sum(g * (g - 1) / 2)
  if (n_pairs < 1) return(list(slope = NA_real_, intercept = NA_real_, n = n, n_pairs = 0))
  big <- (diff(range(y)) + 1) / min(diff(unique(x))) + 1 # bound on |slope|
  below <- function(b) count_inversions(y - b * x)       # number of slopes < b
  # deterministic sample of pairs (i, i + d): offsets spread over the x order
  i <- rep(seq_len(n), times = 7)
  j <- i + rep(pmax(1L, round(n * c(0.03, 0.11, 0.23, 0.37, 0.51, 0.67, 0.83))), each = n)
  keep <- j <= n
  i <- i[keep]
  j <- j[keep]
  dx <- x[j] - x[i]
  s_smp <- ((y[j] - y[i]) / dx)[dx != 0]
  kth <- function(k) { # k-th smallest slope: bracket [lo, hi] with F(lo) < k <= F(hi), F(b) = below(b)
    p <- k / n_pairs
    if (length(s_smp) >= 50) {
      lo <- stats::quantile(s_smp, max(0, p - 0.01), names = FALSE)
      hi <- stats::quantile(s_smp, min(1, p + 0.01), names = FALSE)
    } else {
      lo <- -big
      hi <- big
    }
    w <- max(hi - lo, tol)
    while (below(lo) >= k && lo > -big) {
      lo <- max(-big, lo - w)
      w <- 2 * w
    }
    while ((f_hi <- below(hi)) < k && hi < big) {
      hi <- min(big, hi + w)
      w <- 2 * w
    }
    while (hi - lo > tol * max(1, abs(lo), abs(hi))) {
      mid <- (lo + hi) / 2
      f_mid <- below(mid)
      if (f_mid >= k) {
        hi <- mid
        f_hi <- f_mid
      } else {
        lo <- mid
      }
    }
    list(value = (lo + hi) / 2, f_hi = f_hi)
  }
  k1 <- floor((n_pairs + 1) / 2)
  k2 <- ceiling((n_pairs + 1) / 2)
  s1 <- kth(k1)
  # even number of slopes: the (k1 + 1)-th lies in the same final bracket if F(hi) >= k2
  slope <- if (k1 == k2 || s1$f_hi >= k2) s1$value else (s1$value + kth(k2)$value) / 2
  list(slope = slope, intercept = stats::median(y - slope * x), n = n, n_pairs = n_pairs)
}

#' Range-compression slope with a day-block bootstrap interval
#'
#' Theil-Sen slope of `ln(mod)` on `ln(ref)` (spec 10 A5). Bootstrap: whole UTC days resampled with replacement
#' (G 2026-10-29: values within a day are autocorrelated), percentile interval. Also returned, as a diagnostic not
#' asked by the spec: the reverse slope of `ln(ref)` on `ln(mod)`. Random scatter alone (e.g. noise in the measured
#' values) pulls both slopes below 1; a source that compresses the range gives a forward slope < 1 and a reverse
#' slope > 1.
#'
#' @param mod,ref Source and measured PBLH [m] (pairs with a non-positive or missing value dropped).
#' @param t Times of the pairs [POSIXct UTC] (day blocks).
#' @param opts [measured_compare_options()].
#' @return One-row tibble: `n`, `n_days` [-], `slope`, `slope_lo`, `slope_hi` [-, ln-ln], `intercept` [ln m],
#'   `slope_reverse` [-], `n_boot` [-, successful replicates].
compression_slope <- function(mod, ref, t, opts = measured_compare_options()) {
  ok <- is.finite(mod) & is.finite(ref) & mod > 0 & ref > 0 & !is.na(t)
  lx <- log(ref[ok])
  ly <- log(mod[ok])
  day <- floor(as.numeric(t[ok]) / 86400)
  idx <- split(seq_along(day), day)
  out <- dplyr::tibble(n = sum(ok), n_days = length(idx), slope = NA_real_, slope_lo = NA_real_,
                       slope_hi = NA_real_, intercept = NA_real_, slope_reverse = NA_real_, n_boot = 0L)
  if (sum(ok) < 3 || length(unique(lx)) < 2 || length(unique(ly)) < 2) return(out)
  ts <- theil_sen(lx, ly)
  out$slope <- ts$slope
  out$intercept <- ts$intercept
  out$slope_reverse <- theil_sen(ly, lx)$slope
  if (opts$n_boot > 0 && length(idx) > 1) {
    has_seed <- exists(".Random.seed", envir = globalenv())
    if (has_seed) old <- get(".Random.seed", envir = globalenv())
    on.exit(if (has_seed) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv()))
    set.seed(opts$seed)
    boot <- vapply(seq_len(opts$n_boot), function(b) {
      i <- unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)
      theil_sen(lx[i], ly[i], tol = opts$boot_tol)$slope
    }, numeric(1))
    a <- (1 - opts$conf) / 2
    if (any(!is.na(boot))) {
      q <- stats::quantile(boot, c(a, 1 - a), na.rm = TRUE, names = FALSE)
      out$slope_lo <- q[1]
      out$slope_hi <- q[2]
    }
    out$n_boot <- sum(!is.na(boot))
  }
  out
}

#' Median log ratio by tercile of the measured MLH
#'
#' @param mod,ref Source and measured PBLH [m] (non-positive or missing pairs dropped).
#' @return Tibble, one row per tercile of `ref` (`T1` shallowest; limits at the 1/3 and 2/3 quantiles, upper limit
#'   included): `tercile`, `ref_lo`, `ref_hi` [m, measured range in the tercile], `n` [-], `median_log_ratio` [-,
#'   ln(mod / ref)], `frac_within_factor2` [0-1].
logratio_by_tercile <- function(mod, ref) {
  ok <- is.finite(mod) & is.finite(ref) & mod > 0 & ref > 0
  mod <- mod[ok]
  ref <- ref[ok]
  if (length(ref) < 3) {
    return(dplyr::tibble(tercile = character(), ref_lo = numeric(), ref_hi = numeric(), n = integer(),
                         median_log_ratio = numeric(), frac_within_factor2 = numeric()))
  }
  br <- stats::quantile(ref, c(1 / 3, 2 / 3), names = FALSE)
  k <- findInterval(ref, br, left.open = TRUE) + 1L
  lr <- log(mod / ref)
  dplyr::bind_rows(lapply(1:3, function(j) {
    s <- k == j
    dplyr::tibble(tercile = paste0("T", j), ref_lo = if (any(s)) min(ref[s]) else NA_real_,
                  ref_hi = if (any(s)) max(ref[s]) else NA_real_, n = sum(s),
                  median_log_ratio = if (any(s)) stats::median(lr[s]) else NA_real_,
                  frac_within_factor2 = if (any(s)) mean(abs(lr[s]) <= log(2)) else NA_real_)
  }))
}
