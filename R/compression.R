# Range compression of a PBLH source against a measured MLH (spec 10 A5; G 2026-10-29/30): a hypothesis under test
# (too deep when the layer is shallow, too shallow when it is deep; candidate finding of spec 07). Main evidence,
# where the measured MLH is censored (FI-Hyy): the scatter-only test of the censored tails (censored_tail_null(),
# one-sided and conservative). Supporting evidence: the symmetric slope (reduced major axis of ln(source) on
# ln(measured)) and the Bland-Altman trend of ln(source / measured) against the mean of the two logs, both with a
# day-block bootstrap interval; supporting only, because source and measurement differ in definition, not only in
# noise ("equation error", Warton2006), which the symmetric slope does not separate from a range effect. Forward and
# reverse Theil-Sen slopes are diagnostics only: scatter in the regressor pulls a one-sided slope below 1
# (regression dilution), so a forward slope < 1 alone is no evidence of compression; for the same reason the
# median log ratio by tercile of the measured MLH is affected by regression to the mean.

#' Options of the comparison with a measured MLH (spec 10)
#'
#' @param aerosol_window Apparent solar time window [h, `c(start, end)`, start included, end excluded] for
#'   aerosol-based references (late morning to mid-afternoon: Schween2014, Kotthaus2018); applied to the interval
#'   midpoint (G 2026-10-29).
#' @param n_boot Bootstrap replicates of the symmetric slope and of the Bland-Altman trend [-].
#' @param conf Confidence level of the bootstrap intervals [0-1].
#' @param seed Random seed of the bootstrap and of the censored-tail null [integer].
#' @param n_bins Bins of equal count of the Bland-Altman mean ([bland_altman_bins()]) [-].
#' @param n_sim Simulations of the censored-tail null ([censored_tail_null()]) [-].
#' @param min_pool Minimum number of exact-value log ratios in a season x regime pool of the null; smaller pools
#'   fall back to the regime, then to all exact pairs [-].
#' @return Named list.
measured_compare_options <- function(aerosol_window = c(10, 16), n_boot = 1000L, conf = 0.95, seed = 20261029L,
                                     n_bins = 10L, n_sim = 1000L, min_pool = 30L) {
  list(aerosol_window = aerosol_window, n_boot = as.integer(n_boot), conf = conf, seed = as.integer(seed),
       n_bins = as.integer(n_bins), n_sim = as.integer(n_sim), min_pool = as.integer(min_pool))
}

# Run `expr` with the given seed and restore the caller's random-number state afterwards.
with_seed <- function(seed, expr) {
  has_seed <- exists(".Random.seed", envir = globalenv())
  if (has_seed) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (has_seed) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
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

#' Symmetric (reduced major axis) slope
#'
#' `sign(r) sd(y) / sd(x)` (geometric mean of the slope of y on x and the inverse of the slope of x on y; standardised
#' major axis, Warton2006): treats both variables alike, so scatter in either does not pull it towards 0; equation
#' error (y and x differing in definition, not only in measurement noise) still affects it (Warton2006). Line through
#' the means.
#'
#' @param x,y Numeric vectors (pairs with a non-finite value dropped).
#' @return List: `slope` [unit of y / unit of x], `intercept` [unit of y], `r` [-, Pearson]; `NA` with fewer than 3
#'   pairs or a constant variable.
rma_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(list(slope = NA_real_, intercept = NA_real_, r = NA_real_))
  }
  r <- stats::cor(x, y)
  b <- sign(r) * stats::sd(y) / stats::sd(x)
  list(slope = b, intercept = mean(y) - b * mean(x), r = r)
}

#' Bland-Altman trend of the log ratio
#'
#' Least-squares line of `d = ln(mod) - ln(ref)` on `a = (ln(mod) + ln(ref)) / 2` (BlandAltman1999, regression of
#' the differences on the means). Slope 0 when `sd(ln mod) = sd(ln ref)`, i.e. symmetric slope 1; slope < 0 when the
#' source has the narrower spread (range compression): `slope = 2 (var y - var x) / (var x + var y + 2 cov)`.
#'
#' @param x,y `ln(ref)`, `ln(mod)` [ln m] (pairs with a non-finite value dropped).
#' @return List: `slope` [-], `intercept` [-, d at a = 0], `mean_d` [-, mean log ratio].
bland_altman_fit <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  d <- y[ok] - x[ok]
  a <- (y[ok] + x[ok]) / 2
  if (length(d) < 3 || stats::var(a) == 0) return(list(slope = NA_real_, intercept = NA_real_, mean_d = mean(d)))
  b <- stats::cov(a, d) / stats::var(a)
  list(slope = b, intercept = mean(d) - b * mean(a), mean_d = mean(d))
}

#' Range-compression statistics with a day-block bootstrap interval
#'
#' On `x = ln(ref)`, `y = ln(mod)` (spec 10 A5; G 2026-10-29): symmetric slope ([rma_slope()], main estimate) and
#' Bland-Altman trend ([bland_altman_fit()]), each with a percentile interval from a bootstrap of whole UTC days
#' (values within a day are autocorrelated); forward and reverse Theil-Sen slopes ([theil_sen()]) as point
#' estimates (diagnostics only). Regression dilution: random scatter in `ref` pulls the forward slope below 1, and
#' scatter in `mod` pulls the reverse slope below 1; a source that compresses the range gives a symmetric slope < 1
#' and a negative Bland-Altman slope.
#'
#' @param mod,ref Source and measured PBLH [m] (pairs with a non-positive or missing value dropped).
#' @param t Times of the pairs [POSIXct UTC] (day blocks).
#' @param opts [measured_compare_options()].
#' @param diagnostics Also compute the Theil-Sen diagnostics (`FALSE`: `ts_*` are `NA`; spec 11, many pairs) [logical].
#' @return One-row tibble: `n`, `n_days` [-], `r_log` [-, Pearson on the logs], `slope_sym`, `slope_sym_lo`,
#'   `slope_sym_hi` [-, ln-ln], `intercept_sym` [ln m], `ba_mean` [-, mean ln(mod / ref)], `ba_slope`, `ba_slope_lo`,
#'   `ba_slope_hi` [-], `ba_intercept` [-], `ts_slope`, `ts_slope_reverse` [-, Theil-Sen of ln mod on ln ref and of
#'   ln ref on ln mod], `ts_intercept` [ln m], `n_boot` [-, successful replicates]. Attribute `ba_band`: tibble of
#'   the Bland-Altman line on a grid of the log mean, `a` [ln m], `fit`, `lo`, `hi` [-] (pointwise percentile band).
compression_slope <- function(mod, ref, t, opts = measured_compare_options(), diagnostics = TRUE) {
  ok <- is.finite(mod) & is.finite(ref) & mod > 0 & ref > 0 & !is.na(t)
  lx <- log(ref[ok])
  ly <- log(mod[ok])
  day <- floor(as.numeric(t[ok]) / 86400)
  idx <- split(seq_along(day), day)
  out <- dplyr::tibble(n = sum(ok), n_days = length(idx), r_log = NA_real_, slope_sym = NA_real_,
                       slope_sym_lo = NA_real_, slope_sym_hi = NA_real_, intercept_sym = NA_real_,
                       ba_mean = if (any(ok)) mean(ly - lx) else NA_real_, ba_slope = NA_real_,
                       ba_slope_lo = NA_real_, ba_slope_hi = NA_real_, ba_intercept = NA_real_, ts_slope = NA_real_,
                       ts_slope_reverse = NA_real_, ts_intercept = NA_real_, n_boot = 0L)
  attr(out, "ba_band") <- dplyr::tibble(a = numeric(), fit = numeric(), lo = numeric(), hi = numeric())
  if (sum(ok) < 3 || length(unique(lx)) < 2 || length(unique(ly)) < 2) return(out)
  rm <- rma_slope(lx, ly)
  ba <- bland_altman_fit(lx, ly)
  out$r_log <- rm$r
  out$slope_sym <- rm$slope
  out$intercept_sym <- rm$intercept
  out$ba_slope <- ba$slope
  out$ba_intercept <- ba$intercept
  if (diagnostics) {
    ts <- theil_sen(lx, ly)
    out$ts_slope <- ts$slope
    out$ts_intercept <- ts$intercept
    out$ts_slope_reverse <- theil_sen(ly, lx)$slope
  }
  am <- (lx + ly) / 2
  grid <- seq(stats::quantile(am, 0.01, names = FALSE), stats::quantile(am, 0.99, names = FALSE), length.out = 40)
  band <- dplyr::tibble(a = grid, fit = ba$intercept + ba$slope * grid, lo = NA_real_, hi = NA_real_)
  if (opts$n_boot > 0 && length(idx) > 1) {
    boot <- with_seed(opts$seed, vapply(seq_len(opts$n_boot), function(b) {
      i <- unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)
      f <- bland_altman_fit(lx[i], ly[i])
      c(rma_slope(lx[i], ly[i])$slope, f$slope, f$intercept)
    }, numeric(3)))
    a <- (1 - opts$conf) / 2
    qq <- function(v) stats::quantile(v, c(a, 1 - a), na.rm = TRUE, names = FALSE)
    good <- !is.na(boot[1, ]) & !is.na(boot[2, ])
    if (any(good)) {
      out[c("slope_sym_lo", "slope_sym_hi")] <- as.list(qq(boot[1, good]))
      out[c("ba_slope_lo", "ba_slope_hi")] <- as.list(qq(boot[2, good]))
      lines <- outer(boot[3, good], rep(1, length(grid))) + outer(boot[2, good], grid)
      band$lo <- apply(lines, 2, stats::quantile, a, names = FALSE)
      band$hi <- apply(lines, 2, stats::quantile, 1 - a, names = FALSE)
    }
    out$n_boot <- sum(good)
  }
  attr(out, "ba_band") <- band
  out
}

#' Bland-Altman binned medians
#'
#' Median of `d = ln(mod / ref)` in bins of equal count of `a = (ln mod + ln ref) / 2` (G 2026-10-29: replaces the
#' terciles of the measured MLH, which are affected by regression to the mean).
#'
#' @param mod,ref Source and measured PBLH [m] (non-positive or missing pairs dropped).
#' @param n_bins Number of bins [-] (fewer if `a` has ties at the limits).
#' @return Tibble per bin: `bin` [1..], `a_lo`, `a_hi` [ln m], `gm_m` [m, exp of the median of `a`], `n` [-],
#'   `median_d`, `q25_d`, `q75_d` [-], `frac_within_factor2` [0-1].
bland_altman_bins <- function(mod, ref, n_bins = 10L) {
  ok <- is.finite(mod) & is.finite(ref) & mod > 0 & ref > 0
  d <- log(mod[ok] / ref[ok])
  a <- (log(mod[ok]) + log(ref[ok])) / 2
  if (length(d) < n_bins) {
    return(dplyr::tibble(bin = integer(), a_lo = numeric(), a_hi = numeric(), gm_m = numeric(), n = integer(),
                         median_d = numeric(), q25_d = numeric(), q75_d = numeric(), frac_within_factor2 = numeric()))
  }
  br <- unique(stats::quantile(a, seq(0, 1, length.out = n_bins + 1), names = FALSE))
  k <- findInterval(a, br, rightmost.closed = TRUE, all.inside = TRUE)
  dplyr::bind_rows(lapply(sort(unique(k)), function(j) {
    s <- k == j
    q <- stats::quantile(d[s], c(0.25, 0.5, 0.75), names = FALSE)
    dplyr::tibble(bin = j, a_lo = min(a[s]), a_hi = max(a[s]), gm_m = exp(stats::median(a[s])), n = sum(s),
                  median_d = q[2], q25_d = q[1], q75_d = q[3], frac_within_factor2 = mean(abs(d[s]) <= log(2)))
  }))
}

#' Scatter-only null for the censored tails of a measured MLH
#'
#' FI-Hyy (G 2026-10-29): would scatter alone, without range compression, explain how often the source is above
#' the measurement when the layer is too shallow to measure, and below it when the measurement is only a lower
#' bound? Under the null, `mod = h × exp(ε)` with `h` the true depth and `ε` resampled from the observed exact-value
#' log ratios `ln(mod / ref)` of the same season and regime (pools smaller than `opts$min_pool` fall back to the
#' regime, then to all exact pairs). The true depth of a censored value is unknown, so it is set to the bound:
#' - censored below `z_min` (true `h <= z_min`): simulated share of `mod <= z_min` = share with `ε <= 0`, a LOWER
#'   bound of the share expected under the null (a shallower `h` would give more);
#' - lower bound `L` (qc 3, true `h >= L`): simulated share of `mod < L` = share with `ε < 0`, an UPPER bound.
#' An observed share below the first bound (source too deep for shallow layers), or above the second (too shallow
#' for deep layers), is not explained by scatter. Intervals: simulated share, percentiles over `opts$n_sim`
#' draws; observed share, day-block bootstrap over `opts$n_boot` replicates.
#'
#' @param df Tibble, one row per pair (exact and censored): `timestamp_utc` [POSIXct UTC], `mod` [m], `kind`
#'   (`exact`, `below`, `lower_bound`), `ref` [m, exact value or lower bound], `z_min` [m], `season`, `regime`.
#' @param opts [measured_compare_options()].
#' @return Tibble per tail (`below`, `lower_bound`): `tail`, `n` [-], `n_days` [-], `obs_share`, `obs_lo`,
#'   `obs_hi` [0-1], `sim_share`, `sim_lo`, `sim_hi` [0-1], `sim_bound` (`lower`/`upper`), `n_exact_pool` [-],
#'   `n_fallback_regime`, `n_fallback_all` [-, censored rows whose pool fell back].
censored_tail_null <- function(df, opts = measured_compare_options()) {
  df <- df[!is.na(df$mod) & df$mod > 0 & df$kind %in% c("exact", "below", "lower_bound"), ]
  ex <- df[df$kind == "exact" & df$ref > 0, ]
  eps <- log(ex$mod / ex$ref)
  key_sr <- paste(ex$season, ex$regime)
  a <- (1 - opts$conf) / 2
  tails <- list(below = list(bound = function(x) x$z_min, hit = function(m, b) m <= b, sim = "lower"),
                lower_bound = list(bound = function(x) x$ref, hit = function(m, b) m < b, sim = "upper"))
  dplyr::bind_rows(lapply(names(tails), function(tl) {
    tt <- tails[[tl]]
    x <- df[df$kind == tl, ]
    x <- x[!is.na(tt$bound(x)) & tt$bound(x) > 0, ]
    res <- dplyr::tibble(tail = tl, n = nrow(x), n_days = 0L, obs_share = NA_real_, obs_lo = NA_real_,
                         obs_hi = NA_real_, sim_share = NA_real_, sim_lo = NA_real_, sim_hi = NA_real_,
                         sim_bound = tt$sim, n_exact_pool = length(eps), n_fallback_regime = 0L, n_fallback_all = 0L)
    if (nrow(x) == 0 || length(eps) == 0) return(res)
    b <- tt$bound(x)
    hit <- tt$hit(x$mod, b)
    # pool of each censored row: season x regime, else regime, else all
    pools <- split(seq_along(eps), key_sr)
    pools_r <- split(seq_along(eps), ex$regime)
    pool_of <- lapply(seq_len(nrow(x)), function(i) {
      k <- paste(x$season[i], x$regime[i])
      p <- pools[[k]]
      if (length(p) >= opts$min_pool) return(list(p, 0L, paste("sr", k)))
      k <- as.character(x$regime[i])
      p <- pools_r[[k]]
      if (length(p) >= opts$min_pool) return(list(p, 1L, paste("r", k)))
      list(seq_along(eps), 2L, "all")
    })
    lev <- vapply(pool_of, `[[`, 0L, 2)
    res$n_fallback_regime <- sum(lev == 1L)
    res$n_fallback_all <- sum(lev == 2L)
    day <- floor(as.numeric(x$timestamp_utc) / 86400)
    idx <- split(seq_along(day), day)
    res$n_days <- length(idx)
    res$obs_share <- mean(hit)
    with_seed(opts$seed, {
      # simulated share: model = bound x exp(eps), eps drawn with replacement from the row's pool. The hit
      # depends only on the sign of eps (b exp(eps) <= b iff eps <= 0), so the draws of the rows sharing a pool
      # are one binomial count with p = share of hits in the pool (same distribution, no n x n_sim matrix).
      pool_id <- vapply(pool_of, `[[`, "", 3)
      hits <- numeric(opts$n_sim)
      expct <- 0
      for (g in unique(pool_id)) {
        r <- which(pool_id == g)
        p_hit <- mean(tt$hit(exp(eps[pool_of[[r[1]]][[1]]]), 1))
        hits <- hits + stats::rbinom(opts$n_sim, length(r), p_hit)
        expct <- expct + length(r) * p_hit
      }
      sh <- hits / nrow(x)
      res$sim_share <- expct / nrow(x)
      res[c("sim_lo", "sim_hi")] <- as.list(stats::quantile(sh, c(a, 1 - a), names = FALSE))
      if (opts$n_boot > 0 && length(idx) > 1) {
        bs <- vapply(seq_len(opts$n_boot), function(k) {
          mean(hit[unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)])
        }, numeric(1))
        res[c("obs_lo", "obs_hi")] <- as.list(stats::quantile(bs, c(a, 1 - a), names = FALSE))
      }
    })
    res
  }))
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

#' Median log ratio with a day-block bootstrap interval
#'
#' Spec 11 D: median of `ln(mod / ref)` with a percentile interval from a bootstrap of whole UTC days (as
#' [compression_slope()]).
#'
#' @param mod,ref Source and measured PBLH [m] (non-positive or missing pairs dropped).
#' @param t Times of the pairs [POSIXct UTC] (day blocks).
#' @param opts [measured_compare_options()] (`n_boot`, `conf`, `seed`).
#' @return One-row tibble: `n`, `n_days` [-], `median_log_ratio`, `mlr_lo`, `mlr_hi` [-, natural log].
median_logratio_boot <- function(mod, ref, t, opts = measured_compare_options()) {
  ok <- is.finite(mod) & is.finite(ref) & mod > 0 & ref > 0 & !is.na(t)
  lr <- log(mod[ok] / ref[ok])
  idx <- split(seq_along(lr), floor(as.numeric(t[ok]) / 86400))
  out <- dplyr::tibble(n = length(lr), n_days = length(idx), median_log_ratio = if (length(lr) > 0) stats::median(lr) else NA_real_,
                       mlr_lo = NA_real_, mlr_hi = NA_real_)
  if (opts$n_boot > 0 && length(idx) > 1) {
    b <- with_seed(opts$seed, vapply(seq_len(opts$n_boot), function(k) {
      stats::median(lr[unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)])
    }, numeric(1)))
    a <- (1 - opts$conf) / 2
    out[c("mlr_lo", "mlr_hi")] <- as.list(stats::quantile(b, c(a, 1 - a), names = FALSE))
  }
  out
}
