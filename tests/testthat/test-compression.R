brute_inversions <- function(q) {
  n <- length(q)
  s <- 0
  for (i in seq_len(n - 1)) s <- s + sum(q[i] > q[(i + 1):n])
  s
}
brute_theil_sen <- function(x, y) {
  ij <- utils::combn(length(x), 2)
  dx <- x[ij[2, ]] - x[ij[1, ]]
  ok <- dx != 0
  stats::median((y[ij[2, ]] - y[ij[1, ]])[ok] / dx[ok])
}

test_that("count_inversions equals the brute-force count, ties not counted", {
  set.seed(1)
  for (n in c(1, 2, 3, 7, 16, 33, 100)) {
    q <- round(stats::rnorm(n), 1) # with ties
    expect_equal(count_inversions(q), brute_inversions(q))
  }
  expect_equal(count_inversions(5:1), 10)
  expect_equal(count_inversions(c(2, 2, 2)), 0)
})

test_that("theil_sen equals the median of pairwise slopes (odd/even pairs, ties in x)", {
  set.seed(2)
  for (n in c(2, 5, 8, 40, 101)) {
    x <- round(stats::runif(n, 0, 10), 1)
    y <- 0.7 * x + stats::rnorm(n)
    expect_equal(theil_sen(x, y)$slope, brute_theil_sen(x, y), tolerance = 1e-6)
  }
  # exact line: slope and intercept recovered; outliers do not move the median slope much
  x <- 1:50
  y <- 2 + 0.5 * x
  y[c(5, 40)] <- c(100, -100)
  r <- theil_sen(x, y)
  expect_equal(r$slope, 0.5, tolerance = 1e-6)
  expect_equal(r$intercept, 2, tolerance = 1e-5)
  expect_true(is.na(theil_sen(c(1, 1, 1), c(1, 2, 3))$slope))
})

test_that("rma_slope: symmetric, sign of r, line through the means", {
  set.seed(4)
  x <- stats::rnorm(200)
  y <- 2 + 0.5 * x + stats::rnorm(200, 0, 0.3)
  r <- rma_slope(x, y)
  expect_equal(r$slope, stats::sd(y) / stats::sd(x))
  expect_equal(rma_slope(y, x)$slope, 1 / r$slope)
  expect_equal(r$intercept, mean(y) - r$slope * mean(x))
  expect_equal(rma_slope(x, -y)$slope, -r$slope)
  expect_true(is.na(rma_slope(c(1, 2), c(1, 2))$slope))
  expect_true(is.na(rma_slope(rep(1, 5), 1:5)$slope))
})

test_that("bland_altman_fit equals least squares of the difference on the mean", {
  set.seed(5)
  x <- stats::rnorm(100, 6, 0.6)
  y <- 1 + 0.8 * x + stats::rnorm(100, 0, 0.2)
  f <- bland_altman_fit(x, y)
  co <- stats::coef(stats::lm(I(y - x) ~ I((x + y) / 2)))
  expect_equal(c(f$intercept, f$slope), unname(co))
  expect_equal(f$mean_d, mean(y - x))
  expect_equal(f$slope, 2 * (stats::var(y) - stats::var(x)) / (stats::var(x) + stats::var(y) + 2 * stats::cov(x, y)))
})

test_that("compression_slope: known compression found; scatter alone gives symmetric slope ~1, BA slope ~0", {
  set.seed(3)
  n_days <- 60
  t <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + rep(0:(n_days - 1), each = 6) * 86400 + rep(0:5, n_days) * 3600
  ref <- exp(stats::runif(length(t), log(200), log(2000)))
  mod <- exp(1.5 + 0.7 * log(ref) + stats::rnorm(length(t), 0, 0.05))
  opts <- measured_compare_options(n_boot = 200)
  cs <- compression_slope(mod, ref, t, opts)
  expect_equal(cs$n, 360)
  expect_equal(cs$n_days, 60)
  expect_equal(cs$slope_sym, 0.7, tolerance = 0.03)
  expect_true(cs$slope_sym_lo < cs$slope_sym && cs$slope_sym < cs$slope_sym_hi && cs$slope_sym_hi < 1)
  expect_true(cs$ba_slope_hi < 0)
  expect_equal(cs$ts_slope, 0.7, tolerance = 0.03)
  expect_true(cs$ts_slope_reverse > 1) # compression: reverse slope ~ 1 / 0.7
  expect_equal(cs$n_boot, 200L)
  band <- attr(cs, "ba_band")
  expect_equal(nrow(band), 40)
  expect_true(all(band$lo <= band$fit & band$fit <= band$hi))
  expect_identical(compression_slope(mod, ref, t, opts), cs)
  # regression dilution: same true depth, equal scatter in both -> forward and reverse Theil-Sen < 1,
  # symmetric slope ~ 1, Bland-Altman slope ~ 0
  h <- stats::runif(length(t), log(200), log(2000))
  m2 <- exp(h + stats::rnorm(length(t), 0, 0.4))
  r2 <- exp(h + stats::rnorm(length(t), 0, 0.4))
  c2 <- compression_slope(m2, r2, t, opts)
  expect_true(c2$ts_slope < 0.85 && c2$ts_slope_reverse < 0.85)
  expect_true(c2$slope_sym_lo < 1 && c2$slope_sym_hi > 1)
  expect_true(c2$ba_slope_lo < 0 && c2$ba_slope_hi > 0)
  # the global random state is left untouched
  set.seed(99)
  a <- stats::runif(1)
  set.seed(99)
  invisible(compression_slope(mod, ref, t, opts))
  expect_equal(stats::runif(1), a)
  # non-positive values dropped; too few points -> NA
  expect_equal(compression_slope(c(mod[1:5], 0), c(ref[1:5], 100), t[1:6], opts)$n, 5)
  expect_true(is.na(compression_slope(mod[1:2], ref[1:2], t[1:2], opts)$slope_sym))
})

test_that("bland_altman_bins: equal-count bins of the log mean", {
  ref <- c(100, 200, 400, 800, 1600, 100)
  mod <- c(200, 200, 200, 800, 800, NA)
  b <- bland_altman_bins(mod, ref, n_bins = 5)
  expect_equal(sum(b$n), 5)
  expect_equal(b$median_d, c(log(2), 0, log(0.5), 0, log(0.5)))
  expect_equal(b$gm_m, c(sqrt(2e4), 200, sqrt(8e4), 800, sqrt(1280000)))
  expect_equal(nrow(bland_altman_bins(1, 1)), 0)
})

test_that("censored_tail_null: bounds from the sign of the exact-value log ratios", {
  t0 <- as.POSIXct("2022-06-01", tz = "UTC")
  n_ex <- 100
  ex <- dplyr::tibble(timestamp_utc = t0 + seq_len(n_ex) * 3600, kind = "exact", ref = 500,
                      mod = 500 * exp(c(rep(-0.5, 30), rep(0.5, 70))), z_min = 105, season = "JJA",
                      regime = "stable")
  below <- dplyr::tibble(timestamp_utc = t0 + 1e6 + seq_len(40) * 3600, kind = "below", ref = NA_real_,
                         mod = c(rep(90, 10), rep(300, 30)), z_min = 105, season = "JJA", regime = "stable")
  lb <- dplyr::tibble(timestamp_utc = t0 + 2e6 + seq_len(20) * 7200, kind = "lower_bound", ref = 300,
                      mod = c(rep(200, 15), rep(400, 5)), z_min = 105, season = "MAM", regime = "convective")
  opts <- measured_compare_options(n_boot = 100, n_sim = 500, min_pool = 30)
  r <- censored_tail_null(dplyr::bind_rows(ex, below, lb), opts)
  expect_equal(r$tail, c("below", "lower_bound"))
  expect_equal(r$n, c(40L, 20L))
  expect_equal(r$obs_share, c(10 / 40, 15 / 20))
  # pool: 30 % of eps <= 0; MAM x convective pool empty and regime empty -> all exact pairs
  expect_equal(r$sim_share, c(0.3, 0.3))
  expect_equal(r$n_fallback_all, c(0L, 20L))
  expect_true(all(r$sim_lo <= r$sim_share & r$sim_share <= r$sim_hi))
  expect_true(all(r$obs_lo <= r$obs_share & r$obs_share <= r$obs_hi))
  expect_equal(r$sim_bound, c("lower", "upper"))
  expect_identical(censored_tail_null(dplyr::bind_rows(ex, below, lb), opts), r)
})

test_that("logratio_by_tercile splits the measured range in three", {
  ref <- c(100, 200, 300, 400, 500, 600, 700, 800, 900)
  mod <- c(200, 400, 600, 400, 500, 600, 350, 400, 450)
  lt <- logratio_by_tercile(mod, ref)
  expect_equal(lt$tercile, c("T1", "T2", "T3"))
  expect_equal(lt$n, c(3L, 3L, 3L))
  expect_equal(lt$ref_hi, c(300, 600, 900))
  expect_equal(lt$median_log_ratio, c(log(2), 0, log(0.5)))
  expect_equal(nrow(logratio_by_tercile(1, 1)), 0)
})

test_that("median_logratio_boot: median, interval around it, reproducible, RNG untouched", {
  t <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + rep(0:29, each = 4) * 86400 + rep(0:3, 30) * 3600
  set.seed(7)
  ref <- exp(stats::runif(length(t), log(300), log(2000)))
  mod <- ref * exp(-0.5 + stats::rnorm(length(t), 0, 0.2))
  o <- measured_compare_options(n_boot = 200)
  m <- median_logratio_boot(mod, ref, t, o)
  expect_equal(m$median_log_ratio, stats::median(log(mod / ref)))
  expect_true(m$mlr_lo < m$median_log_ratio && m$median_log_ratio < m$mlr_hi)
  expect_equal(m$n_days, 30L)
  expect_identical(median_logratio_boot(mod, ref, t, o), m)
  set.seed(1)
  a <- stats::runif(1)
  set.seed(1)
  invisible(median_logratio_boot(mod, ref, t, o))
  expect_equal(stats::runif(1), a)
})
