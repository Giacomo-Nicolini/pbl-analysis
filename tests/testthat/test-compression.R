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

test_that("compression_slope recovers a known ln-ln slope; day-block bootstrap is reproducible", {
  set.seed(3)
  n_days <- 60
  t <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + rep(0:(n_days - 1), each = 6) * 86400 + rep(0:5, n_days) * 3600
  ref <- exp(stats::runif(length(t), log(200), log(2000)))
  mod <- exp(1.5 + 0.7 * log(ref) + stats::rnorm(length(t), 0, 0.05))
  opts <- measured_compare_options(n_boot = 200)
  cs <- compression_slope(mod, ref, t, opts)
  expect_equal(cs$n, 360)
  expect_equal(cs$n_days, 60)
  expect_equal(cs$slope, 0.7, tolerance = 0.03)
  expect_true(cs$slope_lo < 0.7 && cs$slope_hi > 0.7)
  expect_true(cs$slope_lo < cs$slope && cs$slope < cs$slope_hi)
  expect_equal(cs$n_boot, 200L)
  expect_true(cs$slope_reverse > 1) # compression: reverse slope ~ 1 / 0.7
  expect_identical(compression_slope(mod, ref, t, opts), cs)
  # the global random state is left untouched
  set.seed(99)
  a <- stats::runif(1)
  set.seed(99)
  invisible(compression_slope(mod, ref, t, opts))
  expect_equal(stats::runif(1), a)
  # non-positive values dropped; too few points -> NA
  expect_equal(compression_slope(c(mod[1:5], 0), c(ref[1:5], 100), t[1:6], opts)$n, 5)
  expect_true(is.na(compression_slope(mod[1:2], ref[1:2], t[1:2], opts)$slope))
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
