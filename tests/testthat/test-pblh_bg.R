test_that("bg_rate: hand-computed case (Eq. B5)", {
  # h = 500 m, wtheta = 0.1 K m/s, u* = 0.3 m/s, L = -20 m, T = 290 K, gamma = 0.01 K/m
  # d1 = 1.4 * 500 - 2 * 2.5 * 0.4 * (-20) = 700 + 40 = 740
  # d2 = 1.2 * 500 - 2.5 * 0.4 * (-20) = 600 + 20 = 620
  # term1 = 500^2 / 740 = 337.837838
  # term2 = 8 * 0.09 * 290 / (0.01 * 9.80665 * 620) = 208.8 / 60.80123 = 3.434141
  # dh/dt = (0.1 / 0.01) / (337.837838 + 3.434141) = 10 / 341.271979 = 0.0293021 m/s
  expect_equal(bg_rate(500, 0.1, 0.3, -20, 290, 0.01), 0.0293021, tolerance = 1e-5)
  # vectorised, including gamma per value
  r <- bg_rate(c(500, 500), 0.1, 0.3, -20, 290, c(0.01, 0.02))
  expect_equal(r[1], 0.0293021, tolerance = 1e-5)
  expect_lt(r[2], r[1])
})

test_that("bg_rate: NA when a denominator is <= 0", {
  # L = 100, h = 100: d1 = 140 - 2 * 2.5 * 0.4 * 100 = 140 - 200 = -60 <= 0
  expect_true(is.na(bg_rate(100, 0.1, 0.3, 100, 290, 0.01)))
  # d1 = 0 exactly: h = 2 * 2.5 * 0.4 * L / 1.4 with L = 70 -> h = 100
  expect_true(is.na(bg_rate(100, 0.1, 0.3, 70, 290, 0.01)))
  expect_true(is.na(bg_rate(NA, 0.1, 0.3, -20, 290, 0.01)))
})

test_that("bg_step reproduces the analytic solution (u* = 0, L = 0)", {
  # dh/dt = (1 + 2A) wtheta / (gamma h)  ->  h(t)^2 = h0^2 + 2 (1 + 2A) wtheta t / gamma
  # h0 = 100, wtheta = 0.1, gamma = 0.01, t = 1800: h^2 = 1e4 + 2.8 * 0.1 * 1800 / 0.01 = 60400
  # h(1800) = sqrt(60400) = 245.764115 m
  # mean over the interval: 2 / (3 a T) ((h0^2 + a T)^1.5 - h0^3), a = 28 -> 183.123711 m
  s <- bg_step(100, 0.1, 0, 0, 290, 0.01)
  expect_true(s$ok)
  expect_lt(abs(s$h_end / 245.764115 - 1), 1e-4)
  expect_lt(abs(s$h_mean / 183.123711 - 1), 1e-3)
})

test_that("bg_step converges: n_sub = 30 vs 300 within 0.5 %", {
  s30 <- bg_step(100, 0.1, 0.3, -20, 290, 0.01, n_sub = 30)
  s300 <- bg_step(100, 0.1, 0.3, -20, 290, 0.01, n_sub = 300)
  expect_lt(abs(s30$h_end / s300$h_end - 1), 0.005)
  expect_lt(abs(s30$h_mean / s300$h_mean - 1), 0.005)
  expect_gt(s30$h_end, 100)
  expect_gt(s30$h_end, s30$h_mean)
})

test_that("bg_step holds h and reports when the rate is undefined", {
  s <- bg_step(100, 0.1, 0.3, 100, 290, 0.01) # L > 0, d1 < 0
  expect_false(s$ok)
  expect_equal(s$h_end, 100)
  expect_equal(s$h_mean, 100)
  # zero heat flux: no growth, still ok
  s <- bg_step(300, 0, 0.3, -1e6, 290, 0.01)
  expect_true(s$ok)
  expect_equal(s$h_end, 300)
})
