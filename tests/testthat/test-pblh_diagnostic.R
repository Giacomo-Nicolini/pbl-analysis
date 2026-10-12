test_that("h_nieuwstadt: hand-computed case (Eq. B1)", {
  # L = 100 m, u* = 0.3 m/s, f = 1e-4 s-1:
  # x = u*/(f L) = 0.3 / 0.01 = 30; sqrt(1 + 2.28 * 30) = sqrt(69.4) = 8.330666
  # h = 100 / 3.8 * (8.330666 - 1) = 26.315789 * 7.330666 = 192.9123 m
  expect_equal(h_nieuwstadt(100, 0.3, 1e-4), 192.9123, tolerance = 1e-6)
  expect_equal(h_nieuwstadt(100, 0.3, -1e-4), 192.9123, tolerance = 1e-6) # |f|
})

test_that("h_nieuwstadt: neutral and strongly stable limits", {
  # large L -> B2. Relative difference ~ b x / 4 with x = u*/(|f| L) (series of sqrt(1 + b x)):
  # u* = 0.1, f = 1e-4, L = 1e6: x = 1e-3 -> 2.28e-3 / 4 = 5.7e-4 < 1e-3
  expect_lt(abs(h_nieuwstadt(1e6, 0.1, 1e-4) / h_neutral(0.1, 1e-4) - 1), 1e-3)
  expect_equal(h_nieuwstadt(1e6, 0.1, 1e-4) / h_neutral(0.1, 1e-4) - 1, -2.28e-3 / 4, tolerance = 1e-2)
  # neutral value: 0.3 * 0.3 / 1e-4 = 900 m
  expect_equal(h_nieuwstadt(Inf, 0.3, 1e-4), 900)
  expect_equal(h_nieuwstadt(1e12, 0.3, 1e-4), 900, tolerance = 1e-6) # switch at |x| < 1e-6
  # small L: h / sqrt(u* L / |f|) -> sqrt(2.28) / 3.8 = 0.39736
  L <- 1e-3
  ratio <- h_nieuwstadt(L, 0.3, 1e-4) / sqrt(0.3 * L / 1e-4)
  expect_equal(ratio, sqrt(2.28) / 3.8, tolerance = 1e-3)
  expect_equal(sqrt(2.28) / 3.8, 0.397, tolerance = 1e-3)
})

test_that("h_nieuwstadt: NA for L <= 0 or missing", {
  expect_true(all(is.na(h_nieuwstadt(c(0, -50, NA), 0.3, 1e-4))))
  expect_equal(is.na(h_nieuwstadt(c(100, -50), c(0.3, 0.3), 1e-4)), c(FALSE, TRUE))
})

test_that("h_neutral and h_zilitinkevich: hand-computed cases", {
  # B2: 0.3 * 0.5 / 1.2e-4 = 1250 m
  expect_equal(h_neutral(0.5, 1.2e-4), 1250)
  # B4: 0.4 * sqrt(0.3 * 100 / 1e-4) = 0.4 * sqrt(3e5) = 0.4 * 547.7226 = 219.0890 m
  expect_equal(h_zilitinkevich(100, 0.3, 1e-4), 219.0890, tolerance = 1e-6)
  # c_s = 0.63: 0.63 * 547.7226 = 345.0652 m
  expect_equal(h_zilitinkevich(100, 0.3, 1e-4, c_s = 0.63), 345.0652, tolerance = 1e-6)
  expect_true(all(is.na(h_zilitinkevich(c(0, -10, NA, Inf), 0.3, 1e-4))))
})

test_that("h_diagnostic selects the formula by stability and option", {
  L <- c(100, -100, NA, 0, 100)
  u <- c(0.3, 0.3, 0.3, 0.3, NA)
  d <- h_diagnostic(L, u, 1e-4)
  expect_equal(d$method, c("B1", "B2", NA, NA, NA))
  expect_equal(d$h_diag_m[1:2], c(h_nieuwstadt(100, 0.3, 1e-4), 900))
  expect_true(all(is.na(d$h_diag_m[3:5])))

  d <- h_diagnostic(L, u, 1e-4, pblh_options(stable_method = "zilitinkevich", c_s = 0.63))
  expect_equal(d$method[1:2], c("B4", "B2"))
  expect_equal(d$h_diag_m[1], h_zilitinkevich(100, 0.3, 1e-4, 0.63))

  d <- h_diagnostic(L, u, 1e-4, pblh_options(stable_method = "neutral"))
  expect_equal(d$method[1:2], c("B2", "B2"))
  expect_equal(d$h_diag_m[1:2], c(900, 900))
})
