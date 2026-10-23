test_that("h_crossings finds morning and afternoon crossings by interpolation", {
  m <- seq(0, 1410, by = 30)
  h <- pmax(-20, 200 * sin(pi * (m - 360) / 720)) # zero at 06:00 and 18:00, max at 12:00
  cr <- h_crossings(m, h, c(0, 0.2))
  expect_equal(cr$morning_min[1], 360, tolerance = 1e-6)
  expect_equal(cr$afternoon_min[1], 1080, tolerance = 1e-6)
  # 20 % of max: sin = 0.2 -> asin(0.2) * 720 / pi after 06:00 (linear interpolation: within 1 min)
  expect_equal(cr$morning_min[2], 360 + asin(0.2) * 720 / pi, tolerance = 1)
  expect_equal(cr$afternoon_min[2], 1080 - asin(0.2) * 720 / pi, tolerance = 1)
})

test_that("h_best_lag recovers a shift (positive = compared cycle earlier)", {
  m <- seq(15, 1425, by = 30)
  f <- function(x) pmax(-20, 200 * sin(pi * (x - 360) / 720))
  ref <- f(m)
  mod_m <- seq(30, 1410, by = 60)
  mod <- f(mod_m + 35) # mod reaches each value 35 min earlier
  expect_equal(h_best_lag(m, ref, mod_m, mod), 35)
  expect_equal(h_best_lag(m, ref, m, ref), 0)
})

test_that("h_timing_table places values per season and gives the ERA5 lag on EC rows", {
  t_mid_ec <- as.POSIXct("2022-06-01", tz = "UTC") + 900 + 1800 * (0:(48 * 10 - 1))
  t_mid_era <- as.POSIXct("2022-06-01", tz = "UTC") + 1800 + 3600 * (0:(24 * 10 - 1))
  f <- function(t) { x <- (as.numeric(t) %% 86400) / 60; pmax(-20, 200 * sin(pi * (x - 360) / 720)) }
  tab <- h_timing_table(list(EC_H = dplyr::tibble(t_mid = t_mid_ec, h = f(t_mid_ec)),
                             ERA5_H = dplyr::tibble(t_mid = t_mid_era, h = f(t_mid_era + 1200))))
  expect_equal(tab$season, c("JJA", "JJA"))
  expect_equal(tab$lag_era5_min, c(20, NA))
  expect_equal(tab$morning_0pct[1], 360, tolerance = 1e-6)
  expect_true(tab$morning_0pct[2] < tab$morning_0pct[1])
})
