test_that("bg_start_heights: h0 from the half-hour before sunrise, else the nearest valid within the window", {
  t <- as.POSIXct("2022-06-01 00:00", tz = "UTC") + (0:15) * 1800
  is_day <- rep(c(FALSE, TRUE), each = 4, times = 2) # two days starting at rows 5 and 13
  out <- dplyr::tibble(
    timestamp_start_utc = t, is_day = is_day,
    day_id = as.Date(rep(c("2022-05-31", "2022-06-01", "2022-06-01", "2022-06-02"), each = 4)),
    regime = c(rep("stable", 4), "convective", "convective", "unstable_diag", "unstable_diag",
               rep("stable", 4), rep("unstable_diag", 4)),
    method = c(rep("B1", 4), "B5", "B5", "B3", "B3", "B1", "B2", NA, NA, rep("B3", 4)),
    h_diag_m = c(50, 60, 70, 80, 300, 400, 500, 600, 90, 95, NA, NA, 200, 200, 200, 200),
    h_bg_m = c(rep(NA, 4), 80, 150, rep(NA, 10)),
    flag_no_h0 = FALSE, flag_bg_stopped_gap = c(rep(FALSE, 4), rep(TRUE, 4), rep(FALSE, 8))
  )
  out$day_id[9:12] <- as.Date("2022-06-02") # night before the second sunrise belongs to the second day
  b <- bg_start_heights(out, h0_search = 4)
  expect_equal(b$sunrise_utc, t[c(5, 13)])
  expect_equal(b$h0_m, c(80, 95))
  expect_equal(b$h0_method, c("B1", "B2"))
  expect_equal(b$h0_lag_halfhours, c(1L, 3L))
  expect_equal(b$n_convective, c(2L, 0L))
  expect_equal(b$first_h_bg_m, c(80, NA))
  expect_equal(b$flag_bg_stopped_gap, c(TRUE, FALSE))
  out$h_diag_m[9:12] <- NA
  expect_true(is.na(bg_start_heights(out, h0_search = 4)$h0_m[2]))
  expect_error(bg_start_heights(out[-3, ]), "30-min grid")
})

test_that("daily_clearness: ratio of sums over daylight with SW_IN present", {
  day <- as.Date(rep(c("2022-06-01", "2022-06-02"), each = 4))
  pot <- c(0, 500, 1000, 500, 0, 500, 1000, 500)
  sw <- c(0, 400, 800, NA, 5, 100, 200, 100)
  k <- daily_clearness(day, sw, pot)
  expect_equal(k$kt, c(1200 / 1500, 400 / 2000))
  expect_equal(k$coverage, c(2 / 3, 1))
  expect_equal(k$n_light, c(3L, 3L))
})
