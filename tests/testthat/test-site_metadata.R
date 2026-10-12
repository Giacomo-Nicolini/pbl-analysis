utc <- function(x) as.POSIXct(x, tz = "UTC")

sites <- read_sites(fixture_path("sites.csv"))
sites_tv <- read_sites_timevarying(fixture_path("sites_timevarying.csv"))

test_that("read_sites parses types and validates required columns", {
  expect_type(sites$lat, "double")
  expect_type(sites$igbp, "character")
  expect_true(is.na(sites$zm_m[sites$site_id == "XX-Mis"]))

  tmp <- withr::local_tempfile(fileext = ".csv")
  x <- readLines(fixture_path("sites.csv"))
  writeLines(c(x, x[2]), tmp)
  expect_error(read_sites(tmp), "Duplicated")
  writeLines(c(x[1], sub("^XX-Tst,50.0", "XX-Tst,NA", x[2])), tmp)
  expect_error(read_sites(tmp), "XX-Tst")
})

test_that("the project metadata files are valid", {
  expect_no_error(read_sites(here::here("config/sites.csv")))
  expect_no_error(read_sites_timevarying(here::here("config/sites_timevarying.csv")))
})

test_that("read_sites_timevarying rejects unknown variables and overlapping zm", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  x <- readLines(fixture_path("sites_timevarying.csv"))
  writeLines(c(x, "XX-Tst,d_m,3,2020-01-01,NA,bad"), tmp)
  expect_error(read_sites_timevarying(tmp), "Unknown")
  writeLines(c(x, "XX-Tst,zm_m,50,2021-06-01,NA,overlap"), tmp)
  expect_error(read_sites_timevarying(tmp), "Overlapping")
})

test_that("zm is a step function in local standard time (UTC+1)", {
  # 2020-01-01 00:00 local (UTC+1) = 2019-12-31 23:00 UTC
  t <- utc(c(
    "2019-12-31 22:30:00", # before first period
    "2019-12-31 23:00:00", # first half-hour of period 1
    "2020-05-31 22:30:00", # last half-hour of period 1
    "2020-05-31 23:00:00", # first half-hour of period 2
    "2021-01-15 12:00:00", # gap between periods 2 and 3
    "2021-12-31 22:30:00", # last half-hour of closed period 3
    "2021-12-31 23:00:00"  # after closed last period
  ))
  h <- get_site_heights("XX-Tst", t, sites, sites_tv)
  expect_equal(h$zm_m, c(30, 30, 30, 35, NA, 40, NA))
  expect_equal(h$zm_source, c(
    "extrapolated_before_first", "timevarying", "timevarying", "timevarying",
    "gap", "timevarying", "gap"
  ))
  expect_equal(h$timestamp_utc, t)
})

test_that("zm falls back to static, then missing", {
  t <- utc(c("2020-01-01 00:00:00", "2024-07-01 12:00:00"))
  h <- get_site_heights("XX-Sta", t, sites, sites_tv)
  expect_equal(h$zm_m, c(3, 3))
  expect_equal(h$zm_source, c("static", "static"))
  h <- get_site_heights("XX-Mis", t, sites, sites_tv)
  expect_true(all(is.na(h$zm_m)))
  expect_equal(h$zm_source, c("missing", "missing"))
})

test_that("canopy height is interpolated and held constant outside the range", {
  # measurements: 20 m at 2020-02-29 23:00 UTC, 22 m at 2020-03-10 23:00 UTC
  t <- utc(c(
    "2020-01-01 00:00:00", # before first measurement
    "2020-02-29 23:00:00", # at first measurement
    "2020-03-05 23:00:00", # midpoint
    "2020-03-10 23:00:00", # at last measurement
    "2021-01-01 00:00:00"  # after last measurement
  ))
  h <- get_site_heights("XX-Tst", t, sites, sites_tv)
  expect_equal(h$canopy_height_m, c(20, 20, 21, 22, 22))
  expect_equal(h$hc_source, c(
    "held_constant", "interpolated", "interpolated", "interpolated", "held_constant"
  ))
})

test_that("canopy height: single measurement, static and missing fallbacks", {
  t <- utc(c("2019-01-01 00:00:00", "2022-01-01 00:00:00"))
  h <- get_site_heights("XX-One", t, sites, sites_tv)
  expect_equal(h$canopy_height_m, c(6, 6))
  expect_equal(h$hc_source, c("held_constant", "held_constant"))
  expect_equal(h$zm_source, c("static", "static"))

  h <- get_site_heights("XX-Sta", t, sites, sites_tv)
  expect_equal(h$canopy_height_m, c(0.5, 0.5))
  expect_equal(h$hc_source, c("static", "static"))

  h <- get_site_heights("XX-Mis", t, sites, sites_tv)
  expect_true(all(is.na(h$canopy_height_m)))
  expect_equal(h$hc_source, c("missing", "missing"))
})

test_that("unknown site aborts", {
  expect_error(get_site_heights("XX-Nop", utc("2020-01-01"), sites, sites_tv), "not found")
})
