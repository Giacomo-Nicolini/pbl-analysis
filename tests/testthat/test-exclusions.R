utc <- function(x) as.POSIXct(x, tz = "UTC")

ex <- read_exclusions(fixture_path("exclusion_periods.csv"))

test_that("read_exclusions parses dates", {
  expect_s3_class(ex$date_start, "Date")
  expect_equal(ex$scope, "all")
  expect_no_error(read_exclusions(here::here("config/exclusion_periods.csv")))
})

test_that("SE-Nor-like period: boundaries in local standard time (UTC+1)", {
  # 2022-10-01 00:00 local = 2022-09-30 23:00 UTC (first excluded half-hour)
  # end of 2023-04-20 local = 2023-04-21 00:00 local = 2023-04-20 23:00 UTC (exclusive)
  d <- data.frame(timestamp_start_utc = utc(c(
    "2022-09-30 22:30:00", # just before: not flagged
    "2022-09-30 23:00:00", # first half-hour: flagged
    "2023-01-15 12:00:00", # inside
    "2023-04-20 22:30:00", # last half-hour: flagged
    "2023-04-20 23:00:00"  # just after: not flagged
  )))
  expect_message(
    out <- flag_exclusions(d, "XX-Nor", ex, utc_offset_h = 1),
    "3 of 5"
  )
  expect_equal(out$excl_flag, c(FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(out$excl_reason[2], "clearcut (synthetic)")
  expect_true(is.na(out$excl_reason[1]))
  expect_equal(nrow(out), 5) # flagged, not dropped

  expect_message(kept <- drop_excluded(out), "Dropped 3")
  expect_equal(nrow(kept), 2)
})

test_that("zero flagged records are reported too", {
  d <- data.frame(timestamp_start_utc = utc("2022-01-01 00:00:00"))
  expect_message(out <- flag_exclusions(d, "XX-Tst", ex, 1), "0 of 1")
  expect_false(out$excl_flag)
})

test_that("unknown scope aborts", {
  ex2 <- ex
  ex2$scope <- "fluxes_only"
  d <- data.frame(timestamp_start_utc = utc("2022-01-01 00:00:00"))
  expect_error(flag_exclusions(d, "XX-Nor", ex2, 1), "Unknown")
})
