f_ms <- function() fixture_path("meteoswiss", "meteoswiss_test.nc")

test_that("read_meteoswiss_stratfinder reads profiles, fill values and attributes", {
  r <- read_meteoswiss_stratfinder(f_ms())
  expect_equal(nrow(r), 63)
  expect_equal(r$time_utc[1], as.POSIXct("2022-06-01 00:00:30", tz = "UTC"))
  expect_equal(sum(is.na(r$mlh_m)), 2)
  expect_equal(sort(unique(r$quality_mlh)), c(0L, 1L, 9L))
  m <- attr(r, "meteoswiss_meta")
  expect_equal(m$station_alt, 1590)
  expect_equal(m$licence, "CC BY-NC 4.0")
  expect_equal(m$search_upper_limit, 4000)
  expect_error(read_meteoswiss_stratfinder("nothing.nc"), "not found")
})

test_that("meteoswiss_halfhour: quality 0 only, distinct minutes, coverage -> qc, bit 8 share", {
  h <- meteoswiss_halfhour(read_meteoswiss_stratfinder(f_ms()))
  expect_equal(h$timestamp_start_utc, as.POSIXct("2022-06-01", tz = "UTC") + c(0, 1800, 3600))
  expect_equal(h$n_samples, c(28L, 20L, 10L)) # quality-1 profiles not counted; two profiles in one minute once
  expect_equal(h$qc_flag, c(0L, 1L, NA))
  expect_equal(h$pblh_m, c(227, 500, NA)) # quality-1 values (5000 m) not used; 10 of 30 minutes: no value
  expect_equal(h$x_share_lowest, c(7 / 28, 0, 0))
  expect_equal(h$x_ablh_m, c(900, NA, 1200))
  expect_equal(h$x_segment, c(0L, 0L, 1L))
  expect_equal(h$x_n_profiles, c(30L, 21L, 12L))
})

test_that("standardize_meteoswiss_mlh writes a valid standard table", {
  r <- read_meteoswiss_stratfinder(f_ms())
  ex <- read_exclusions(fixture_path("exclusion_periods.csv"))
  s <- standardize_meteoswiss_mlh(r, "XX-Tst", ex, 1)
  expect_equal(s$data$source_id[1], "meteoswiss_stratfinder_mlh")
  expect_equal(s$meta$method_class, "aerosol")
  expect_equal(s$meta$time_reference, "interval_mean")
  expect_equal(s$meta$height_reference, "agl")
  expect_equal(s$meta$instrument_altitude_m, 1590)
  expect_match(s$meta$citation, "Kotthaus2020")
  expect_match(s$meta$licence, "CC BY-NC 4.0")
  expect_true(validate_standard_pblh(s$data, s$meta))
  dir <- tempfile("std_")
  p <- write_standard_pblh(s$data, s$meta, dir)
  expect_equal(read_standard_pblh(p[["data"]])$pblh_m, s$data$pblh_m)
})
