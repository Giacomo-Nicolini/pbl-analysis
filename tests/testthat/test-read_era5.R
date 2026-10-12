utc <- function(x) as.POSIXct(x, tz = "UTC")

era_root <- normalizePath(fixture_path("era5"))
cfg <- list(era5_dir = file.path(era_root, "{site}"))

find_q <- function() suppressWarnings(suppressMessages(find_era5_files("XX-Tst", cfg)))

test_that("file names are parsed, non-BLH files ignored, ERA5 preferred over ERA5T", {
  expect_message(
    expect_warning(f <- find_era5_files("XX-Tst", cfg), "keeping ERA5"),
    "ignoring 2 files"
  )
  expect_equal(nrow(f), 3)
  expect_false(any(grepl("2m_temperature|README", f$file)))
  c2020 <- f[f$year == 2020 & f$lon_offset == 0 & f$lat_offset == 0, ]
  expect_equal(nrow(c2020), 1)
  expect_equal(c2020$product, "ERA5")
  expect_equal(f$product[f$year == 2021], "ERA5T")
  nb <- f[f$lon_offset != 0, ]
  expect_equal(nb$lon_offset, 1L)
  expect_equal(nb$lat_offset, -1L)
  expect_type(f$year, "integer")
})

test_that("offsets selection", {
  f <- find_q()
  d <- suppressWarnings(read_era5_blh(f, offsets = "center"))
  expect_equal(nrow(d), 8)
  expect_true(all(d$lon_offset == 0 & d$lat_offset == 0))
  d_all <- suppressWarnings(read_era5_blh(f, offsets = "all"))
  expect_equal(nrow(d_all), 12)
  expect_setequal(unique(d_all$lon_offset), c(0L, 1L))
})

test_that("ERA5 output columns, UTC times, and reporting of gaps / negative values", {
  f <- find_q()
  w <- testthat::capture_warnings(d <- read_era5_blh(f))
  expect_true(any(grepl("gap", w)))
  expect_true(any(grepl("1 negative", w)))
  expect_named(d, c("timestamp_utc", "blh_m", "product", "lon_offset", "lat_offset",
                    "grid_lat", "grid_lon", "site_lat", "site_lon", "grid_dist_km"))
  expect_equal(attr(d$timestamp_utc, "tzone"), "UTC")
  expect_equal(d$timestamp_utc[1], utc("2020-01-01 00:00:00"))
  expect_equal(d$blh_m[1:2], c(500, 520)) # ERA5, not the ERA5T values
  expect_equal(sum(d$blh_m < 0), 1)       # kept, not removed
  expect_equal(d$grid_dist_km[1], 0)
})

test_that("time != datetime warns; duplicated timestamps abort", {
  f <- find_q()
  f <- f[f$year == 2020 & f$lon_offset == 0 & f$lat_offset == 0, ]
  tmp <- withr::local_tempdir()
  x <- readr::read_csv(f$path, col_types = readr::cols(.default = "c"))
  x$datetime[2] <- "2020-01-01 05:00:00"
  p <- file.path(tmp, f$file)
  readr::write_csv(x, p)
  f$path <- p
  expect_warning(read_era5_blh(f), "datetime")

  x <- rbind(x, x[1, ])
  readr::write_csv(x, p)
  expect_error(suppressWarnings(read_era5_blh(f)), "duplicated")
})

test_that("haversine distance against hand-computed values", {
  # 1 degree of latitude along a meridian = 6371 * pi / 180 = 111.19493 km
  expect_equal(haversine_km(50, 10, 51, 10), 111.19493, tolerance = 1e-6)
  # 1 degree of longitude on the equator = same arc
  expect_equal(haversine_km(0, 0, 0, 1), 111.19493, tolerance = 1e-6)
  # 1 degree of longitude at 60 N: a = cos(60)^2 sin(0.5 deg)^2,
  # d = 2R asin(cos(60) sin(0.5 deg)) = 2 * 6371 * asin(0.5 * 0.0087265) = 55.59693 km
  expect_equal(haversine_km(60, 0, 60, 1), 55.59693, tolerance = 1e-6)
  expect_equal(haversine_km(51.3, 4.5, 51.3, 4.5), 0)
})
