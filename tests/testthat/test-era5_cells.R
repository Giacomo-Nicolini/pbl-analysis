utc <- function(x) as.POSIXct(x, tz = "UTC")

# Synthetic 3x3 block around a tower cell at (lat0, lon0), field f(lat, lon, t).
make_block <- function(f, lat0 = 50, lon0 = 10, times = utc("2020-01-01 00:00:00") + 3600 * 0:2) {
  g <- expand.grid(cell_row = -1:1, cell_col = -1:1, timestamp_utc = times)
  g$grid_lat <- lat0 + 0.25 * g$cell_row
  g$grid_lon <- lon0 + 0.25 * g$cell_col
  g$blh_m <- f(g$grid_lat, g$grid_lon, as.numeric(g$timestamp_utc - times[1]) / 3600)
  dplyr::as_tibble(g)
}

test_that("bilinear reproduces a field linear in lat and lon exactly; weights sum to 1", {
  f <- function(lat, lon, h) 1000 + 200 * (lat - 50) - 120 * (lon - 10) + 10 * h
  d <- make_block(f)
  for (pt in list(c(50.0732, 10.1914), c(49.8811, 9.9043), c(50.2, 9.76))) {
    b <- era5_bilinear(d, pt[1], pt[2], "blh_m")
    expect_equal(b$blh_m_bilin, f(pt[1], pt[2], 0:2), tolerance = 1e-10)
    expect_equal(sum(attr(b, "weights")$weight), 1)
  }
  w <- attr(era5_bilinear(d, 49.8811, 9.9043, "blh_m"), "weights")
  expect_setequal(paste(w$cell_row, w$cell_col), c("0 0", "-1 0", "0 -1", "-1 -1"))
})

test_that("bilinear at the cell centre returns the tower cell; missing values propagate", {
  d <- make_block(function(lat, lon, h) lat * 10 + lon + h)
  b <- era5_bilinear(d, 50, 10, "blh_m")
  expect_equal(b$blh_m_bilin, d$blh_m[d$cell_row == 0 & d$cell_col == 0])
  d$blh_m[d$cell_row == 1 & d$cell_col == 1 & d$timestamp_utc == utc("2020-01-01 01:00:00")] <- NA
  b2 <- era5_bilinear(d, 50.1, 10.1, "blh_m")
  expect_equal(is.na(b2$blh_m_bilin), c(FALSE, TRUE, FALSE))
})

test_that("spread: hand-computed sd and cv on 9 values; n9 < 9 handled", {
  v <- c(100, 200, 300, 400, 500, 600, 700, 800, 900)
  d <- make_block(function(lat, lon, h) 0, times = utc("2020-01-01 00:00:00"))
  d$blh_m <- v
  s <- era5_spread(d, "blh_m")
  # mean 500; sum of squares 2 * (400^2 + 300^2 + 200^2 + 100^2) = 600000; sd = sqrt(600000 / 8) = 273.8613
  expect_equal(s$mean9, 500)
  expect_equal(s$sd9, sqrt(600000 / 8))
  expect_equal(s$cv9, sqrt(600000 / 8) / 500)
  expect_equal(c(s$min9, s$max9, s$range9, s$n9), c(100, 900, 800, 9))

  d$blh_m[c(1, 9)] <- NA # 7 values left: 200..800, mean 500, sd = sqrt(280000 / 6)
  s2 <- era5_spread(d, "blh_m")
  expect_equal(s2$n9, 7)
  expect_equal(s2$sd9, sqrt(280000 / 6))

  d$blh_m <- NA
  expect_equal(era5_spread(d, "blh_m")$n9, 0L)
})

test_that("era5_cell_metadata on the fixture", {
  cfg <- list(era5_static_dir = normalizePath(fixture_path("era5_2026", "static")))
  sites <- read_sites(fixture_path("sites.csv"))
  m <- era5_cell_metadata("XX-Tst", cfg, sites)
  expect_equal(nrow(m), 1)
  expect_equal(c(m$cell_lat, m$cell_lon, m$cell_dist_km), c(50, 10, 0))
  expect_equal(m$orog_m, 100, tolerance = 1e-9)
  expect_equal(m$dalt_m, 0, tolerance = 1e-9) # site altitude 100 m
  # orography 100 + 10 row + 5 col: min 85, max 115; lsm 0.6..1; sdor 5 + row + col: max 7
  expect_equal(c(m$orog_min_m, m$orog_max_m), c(85, 115), tolerance = 1e-9)
  expect_equal(c(m$lsm_min, m$lsm_max, m$sdor_max_m, m$n_cells), c(0.6, 1, 7, 9))
  expect_error(era5_cell_metadata("XX-Nope", cfg, sites), "not found")
})
