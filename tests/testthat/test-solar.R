utc <- function(x) as.POSIXct(x, tz = "UTC")
mins <- function(a, b) abs(as.numeric(difftime(a, b, units = "mins")))

test_that("equation of time extremes (textbook values)", {
  # EoT maximum ~ +16.4 min around 3 Nov, minimum ~ -14.2 min around 11-12 Feb
  # (e.g. Meeus 1991, ch. 28; NOAA solar calculator)
  nov <- solar_params(utc(seq(utc("2024-10-25 12:00:00"), by = "day", length.out = 20)))$eot
  feb <- solar_params(utc(seq(utc("2024-02-01 12:00:00"), by = "day", length.out = 25)))$eot
  expect_equal(max(nov), 16.4, tolerance = 0.2 / 16.4)
  expect_equal(min(feb), -14.2, tolerance = 0.2 / 14.2)
})

test_that("solstice and equinox declination", {
  # June solstice 2024-06-20 20:51 UTC: declination = +obliquity ~ 23.44 deg
  expect_equal(solar_params(utc("2024-06-20 20:51:00"))$decl, 23.44, tolerance = 0.01 / 23.44)
  # December solstice 2024-12-21 09:20 UTC: ~ -23.44 deg
  expect_equal(solar_params(utc("2024-12-21 09:20:00"))$decl, -23.44, tolerance = 0.01 / 23.44)
  # March equinox 2024-03-20 03:06 UTC: ~ 0 deg
  expect_lt(abs(solar_params(utc("2024-03-20 03:06:00"))$decl), 0.01)
})

test_that("solar noon follows the equation of time", {
  # Greenwich, 3 Nov 2024: noon = 12:00 - EoT (~16.4 min) ~ 11:43.6 UTC
  st <- sun_times(as.Date("2024-11-03"), 51.48, 0)
  expect_lt(mins(st$solar_noon_utc, utc("2024-11-03 11:43:36")), 0.5)
  # elevation maximal at solar noon: 90 - lat + decl
  sp <- solar_params(st$solar_noon_utc)
  expect_equal(solar_elevation(st$solar_noon_utc, 51.48, 0), 90 - 51.48 + sp$decl, tolerance = 1e-4)
})

test_that("geometric sunrise/sunset: equinox day length and symmetry", {
  # Equator at equinox: geometric day length 12 h (within ~1 min)
  st <- sun_times(as.Date("2024-03-20"), 0, 0)
  expect_lt(abs(as.numeric(difftime(st$sunset_utc, st$sunrise_utc, units = "mins")) - 720), 1.5)
  # sunrise and sunset symmetric around solar noon (within 1 min)
  st <- sun_times(as.Date("2024-06-21"), 51.31, 4.52)
  d1 <- mins(st$solar_noon_utc, st$sunrise_utc)
  d2 <- mins(st$sunset_utc, st$solar_noon_utc)
  expect_lt(abs(d1 - d2), 1)
  # solar elevation is ~0 at the computed geometric sunrise/sunset
  expect_lt(abs(solar_elevation(st$sunrise_utc, 51.31, 4.52)), 0.05)
  expect_lt(abs(solar_elevation(st$sunset_utc, 51.31, 4.52)), 0.05)
  # summer solstice day length at 51.31 N: 2 * acos(-tan(lat) tan(23.44)) / 15 h = 16.27 h
  expected <- 2 * acos(-tan(51.31 * pi / 180) * tan(23.44 * pi / 180)) * 180 / pi / 15 * 60
  expect_lt(abs(d1 + d2 - expected), 3)
})

test_that("refraction zenith gives earlier sunrise and later sunset", {
  g <- sun_times(as.Date("2024-03-21"), 51.31, 4.52)
  r <- sun_times(as.Date("2024-03-21"), 51.31, 4.52, zenith_deg = 90.833)
  expect_gt(mins(g$sunrise_utc, r$sunrise_utc), 3)
  expect_lt(r$sunrise_utc, g$sunrise_utc)
  expect_gt(r$sunset_utc, g$sunset_utc)
})

test_that("polar day and night do not crash", {
  st <- sun_times(as.Date(c("2024-06-21", "2024-12-21", "2024-03-21")), 78.2, 15.6)
  expect_equal(st$polar, c("day", "night", NA))
  expect_true(all(is.na(st$sunrise_utc[1:2])))
  expect_false(is.na(st$sunrise_utc[3]))

  t <- seq(utc("2024-06-21 00:00:00"), by = 1800, length.out = 48)
  expect_true(all(is_day(t, t + 1800, 78.2, 15.6)))
  t <- seq(utc("2024-12-21 00:00:00"), by = 1800, length.out = 48)
  expect_false(any(is_day(t, t + 1800, 78.2, 15.6)))
})

test_that("is_day uses the interval midpoint", {
  st <- sun_times(as.Date("2024-03-21"), 51.31, 4.52)
  # interval whose midpoint is 5 min after sunrise is day; 5 min before is night
  s1 <- st$sunrise_utc + 5 * 60 - 900
  s0 <- st$sunrise_utc - 5 * 60 - 900
  expect_true(is_day(s1, s1 + 1800, 51.31, 4.52))
  expect_false(is_day(s0, s0 + 1800, 51.31, 4.52))
})
