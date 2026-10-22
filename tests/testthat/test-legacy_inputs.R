utc <- function(x) as.POSIXct(x, tz = "UTC")

test_that("eddypro_ambient: hand-computed case (EddyPro engine, no biomet)", {
  # Ts = 293.15 K, chi = 10 mmol mol-1, altitude 16 m:
  # AP = 101300 exp(-16 / 8200) = 101102.5342 Pa
  # Va = 8.314 * 293.15 / AP = 0.02410671 m3 mol-1; rho_w = (10 / Va) 18.02e-3 1e-3 = 0.00747510 kg m-3
  # pass 1 at Ts: e = rho_w 461.5 Ts = 1011.30 Pa; rho_d = (AP - e) / (287.04 Ts) = 1.189498; Q = 0.00624500
  #   Ta1 = Ts / (1 + 0.51 Q) = 292.219296 K
  # iteration at Ta1: es = 2200.2137 Pa, e = 1008.0858 Pa, RH = 45.8176 %; rho_d = 1.193325; rho_a = 1.200800
  #   Cpd(Ta1) = 1005.529113; Q = 0.00622510; Ta = 292.222253 K
  # Cpv(RH, Ta) = 1874.727833; cp = (Cpv rho_w + Cpd rho_d) / rho_a = 1010.939961 J kg-1 K-1
  a <- eddypro_ambient(293.15, 10, 16)
  expect_equal(a$AP, 101102.5342, tolerance = 1e-9)
  expect_equal(a$AT, 292.222253, tolerance = 1e-8)
  expect_equal(a$rho, 1.200800, tolerance = 1e-6)
  expect_equal(a$cp, 1010.939961, tolerance = 1e-8)
})

test_that("eddypro_ambient: no H2O -> dry-air fallback; missing Ts -> NA; RH quirk", {
  a <- eddypro_ambient(c(293.15, 293.15, NA), c(NA, 0, 10), 16)
  # AT = Ts; rho = AP / (287.04 Ts) = 1.201516; cp = 1005 + (20 + 23.12)^2 / 3364 = 1005.552715
  expect_equal(a$AT[1:2], c(293.15, 293.15))
  expect_equal(a$rho[1], 101300 * exp(-16 / 8200) / (287.04 * 293.15))
  expect_equal(a$cp[1], 1005 + (293.15 - 273.15 + 23.12)^2 / 3364)
  expect_true(is.na(a$AT[3]) && is.na(a$rho[3]) && is.na(a$cp[3]))
  expect_equal(a$AP[3], 101300 * exp(-16 / 8200))
  # supersaturated (RH ~ 650 % at 0 degC, chi 40 mmol mol-1): the engine passes its error value -9999 into Cpv
  # (Cpv ~ 1859 - 0.13 * 9999 ~ 560 instead of ~ 1860), so cp is ~ 998.6 instead of ~ 1010 J kg-1 K-1
  s <- eddypro_ambient(273.15, 40, 0)
  expect_equal(s$cp, 998.619705, tolerance = 1e-8)
  expect_lt(s$cp, 1005)
})

test_that("locf: last observation carried forward, leading NA kept", {
  expect_equal(locf(c(NA, 1, NA, NA, 3, NA)), c(NA, 1, 1, 1, 3, 3))
  expect_equal(locf(c(NA_real_, NA_real_)), c(NA_real_, NA_real_))
  expect_equal(locf(c(2, NA)), c(2, 2))
})

test_that("legacy_inputs: runs split, trimmed to first data row, LOCF within each run", {
  t <- utc("2025-12-31 20:00:00") + 1800 * 0:11 # site UTC+1: local 21:00 .. 02:30 on 2026-01-01
  loc_s <- t + 3600
  d <- dplyr::tibble(
    timestamp_start_utc = t, timestamp_end_utc = t + 1800,
    TIMESTAMP_START = format(loc_s, "%Y%m%d%H%M"), TIMESTAMP_END = format(loc_s + 1800, "%Y%m%d%H%M"),
    MO_LENGTH = c(NA, NA, -50, NA, NA, 30, NA, NA, NA, 40, NA, NA),
    USTAR = c(NA, NA, 0.3, 0.3, 0.3, 0.3, NA, NA, 0.2, 0.2, 0.2, 0.2),
    H = c(NA, NA, 10, 10, 10, -5, NA, NA, -3, -3, -3, -3),
    T_SONIC_K = c(NA, NA, rep(280, 4), NA, NA, rep(279, 4)),
    H2O = 6, zm_m = 30, canopy_height_m = 21
  )
  site <- dplyr::tibble(site_id = "XX-Tst", altitude_m = 100)
  x <- legacy_inputs(d, site)
  # run 1: local starts before 2026-01-01 00:00 (rows 1-6, trimmed to row 3); run 2: rows 7-12, trimmed to row 9
  expect_equal(x$run_id, c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L))
  expect_equal(x$MOL, c(-50, -50, -50, 30, NA, 40, 40, 40)) # LOCF does not cross runs
  expect_equal(x$TIMESTAMP[1], format(loc_s[3] + 1800, "%Y-%m-%d %H:%M:%S"))
  expect_equal(x$d, rep(14, 8))
  expect_equal(x$AT[1], eddypro_ambient(280, 6, 100)$AT)
  one <- legacy_inputs(d, site, split_local = NULL)
  expect_equal(unique(one$run_id), 1L)
  expect_equal(nrow(one), 10)
  expect_equal(legacy_inputs(d, site, ambient = "AT_sonic")$AT[1], 280)
  expect_equal(attr(legacy_inputs(d, site, ambient = "rho_cp_legacy"), "pass_null"), c("rho", "cp"))
  expect_error(legacy_inputs(dplyr::select(d, -H2O), site), "H2O")
})

test_that("eddypro_ambient returns the engine RH; amb_na variant blanks AT, rho, cp", {
  a <- eddypro_ambient(c(293.15, 293.15, 273.15), c(10, NA, 40), 16)
  expect_equal(a$RH[1], 45.8176, tolerance = 1e-5) # hand case above
  expect_true(is.na(a$RH[2]))                       # no H2O
  expect_equal(a$RH[3], -9999)                      # supersaturated: engine error value
  t <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC") + 1800 * 0:2
  d <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t + 1800,
                     TIMESTAMP_START = format(t + 3600, "%Y%m%d%H%M"), TIMESTAMP_END = format(t + 5400, "%Y%m%d%H%M"),
                     MO_LENGTH = -50, USTAR = 0.4, H = 100, T_SONIC_K = c(293.15, 293.15, 273.15), H2O = c(10, NA, 40),
                     zm_m = 30, canopy_height_m = 20)
  x <- legacy_inputs(d, dplyr::tibble(altitude_m = 16), ambient = "amb_na")
  expect_equal(is.na(x$AT), c(FALSE, TRUE, TRUE))
  expect_equal(is.na(x$rho), c(FALSE, TRUE, TRUE))
  expect_equal(is.na(x$cp), c(FALSE, TRUE, TRUE))
})

test_that("tau_na: USTAR set to NA where L2 TAU is missing", {
  t <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC") + 1800 * 0:2
  d <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t + 1800,
                     TIMESTAMP_START = format(t + 3600, "%Y%m%d%H%M"), TIMESTAMP_END = format(t + 5400, "%Y%m%d%H%M"),
                     MO_LENGTH = -50, USTAR = c(0.4, 0.5, 0.6), TAU = c(-0.2, NA, -0.4), H = 100, T_SONIC_K = 293,
                     H2O = 10, zm_m = 30, canopy_height_m = 20)
  site <- dplyr::tibble(altitude_m = 16)
  expect_equal(legacy_inputs(d, site, tau_na = TRUE)$ustar, c(0.4, NA, 0.6))
  expect_equal(legacy_inputs(d, site)$ustar, c(0.4, 0.5, 0.6))
  expect_error(legacy_inputs(dplyr::select(d, -TAU), site, tau_na = TRUE), "TAU")
})

test_that("legacy_inputs_cleanset: cleanset columns, LOCF, trimmed to first data row", {
  cs <- data.frame(TIMESTAMP_START = c("202010260000", "202010260030", "202010260100", "202010260130"),
                   TIMESTAMP_END = c("202010260030", "202010260100", "202010260130", "202010260200"),
                   MO_LENGTH = c(NA, 50, NA, -20), USTAR = c(NA, 0.3, 0.3, 0.4), H = c(NA, -10, NA, 5),
                   AT = c(NA, 285, 285, 286), AP = c(NA, 101103, 101103, 101103), RHO = c(NA, 1.2, 1.2, 1.2),
                   CP = c(NA, 1010, 1010, 1010))
  x <- legacy_inputs_cleanset(cs, utc_offset_h = 1)
  expect_equal(nrow(x), 3)
  expect_equal(x$MOL, c(50, 50, -20))
  expect_equal(x$TIMESTAMP[1], "2020-10-26 01:00:00")
  expect_equal(x$timestamp_start_utc[1], as.POSIXct("2020-10-25 23:30:00", tz = "UTC"))
  expect_equal(x$rho, c(1.2, 1.2, 1.2))
  expect_error(legacy_inputs_cleanset(cs[, -9], 1), "CP")
})
