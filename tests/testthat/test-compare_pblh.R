utc <- function(x) as.POSIXct(x, tz = "UTC")

# Half-hours 00:00-00:30 ... 02:30-03:00 UTC.
hh <- function(pblh, regime = NULL, qc = NULL) {
  s <- utc("2020-06-01 00:00:00") + 1800 * (seq_along(pblh) - 1)
  d <- dplyr::tibble(timestamp_start_utc = s, timestamp_end_utc = s + 1800, pblh_m = pblh)
  if (!is.null(regime)) d$x_regime <- regime
  if (!is.null(qc)) d$qc_flag <- as.integer(qc)
  d
}

test_that("hour value = mean of the half-hour ending at t and the one starting at t", {
  d <- hh(c(100, 200, 300, 400, 500, 600))
  hours <- utc("2020-06-01 00:00:00") + 3600 * 0:3
  m <- match_halfhours_to_hours(d, hours)
  # 00:00: only the half-hour starting at 00:00 exists -> NA with require_both
  # 01:00: (200 + 300) / 2; 02:00: (400 + 500) / 2; 03:00: only 600 -> NA
  expect_equal(m$pblh_m, c(NA, 250, 450, NA))
  expect_equal(m$n_samples, c(1L, 2L, 2L, 1L))
  m2 <- match_halfhours_to_hours(d, hours, require_both = FALSE)
  expect_equal(m2$pblh_m, c(100, 250, 450, 600))
})

test_that("require_both with a missing half-hour; regime and qc rules", {
  d <- hh(c(100, 200, NA, 400, 500, 600),
          regime = c("stable", "stable", "stable", "convective", "stable", "stable"),
          qc = c(0, 1, 0, 0, 3, 0))
  hours <- utc("2020-06-01 00:00:00") + 3600 * 1:2
  m <- match_halfhours_to_hours(d, hours)
  expect_equal(m$pblh_m, c(NA, 450))
  expect_equal(m$x_regime, c(NA, "mixed"))
  expect_equal(m$qc_flag, c(NA, 3L))
  m2 <- match_halfhours_to_hours(d, hours, require_both = FALSE)
  expect_equal(m2$pblh_m, c(200, 450))
  expect_equal(m2$x_regime, c("stable", "mixed")) # regime of the valid half-hour
  expect_equal(m2$qc_flag, c(1L, 3L))
  expect_equal(m2$n_samples, c(1L, 2L))
  d$x_regime[4:5] <- "convective"
  expect_equal(match_halfhours_to_hours(d, hours)$x_regime, c(NA, "convective"))
})

test_that("accumulation matching uses the half-hours ending at t - 30 min and t", {
  s <- utc("2020-06-01 00:00:00") + 1800 * 0:3
  d <- dplyr::tibble(timestamp_start_utc = s, timestamp_end_utc = s + 1800, H = c(10, 20, 30, 40))
  # interval (00:00, 01:00]: half-hours ending 00:30 (10) and 01:00 (20)
  m <- match_halfhours_to_accum(d, utc("2020-06-01 01:00:00") + 3600 * 0:1, "H")
  expect_equal(m$value, c(15, 35))
  expect_equal(m$n_samples, c(2L, 2L))
  d$H[3] <- NA
  expect_equal(match_halfhours_to_accum(d, utc("2020-06-01 02:00:00"), "H")$value, NA_real_)
  expect_equal(match_halfhours_to_accum(d, utc("2020-06-01 02:00:00"), "H", require_both = FALSE)$value, 40)
})

test_that("pblh_metrics: hand-computed values", {
  mod <- c(200, 300, 500, 800, NA)
  ref <- c(100, 300, 1000, 400, 50)
  m <- pblh_metrics(mod, ref)
  # errors: 100, 0, -500, 400 -> bias 0, mae 250, rmse sqrt((1e4 + 0 + 2.5e5 + 1.6e5) / 4) = sqrt(105000)
  expect_equal(m$n, 4L)
  expect_equal(m$bias, 0)
  expect_equal(m$mae, 250)
  expect_equal(m$rmse, sqrt(105000))
  expect_equal(m$r_pearson, stats::cor(mod[1:4], ref[1:4]))
  # ranks mod 1 2 3 4, ref 1 2 4 3 -> rho = 1 - 6 * 2 / (4 * 15) = 0.8
  expect_equal(m$r_spearman, 0.8)
  # log ratios: log 2, 0, -log 2, log 2 -> median log(2) / 2; deviations log2/2 * (1, 1, 3, 1)
  # -> median absolute deviation log(2)/2, scaled: 1.4826 * 0.693147 / 2 = 0.51383
  expect_equal(m$median_log_ratio, log(2) / 2)
  expect_equal(m$mad_log_ratio, 1.4826 * log(2) / 2)
  expect_equal(round(m$mad_log_ratio, 5), 0.51383)
  expect_equal(m$frac_within_factor2, 1)
  expect_equal(m$n_log_excluded, 0L)
})

test_that("log metrics skip non-positive values and report them", {
  m <- pblh_metrics(c(100, 0, -5, 400), c(100, 200, 300, 100))
  expect_equal(m$n, 4L)
  expect_equal(m$n_log, 2L)
  expect_equal(m$n_log_excluded, 2L)
  expect_equal(m$median_log_ratio, log(4) / 2)
  expect_equal(m$frac_within_factor2, 0.5)
  e <- pblh_metrics(numeric(0), numeric(0))
  expect_equal(e$n, 0L)
  expect_true(is.na(e$bias) && is.na(e$median_log_ratio))
})

test_that("pair_table, metrics by stratum and long format", {
  t <- utc("2020-01-01 00:00:00") + 3600 * c(0, 1, 2, 24 * 180)
  mod <- dplyr::tibble(timestamp_utc = t, pblh_m = c(100, 200, NA, 400),
                       x_regime = c("stable", "stable", "convective", "convective"), qc_flag = 0L)
  ref <- dplyr::tibble(timestamp_utc = t, pblh_m = c(50, 400, 300, 400), x_blh_cv9 = 0.1)
  pt <- pair_table(mod, ref, "p1")
  expect_equal(nrow(pt), 3)
  expect_equal(as.character(pt$season), c("DJF", "DJF", "JJA")) # 2020-06-29
  expect_equal(pt$hour_utc, c(0L, 1L, 0L))
  by <- pblh_metrics_by(pt, by = "regime")
  expect_equal(by$regime, c("convective", "stable"))
  expect_equal(by$n, c(1L, 2L))
  expect_equal(by$bias, c(0, (50 - 200) / 2))
  lg <- pblh_metrics_long(pt)
  expect_setequal(unique(lg$stratum), c("all", "regime=convective", "regime=stable", "season=DJF",
                                        "season=JJA", "hour_utc=0", "hour_utc=1"))
  expect_equal(lg$value[lg$stratum == "all" & lg$metric == "n"], 3)
  expect_equal(lg$value[lg$stratum == "regime=stable" & lg$metric == "mae"], 125)
})

test_that("diurnal_summary: median, quartiles, n", {
  df <- dplyr::tibble(season = "JJA", hour_utc = c(0L, 0L, 0L, 0L, 1L), a = c(1, 2, 3, 4, 10), b = c(NA, 1, 1, 1, 2))
  s <- diurnal_summary(df, c("a", "b"))
  r <- s[s$variable == "a" & s$hour_utc == 0, ]
  expect_equal(c(r$median, r$p25, r$p75, r$n), c(2.5, 1.75, 3.25, 4))
  expect_equal(s$n[s$variable == "b" & s$hour_utc == 0], 3)
})

test_that("EC H from model output; cloud classes; season", {
  out <- dplyr::tibble(timestamp_start_utc = utc("2020-06-01 12:00:00"), timestamp_end_utc = utc("2020-06-01 12:30:00"),
                       wtheta = 0.1, p_Pa = 100000, T_K = 300)
  # rho = 100000 / (287.05 * 300) = 1.161236 kg m-3; H = 0.1 * 1.161236 * 1005 = 116.704 W m-2
  expect_equal(ec_h_from_pblhmod(out)$H_ec_Wm2, 0.1 * 100000 / (287.05 * 300) * 1005)
  expect_equal(round(ec_h_from_pblhmod(out)$H_ec_Wm2, 3), 116.704)
  expect_equal(as.character(cloud_class(c(NA, 500, 2000, 2000.1))),
               c("no_cloud", "cbh_le_2000", "cbh_le_2000", "cbh_gt_2000"))
  expect_equal(as.character(season_of(utc(c("2020-12-31 23:00:00", "2021-03-01 00:00:00", "2021-08-31 00:00:00",
                                            "2021-11-30 00:00:00")))), c("DJF", "MAM", "JJA", "SON"))
})

test_that("gamma_850_700 from the pressure-level fixture", {
  cfg <- list(era5_pl_dir = file.path(normalizePath(fixture_path("era5_2026", "pl")), "{site}"))
  pl <- read_era5_pressure_levels(suppressMessages(find_era5_files("XX-Tst", cfg, "pressure_level")), cells = "center")
  g <- gamma_850_700(pl)
  expect_equal(nrow(g), 8)
  # first time: theta_700 = 272 (1000/700)^k, theta_850 = 280 (1000/850)^k, k = 287.05/1005; dz = 3000 - 1500 m
  k <- 287.05 / 1005
  expect_equal(g$gamma_850_700[1], (272 * (1000 / 700)^k - 280 * (1000 / 850)^k) / 1500, tolerance = 1e-9)
  gm <- gamma_monthly(g)
  expect_equal(gm$hour_utc, c(0L, 12L))
  expect_equal(gm$n, c(2L, 2L))
})

test_that("hour_bin_centred puts interval midpoints in bins centred on the full hour", {
  t0 <- as.POSIXct("2022-06-01 11:30:00", tz = "UTC")
  # half-hours starting 11:30, 12:00 (midpoints 11:45, 12:15) -> 12; 12:30 (12:45) -> 13; 23:30 (23:45) -> 0
  expect_equal(hour_bin_centred(t0 + c(0, 1800, 3600) + 900), c(12L, 12L, 13L))
  expect_equal(hour_bin_centred(as.POSIXct("2022-06-01 23:30:00", tz = "UTC") + 900), 0L)
  # instantaneous value at 12:00 -> 12
  expect_equal(hour_bin_centred(as.POSIXct("2022-06-01 12:00:00", tz = "UTC")), 12L)
})
