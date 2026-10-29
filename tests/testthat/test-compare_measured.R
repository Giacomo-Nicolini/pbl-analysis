test_that("measured_interval codes exact, below z_min and lower bounds", {
  m <- measured_interval(c(300, NA, 500, 400, NA, NA), c(0L, 1L, 3L, 2L, NA, 0L), 105)
  expect_equal(m$kind, c("exact", "below", "lower_bound", "exact", NA, NA))
  expect_equal(m$left, c(300, 0, 500, 400, NA, NA))
  expect_equal(m$right, c(300, 105, NA, 400, NA, NA))
})

test_that("censored_median: Turnbull median with censored values", {
  # exact only: ordinary median
  expect_equal(censored_median(c(100, 200, 300), c(100, 200, 300)), 200)
  # two values in (0, 105], three exact: mass 2/5 below 105, then 1/5 at 300, 400, 500 -> median 300
  expect_equal(censored_median(c(0, 0, 300, 400, 500), c(105, 105, 300, 400, 500)), 300)
  # the plain median of the exact values would be 400: censoring moves the median down
  expect_equal(stats::median(c(300, 400, 500)), 400)
  # three lower bounds at 300 above two exact values: the curve never reaches 0.5 below 300 -> not estimable
  expect_true(is.na(censored_median(c(100, 200, 300, 300, 300), c(100, 200, NA, NA, NA))))
  # exact + right-censored only: same as Kaplan-Meier (S = 0.8, 0.6, 0.4 at 100, 200, 400 -> median 400)
  km <- survival::survfit(survival::Surv(c(100, 200, 300, 400, 900), c(1, 1, 0, 1, 0)) ~ 1)
  expect_equal(censored_median(c(100, 200, 300, 400, 900), c(100, 200, NA, 400, NA)),
               unname(summary(km)$table["median"]))
  # mixed: 2 below 105, exact 300, 400, 500, 2 lower bounds at 600 -> S 5/7, 4/7, 3/7 -> median 400
  expect_equal(censored_median(c(0, 0, 300, 400, 500, 600, 600), c(105, 105, 300, 400, 500, NA, NA)), 400)
  # lower bounds above the median do not change it
  expect_equal(censored_median(c(100, 200, 300, 900, 900), c(100, 200, 300, NA, NA)), 300)
  expect_true(is.na(censored_median(100, 100)))
})

test_that("measured_diurnal returns censored and exact medians per group", {
  d <- dplyr::tibble(
    season = "JJA", hour_utc = c(rep(12L, 5), rep(0L, 3)),
    pblh_m = c(NA, NA, 300, 400, 500, 100, 200, 300), qc_flag = c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L)
  )
  d <- dplyr::bind_cols(d, measured_interval(d$pblh_m, d$qc_flag, 105))
  r <- measured_diurnal(d)
  expect_equal(r$hour_utc, c(0L, 12L))
  expect_equal(r$median_censored, c(200, 300))
  expect_equal(r$median_exact, c(200, 400))
  expect_equal(r$n_below, c(0L, 2L))
})

test_that("match_measured_to_hours keeps the censoring", {
  t0 <- as.POSIXct("2022-06-01 10:00:00", tz = "UTC")
  st <- t0 + 1800 * (0:7)
  meas <- dplyr::tibble(
    timestamp_start_utc = st, timestamp_end_utc = st + 1800,
    pblh_m = c(400, 600, 400, 800, NA, NA, NA, 500), qc_flag = c(0L, 0L, 0L, 3L, 1L, 1L, 1L, 0L),
    x_censor_below_m = 105
  )
  # hours at 10:30, 11:30, 12:30, 13:30: halves (1,2), (3,4), (5,6), (7,8)
  h <- match_measured_to_hours(meas, t0 + 1800 + 3600 * (0:3))
  expect_equal(h$qc_flag, c(0L, 3L, 1L, NA))
  expect_equal(h$pblh_m, c(500, 600, NA, NA))
  expect_equal(h$x_censor_below_m, rep(105, 4))
})

test_that("lower-bound and censored-below summaries", {
  df <- dplyr::tibble(mod = c(100, 300, 900, 50), bound = c(200, 200, 300, NA), regime = c("a", "a", "b", "b"))
  s <- lower_bound_summary(df)
  expect_equal(s$n, 3)
  expect_equal(s$frac_mod_below_bound, 1 / 3)
  expect_equal(s$median_ratio, 1.5)
  sb <- lower_bound_summary(df, "regime")
  expect_equal(sb$frac_mod_below_bound, c(0.5, 0))
  c0 <- censored_below_summary(dplyr::tibble(mod = c(50, 105, 300), z_min = 105))
  expect_equal(c0$frac_mod_below_zmin, 2 / 3)
  expect_equal(c0$median_mod, 105)
})

# --- spec 10: hourly interval-mean references ---
hh_std <- function() {
  s <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + (0:5) * 1800
  dplyr::tibble(timestamp_start_utc = s, timestamp_end_utc = s + 1800, pblh_m = c(100, 300, 500, NA, 700, 900),
                qc_flag = 0L, x_regime = c("convective", "convective", "convective", "convective", "stable",
                                           "convective"))
}

test_that("match_halfhours_to_intervals: both half-hours inside the hour, regime of both", {
  h <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + (0:2) * 3600
  m <- match_halfhours_to_intervals(hh_std(), h)
  expect_equal(m$timestamp_utc, h)
  expect_equal(m$pblh_m, c(200, NA, 800))           # 10:00-11:00 = mean(100, 300); 11:00 needs 11:30 (NA)
  expect_equal(m$x_regime, c("convective", NA, "mixed"))
  expect_equal(m$n_samples, c(2L, 1L, 2L))
})

test_that("match_instants_to_intervals: mean of the values at t and t + 1 h, both required", {
  t <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + (0:3) * 3600
  e <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t, pblh_m = c(400, 600, NA, 1000))
  m <- match_instants_to_intervals(e, t[1:3])
  expect_equal(m$pblh_m, c(500, NA, NA))
  expect_equal(m$n_samples, c(2L, 1L, 1L))
})

test_that("measured_values: uncensored sources are exact where pblh_m is present", {
  s <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + (0:2) * 3600
  meas <- dplyr::tibble(timestamp_start_utc = s, timestamp_end_utc = s + 3600, pblh_m = c(500, NA, 800),
                        qc_flag = c(0L, NA, 1L))
  v <- measured_values(meas)
  expect_equal(v$kind, c("exact", NA, "exact"))
  expect_equal(v$right, c(500, NA, 800))
  # censored source: coded by measured_interval()
  meas$x_censor_below_m <- 105
  meas$qc_flag <- c(0L, 1L, 3L)
  expect_equal(measured_values(meas)$kind, c("exact", "below", "lower_bound"))
})

test_that("measured_time_base: hourly reference puts PBLHmod and ERA5 on its intervals", {
  s <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + (0:2) * 3600
  meas <- dplyr::tibble(timestamp_start_utc = s, timestamp_end_utc = s + 3600, pblh_m = c(500, 600, 800),
                        qc_flag = 0L)
  t <- as.POSIXct("2022-06-01 10:00", tz = "UTC") + (0:3) * 3600
  era <- list(era5_blh_center = dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t,
                                              pblh_m = c(400, 600, 800, 1000)))
  tb <- measured_time_base(meas, hh_std(), era)
  expect_equal(tb$base, "hour")
  expect_equal(tb$mid_s, 1800)
  expect_equal(tb$mod$pblh_m, c(200, NA, 800))
  expect_equal(tb$era$era5_blh_center$pblh_m, c(500, 700, 900))
  expect_equal(tb$era$era5_blh_center$x_regime, c("convective", NA, "mixed"))
  expect_equal(tb$era_meas$pblh_m, meas$pblh_m)
  meas$timestamp_end_utc <- meas$timestamp_start_utc + 600
  expect_error(measured_time_base(meas, hh_std(), era), "not supported")
})

test_that("keep_regime_by_distance follows the distance rule", {
  r <- c("convective", "stable", NA, "mixed")
  expect_equal(keep_regime_by_distance(r, "all_regimes"), rep(TRUE, 4))
  expect_equal(keep_regime_by_distance(r, "daytime_convective"), c(TRUE, FALSE, FALSE, FALSE))
})
