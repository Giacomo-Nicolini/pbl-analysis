test_that("theta and level heights: Poisson's equation and Phi / g (hand case)", {
  k <- pblh_constants()
  # T = 280 K at 850 hPa: exponent 287.05 / 1005 = 0.28562; ln(1000 / 850) = 0.16252; theta = 280 e^0.046419
  # = 293.30 K
  expect_equal(round(280 * (1000 / 850)^(k$r_d / k$cp), 2), 293.30)
  # gamma_profiles uses height above the site = Phi / g - altitude
  t0 <- as.POSIXct("2022-06-01", tz = "UTC")
  p <- c(1000, 925, 850, 700, 600)
  z <- c(100, 800, 1500, 3000, 4200) # m a.s.l.
  pl <- dplyr::tibble(cell_row = 0L, cell_col = 0L, timestamp_utc = t0, p_hPa = p, zg_m = z,
                      theta_K = 290 + 0.004 * z)
  g <- gamma_profiles(pl, altitude_m = 0)
  expect_equal(g$gamma, 0.004)
  expect_equal(g$n_levels, 3L) # 800, 1500, 3000 m in [500, 3500]
  expect_equal(g$levels, "925;850;700")
  # site at 350 m: 925 hPa at 450 m above the site drops out; above an orography of 0 m it stays
  g2 <- gamma_profiles(pl, altitude_m = 350, orography_m = 0)
  expect_equal(g2$levels, "850;700")
  expect_equal(g2$levels_orog, "925;850;700")
  expect_equal(g2$gamma_orog, 0.004)
})

test_that("gamma_fit_profile: fit layer, levels below ground excluded, 2 and 3 levels, plausibility", {
  o <- gamma_options()
  z <- c(-120, 300, 800, 1500, 3000, 4200)
  p <- c(850, 800, 750, 700, 600, 500)
  th <- 300 + 0.005 * z
  f <- gamma_fit_profile(z, th, p, o)
  expect_equal(f$gamma, 0.005)
  expect_equal(f$n_levels, 3L)
  expect_equal(f$r2, 1)
  expect_equal(f$flag, "ok")
  # levels below the ground are never used, even with a fit layer from 0 m
  f0 <- gamma_fit_profile(z, th, p, gamma_options(zmin_agl = 0, zmax_agl = 1000))
  expect_equal(f0$n_levels, 2L) # 300, 800 m (not -120 m)
  expect_equal(f0$levels, "800;750")
  # 2 levels: slope of the two points, R2 = 1
  f2 <- gamma_fit_profile(c(800, 1500), c(300, 302.1), c(925, 850), o)
  expect_equal(f2$gamma, 2.1 / 700)
  expect_equal(f2$n_levels, 2L)
  # 3 levels, not collinear: least squares
  f3 <- gamma_fit_profile(c(800, 1500, 3000), c(300, 303, 309), c(925, 850, 700), o)
  expect_equal(f3$gamma, unname(stats::coef(stats::lm(c(300, 303, 309) ~ c(800, 1500, 3000)))[2]))
  expect_true(f3$r2 < 1)
  # too few levels; implausible gamma -> NA with the raw value kept
  expect_equal(gamma_fit_profile(c(800, 4000), c(300, 310), c(925, 600), o)$flag, "too_few_levels")
  fi <- gamma_fit_profile(c(800, 1500), c(300, 300.5), c(925, 850), o) # 0.0007 < 0.001
  expect_true(is.na(fi$gamma))
  expect_equal(fi$gamma_raw, 0.5 / 700)
  expect_equal(fi$flag, "implausible")
  expect_true(is.na(gamma_fit_profile(c(800, 1500), c(300, 330), c(925, 850), o)$gamma)) # 0.043 > 0.02
})

test_that("gamma_daily: the profile closest to and not after sunrise, by site and season", {
  t0 <- as.POSIXct("2022-06-18", tz = "UTC")
  times <- c(t0 + (0:15) * 6 * 3600, as.POSIXct("2022-12-18", tz = "UTC") + (0:15) * 6 * 3600)
  prof <- dplyr::tibble(timestamp_utc = times, gamma = seq_along(times) / 1e4, gamma_raw = NA_real_,
                        n_levels = 3L, r2 = 1, levels = "925;850;700", flag = "ok", gamma_orog = NA_real_,
                        levels_orog = NA_character_)
  hour_of <- function(x) as.integer(format(x, "%H", tz = "UTC"))
  sites <- list(FI = c(61.85, 24.29), FR = c(48.84, 1.95), CH = c(46.82, 9.86))
  for (s in names(sites)) {
    g <- gamma_daily(prof, as.Date(c("2022-06-21", "2022-12-21")), sites[[s]][1], sites[[s]][2])
    expect_true(all(g$profile_utc <= g$sunrise_utc))
    expect_true(all(g$sunrise_utc - g$profile_utc < 6 * 3600))
    expect_equal(hour_of(g$profile_utc), c(0L, 6L)) # summer: 00 UTC; winter: 06 UTC
  }
  # forced UTC hour; a date without a profile
  g <- gamma_daily(prof, as.Date(c("2022-06-21", "2022-08-01")), 61.85, 24.29, gamma_options(profile_hour = 12))
  expect_equal(hour_of(g$profile_utc[1]), 12L)
  expect_equal(g$flag[2], "no_profile")
  expect_true(is.na(g$gamma[2]))
})

test_that("gamma_climatology and gamma_for_rows: monthly medians, daily with fallback", {
  daily <- dplyr::tibble(date = as.Date(c("2022-06-01", "2022-06-02", "2022-06-03", "2022-07-01")),
                         gamma = c(0.004, NA, 0.006, 0.003), flag = c("ok", "implausible", "ok", "ok"))
  cl <- gamma_climatology(daily)
  expect_equal(nrow(cl), 12)
  expect_equal(cl$median[6], 0.005)
  expect_equal(cl$n_days[6], 2L)
  expect_equal(cl$n_flagged[6], 1L)
  cl$median[is.na(cl$median)] <- 0.005
  day <- as.Date(c("2022-06-01", "2022-06-02", "2022-06-02", "2022-07-01", "2023-01-05"))
  g <- gamma_for_rows(day, "era5_daily", daily, cl)
  expect_equal(g$gamma, c(0.004, 0.005, 0.005, 0.003, 0.005))
  expect_equal(g$source_used, c("era5_daily", "era5_clim_fallback", "era5_clim_fallback", "era5_daily",
                                "era5_clim_fallback"))
  gc <- gamma_for_rows(day, "era5_clim", daily, cl)
  expect_equal(gc$gamma, c(0.005, 0.005, 0.005, 0.003, 0.005))
  expect_true(all(gc$source_used == "era5_clim"))
  # a month without a climatology value: constant gamma, flagged (G 2026-10-30)
  cl$median[1] <- NA
  expect_warning(gf <- gamma_for_rows(day, "era5_clim", daily, cl, gamma_const = 0.01), "missing")
  expect_equal(gf$gamma, c(0.005, 0.005, 0.005, 0.003, 0.01))
  expect_equal(gf$source_used[5], "constant_fallback")
  expect_warning(gd <- gamma_for_rows(day, "era5_daily", daily, cl, gamma_const = 0.02), "missing")
  expect_equal(gd$gamma[5], 0.02)
  expect_equal(gd$source_used, c("era5_daily", "era5_clim_fallback", "era5_clim_fallback", "era5_daily",
                                 "constant_fallback"))
})

test_that("run_pblh_model gamma_source: constant unchanged, ERA5 tables per day, bg_only keeps N constant", {
  site <- test_site_meta()
  d <- make_diurnal(3)
  ref <- run_pblh_model(d, site, pblh_options())
  expect_true(all(ref$gamma_used == 0.01) && all(ref$gamma_source_used == "constant"))
  # climatology equal to the constant: same output
  cl <- dplyr::tibble(month = 1:12, median = 0.01)
  daily <- dplyr::tibble(date = as.Date(character()), gamma = numeric())
  e <- list(daily = daily, clim = cl)
  same <- run_pblh_model(d, site, pblh_options(gamma_source = "era5_clim"), gamma_era5 = e)
  expect_equal(same$pblh_m, ref$pblh_m, tolerance = 1e-12)
  expect_true(all(same$gamma_source_used == "era5_clim"))
  # daily gamma on one day only: smaller gamma -> deeper B&G that day; other days fall back to the climatology
  days <- sort(unique(ref$day_id[ref$regime %in% "convective"]))
  e$daily <- dplyr::tibble(date = days[2], gamma = 0.0025)
  dl <- run_pblh_model(d, site, pblh_options(gamma_source = "era5_daily"), gamma_era5 = e)
  on2 <- dl$day_id == days[2]
  expect_true(all(dl$gamma_used[on2] == 0.0025) && all(dl$gamma_source_used[on2] == "era5_daily"))
  expect_true(all(dl$gamma_source_used[!on2] == "era5_clim_fallback"))
  conv2 <- on2 & dl$regime %in% "convective"
  expect_true(max(dl$pblh_m[conv2]) > max(ref$pblh_m[conv2]))
  expect_equal(dl$pblh_m[!on2], ref$pblh_m[!on2], tolerance = 1e-12)
  # bg_only: B&G gamma from ERA5, N (Eq. B3, B2/B3 switch) from the constant gamma
  b <- run_pblh_model(d, site, pblh_options(gamma_source = "era5_daily", gamma_scope = "bg_only"), gamma_era5 = e)
  expect_true(all(b$gamma_n_used == 0.01))
  expect_equal(b$h_bg_m, dl$h_bg_m, tolerance = 1e-12)
  expect_error(run_pblh_model(d, site, pblh_options(gamma_source = "era5_clim")), "gamma_era5")
  expect_error(pblh_options(gamma_source = "era5"), "gamma_source")
})
