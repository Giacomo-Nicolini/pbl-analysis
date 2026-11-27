# Spec 17: EC-derived gamma (theta conversion, Methods A and B, window rules, flags, climatology, model options).

test_that("dry potential temperature on a hand case", {
  # 300 K at 900 hPa: 300 (1e5 / 9e4)^(287.05 / 1005) = 309.1652 K
  expect_equal(theta_dry(300, 90000), 300 * (1e5 / 9e4)^(287.05 / 1005), tolerance = 1e-12)
  expect_equal(theta_dry(300, 90000), 309.1652, tolerance = 1e-6)
  expect_equal(theta_dry(285, 1e5), 285)
})

# h0 as run_pblh_model(): diagnostic height of the half-hour before sunrise, else the nearest valid within 2 h
model_h0 <- function(pm, s) {
  cand <- seq.int(s - 1, max(1, s - 4))
  cand <- cand[cand >= 1 & !is.na(pm$h_diag_m[cand])]
  pm$h_diag_m[cand[1]]
}

# Synthetic convective forcing of one window (rows of 30 min)
syn_forcing <- function(n = 12, w_max = 0.15, ustar = 0.35, T_K = 295) {
  ph <- seq(0.15, 0.95, length.out = n)
  w <- w_max * sin(pi * ph)                      # K m s-1
  L <- -ustar^3 * T_K / (0.4 * 9.80665 * w)      # Obukhov length of the kinematic flux [m]
  list(w = w, ustar = rep(ustar, n), L = L, T_K = rep(T_K, n))
}

test_that("Method B recovers gamma on a synthetic encroachment day within 1 %", {
  A <- pblh_constants()$bg_A
  f <- syn_forcing()
  for (g in c(0.003, 0.006, 0.012)) {
    q_mid <- c(0, cumsum(f$w * 1800))[seq_along(f$w)] + f$w * 900
    h <- sqrt(2 * (1 + 2 * A) * q_mid / g)       # h0 = 0
    theta <- 290 + g * (1 + A) / (1 + 2 * A) * h
    fb <- gamma_ec_fit_b(theta, f$w)
    expect_equal(fb$gamma, g, tolerance = 0.01)
    expect_gt(fb$r2, 0.999)
  }
  # the conversion factor (1 + 2A) / (1 + A)^2 is 0.97 for A = 0.2
  expect_equal((1 + 2 * A) / (1 + A)^2, 0.9722, tolerance = 1e-4)
  expect_true(is.na(gamma_ec_fit_b(rep(290, 12), f$w)$gamma) || gamma_ec_fit_b(rep(290, 12), f$w)$slope == 0)
})

test_that("Method A recovers gamma on a synthetic bg_step day (mechanical terms on, h0 > 0)", {
  f <- syn_forcing()
  for (g in c(0.003, 0.0075, 0.015)) {
    ph <- bg_theta_path(150, f$w, f$ustar, f$L, f$T_K, g)
    fa <- gamma_ec_fit_a(291 + ph$c_mean, 150, f$w, f$ustar, f$L, f$T_K)
    expect_equal(fa$gamma, g, tolerance = 0.02)
    expect_equal(fa$theta0, 291, tolerance = 1e-3)
    expect_gt(fa$r2, 0.999)
  }
})

test_that("Method A with 0.1 K noise: median bias over 200 synthetic days < 5 %", {
  set.seed(17)
  g_true <- 0.006
  rel <- vapply(seq_len(200), function(k) {
    f <- syn_forcing(w_max = runif(1, 0.08, 0.2), ustar = runif(1, 0.2, 0.5))
    h0 <- runif(1, 80, 300)
    ph <- bg_theta_path(h0, f$w, f$ustar, f$L, f$T_K, g_true)
    th <- 290 + ph$c_mean + stats::rnorm(length(f$w), 0, 0.1)
    fa <- gamma_ec_fit_a(th, h0, f$w, f$ustar, f$L, f$T_K)
    fa$gamma_opt / g_true - 1
  }, numeric(1))
  expect_lt(abs(stats::median(rel)), 0.05)
})

test_that("optimum at a bound gives NA and the bound side", {
  f <- syn_forcing()
  ph <- bg_theta_path(150, f$w, f$ustar, f$L, f$T_K, 0.05) # steeper than gamma_max = 0.02
  fa <- gamma_ec_fit_a(290 + ph$c_mean, 150, f$w, f$ustar, f$L, f$T_K)
  expect_true(fa$at_bound)
  expect_true(is.na(fa$gamma))
  expect_equal(fa$bound, "high")
  ph2 <- bg_theta_path(150, f$w, f$ustar, f$L, f$T_K, 0.0003) # below gamma_min = 0.001
  fa2 <- gamma_ec_fit_a(290 + ph2$c_mean, 150, f$w, f$ustar, f$L, f$T_K)
  expect_equal(fa2$bound, "low")
  expect_true(is.na(fa2$gamma))
})

test_that("the window integrator reproduces the B&G heights of run_pblh_model()", {
  site <- test_site_meta()
  out <- run_pblh_model(make_diurnal(3), site, pblh_options())
  d2 <- sort(unique(out$day_id[out$regime %in% "convective"]))[2] # a day with a night before it
  rd <- which(out$day_id == d2)
  conv <- rd[out$regime[rd] %in% "convective"]
  p <- conv[out$wtheta[conv] > 0][1]
  w <- p:max(conv)
  h0 <- model_h0(out, rd[out$is_day[rd] %in% TRUE][1])
  expect_false(is.na(h0))
  ph <- bg_theta_path(h0, out$wtheta[w], out$ustar[w], out$L[w], out$T_K[w], 0.01)
  expect_equal(ph$h_mean, out$h_bg_m[w], tolerance = 1e-9)
})

# Synthetic model output whose temperature follows the slab model with a known gamma (p = 1e5 Pa: theta = T)
syn_pm <- function(g_true = 0.006, n_days = 3, site = test_site_meta()) {
  pm <- run_pblh_model(make_diurnal(n_days, h_max = 250), site, pblh_options())
  for (d in unique(pm$day_id[pm$regime %in% "convective"])) {
    rd <- which(pm$day_id == d)
    conv <- rd[pm$regime[rd] %in% "convective"]
    p <- conv[pm$wtheta[conv] > 0][1]
    w <- p:max(conv)
    h0 <- model_h0(pm, rd[pm$is_day[rd] %in% TRUE][1])
    th <- rep(290, length(w))
    for (it in 1:4) th <- 290 + bg_theta_path(h0, pm$wtheta[w], pm$ustar[w], pm$L[w], th, g_true)$c_mean # T in Eq. B5
    pm$T_K[w] <- th
    after <- rd[rd > max(w)]
    pm$T_K[after] <- th[length(th)] - 0.5 # cooler after the window
  }
  pm$p_Pa <- 1e5
  pm
}

test_that("gamma_ec_daily recovers gamma end to end; window rules and flags", {
  site <- test_site_meta()
  pm <- syn_pm(0.006, site = site)
  n <- nrow(pm)
  sw_pot <- rep(800, n)
  dl <- gamma_ec_daily(pm, sw_in = sw_pot * 0.8, sw_in_pot = sw_pot, site, gamma_ec_options(end_offset_h = 12))
  fit <- dl[!is.na(dl$gamma_A), ]
  expect_gt(nrow(fit), 0)
  expect_equal(fit$gamma_A, rep(0.006, nrow(fit)), tolerance = 0.02)
  expect_true(all(fit$used))
  expect_equal(dl$gamma, ifelse(dl$used, dl$gamma_A, NA_real_))
  expect_true(all(abs(fit$clear_index - 0.8) < 1e-12))
  # noon + offset: the window ends at the last half-hour whose midpoint is at or before solar noon + 3 h
  d3 <- gamma_ec_daily(pm, sw_pot * 0.8, sw_pot, site, gamma_ec_options(end_offset_h = 3))
  k <- which(!is.na(d3$window_end_utc))[1]
  noon <- sun_times(d3$date[k], site$lat, site$lon)$solar_noon_utc
  expect_lte(as.numeric(d3$window_end_utc[k]) + 900, as.numeric(noon) + 3 * 3600)
  expect_gt(as.numeric(d3$window_end_utc[k]) + 1800 + 900, as.numeric(noon) + 3 * 3600)
  # maximum theta earlier: the window ends there
  pm2 <- pm
  rd <- which(pm2$day_id == d3$date[k] & pm2$regime %in% "convective")
  pk <- rd[6]
  pm2$T_K[pk] <- pm2$T_K[pk] + 50
  d4 <- gamma_ec_daily(pm2, sw_pot * 0.8, sw_pot, site, gamma_ec_options(end_offset_h = 12))
  expect_equal(d4$window_end_utc[d4$date == d3$date[k]], as_utc_time(pm2$timestamp_start_utc[pk]))
  # cloudy, unknown clear index, gap, T_SONIC rows, flat theta
  dc <- gamma_ec_daily(pm, sw_pot * 0.3, sw_pot, site)
  expect_true(all(grepl("cloudy", dc$flags[!is.na(dc$window_start_utc)])) && !any(dc$used))
  du <- gamma_ec_daily(pm, rep(NA_real_, n), sw_pot, site)
  expect_true(all(grepl("clear_unknown", du$flags[!is.na(du$window_start_utc)])) && !any(du$used))
  pg <- pm
  pg$ustar[rd[3]] <- NA
  dg <- gamma_ec_daily(pg, sw_pot * 0.8, sw_pot, site)
  expect_true(grepl("gap", dg$flags[dg$date == d3$date[k]]) && is.na(dg$gamma_A[dg$date == d3$date[k]]))
  ps <- pm
  ps$T_source <- "T_SONIC"
  ds <- gamma_ec_daily(ps, sw_pot * 0.8, sw_pot, site)
  expect_true(all(grepl("few_rows", ds$flags[!is.na(ds$window_start_utc)])))
  expect_gt(sum(ds$n_tsonic), 0)
  pf <- pm
  pf$T_K <- 290
  df <- gamma_ec_daily(pf, sw_pot * 0.8, sw_pot, site)
  expect_true(all(grepl("rise_lt_min", df$flags[!is.na(df$window_start_utc)])) && !any(df$used))
  # min_n larger than the window: few_rows
  dn <- gamma_ec_daily(pm, sw_pot * 0.8, sw_pot, site, gamma_ec_options(min_n = 40))
  expect_true(all(grepl("few_rows", dn$flags[!is.na(dn$window_start_utc)])))
  # model flags carried over
  pn <- pm
  pn$flag_no_h0[pn$day_id == d3$date[k]] <- TRUE
  dh <- gamma_ec_daily(pn, sw_pot * 0.8, sw_pot, site, gamma_ec_options(end_offset_h = 12))
  expect_true(grepl("no_h0", dh$flags[dh$date == d3$date[k]]) && !dh$used[dh$date == d3$date[k]])
  fc <- gamma_ec_flag_counts(dc)
  expect_equal(fc$n_days[fc$flag == "cloudy"], sum(grepl("cloudy", dc$flags)))
})

test_that("climatology layout, minimum days and the ec_clim / ec_daily options of run_pblh_model()", {
  daily <- dplyr::tibble(date = as.Date("2024-06-01") + 0:39, gamma_A = rep(c(0.004, 0.006), 20),
                         used = c(rep(TRUE, 30), rep(FALSE, 10)))
  daily$gamma <- ifelse(daily$used, daily$gamma_A, NA_real_)
  cl <- gamma_ec_climatology(daily, min_days = 5)
  expect_named(cl, c("month", "median", "p25", "p75", "n_days", "n_fitted", "n_flagged"))
  expect_equal(nrow(cl), 12)
  expect_equal(cl$n_days[cl$month == 6], 30)
  expect_equal(cl$n_days[cl$month == 7], 0)
  expect_equal(cl$n_flagged[cl$month == 7], 10)
  expect_equal(cl$median[cl$month == 6], 0.005)
  expect_true(is.na(cl$median[cl$month == 7]))
  expect_true(is.na(gamma_ec_climatology(daily[1:4, ], min_days = 5)$median[6]))
  # model options: constant identity, climatology and daily values, fallbacks
  site <- test_site_meta()
  d <- make_diurnal(3)
  ref <- run_pblh_model(d, site, pblh_options())
  e <- list(daily = dplyr::tibble(date = as.Date(character()), gamma = numeric()),
            clim = dplyr::tibble(month = 1:12, median = 0.01))
  same <- run_pblh_model(d, site, pblh_options(gamma_source = "ec_clim"), gamma_era5 = e)
  expect_equal(same$pblh_m, ref$pblh_m, tolerance = 1e-12)
  expect_true(all(same$gamma_source_used == "ec_clim"))
  days <- sort(unique(ref$day_id[ref$regime %in% "convective"]))
  e$daily <- dplyr::tibble(date = days[2], gamma = 0.003)
  dl <- run_pblh_model(d, site, pblh_options(gamma_source = "ec_daily"), gamma_era5 = e)
  on2 <- dl$day_id == days[2]
  expect_true(all(dl$gamma_source_used[on2] == "ec_daily") && all(dl$gamma_source_used[!on2] == "ec_clim_fallback"))
  e$clim$median[6] <- NA
  expect_warning(fb <- run_pblh_model(d, site, pblh_options(gamma_source = "ec_clim"), gamma_era5 = e),
                 "EC-derived gamma climatology missing")
  m6 <- format(fb$day_id, "%m") == "06" # the night before the first sunrise belongs to 31 May
  expect_true(all(fb$gamma_source_used[m6] == "constant_fallback") && all(fb$gamma_source_used[!m6] == "ec_clim"))
  expect_equal(fb$pblh_m, ref$pblh_m, tolerance = 1e-12)
  # a 17a monthly table is read as the clim table
  tmp <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(dplyr::mutate(gamma_ec_climatology(daily), site_id = "XX-Tst", .before = 1), tmp)
  rcl <- readr::read_csv(tmp, show_col_types = FALSE)
  g <- gamma_for_rows(as.Date("2024-06-15"), "ec_clim", e$daily, rcl)
  expect_equal(g$gamma, 0.005)
  expect_error(pblh_options(gamma_source = "ec"), "gamma_source")
})

test_that("gamma benchmark: zero crossing and bootstrap interval", {
  expect_equal(zero_crossing(c(0.002, 0.004, 0.006), c(0.4, 0.1, -0.2)), 0.004 + 0.002 / 3)
  expect_true(is.na(zero_crossing(1:3, c(1, 2, 3))))
  expect_equal(zero_crossing(1:3, c(1, 0, -1)), 2)
  grid <- stats::setNames(c(0.002, 0.004, 0.006, 0.008), sprintf("gamma_g%03d", c(2, 4, 6, 8)))
  wide <- dplyr::tibble(site_id = "XX-Tst", case = names(grid), reference = "meas", stratum_set = "window_convective",
                        median_log_ratio = c(0.3, 0.1, -0.1, -0.3), mlr_lo = c(0.2, 0, -0.2, -0.4),
                        mlr_hi = c(0.4, 0.2, 0, -0.2), n_days = 50)
  wide <- dplyr::bind_rows(wide, dplyr::mutate(wide, stratum_set = "convective", median_log_ratio = 9))
  b <- gamma_benchmark_table(wide, c(`XX-Tst` = "window_convective"), grid)
  expect_equal(b$best$gamma_best, 0.005)
  expect_equal(b$best$gamma_best_lo, 0.004)
  expect_equal(b$best$gamma_best_hi, 0.006)
  expect_equal(nrow(b$curve), 4)
  tab <- dplyr::mutate(b$best, gamma_ec_median = 0.0055, gamma_ec_p25 = 0.004, gamma_ec_p75 = 0.007, gamma_site = 0.005)
  expect_s3_class(ggplot2::ggplot_build(plot_gamma_benchmark(tab, 0.0046)), "ggplot_built")
})

test_that("energy-balance ratio and write-if-changed", {
  t <- as.POSIXct("2024-07-01", tz = "UTC") + 1800 * 0:3
  d <- dplyr::tibble(timestamp_start_utc = t, fn_SW_IN_POT = c(0, 500, 500, 500), fn_H_F_MDS = 100,
                     fn_H_F_MDS_QC = c(0, 0, 0, 1), fn_LE_F_MDS = 200, fn_LE_F_MDS_QC = 0, fn_NETRAD = 450,
                     fn_G_F_MDS = c(50, 50, NA, 50))
  e <- energy_balance_ratio(d)
  expect_equal(e$n_rows, 1)
  expect_equal(e$ebr, 300 / 400)
  expect_equal(e$n_rows_no_g, 2)
  expect_equal(e$ebr_no_g, 600 / 900)
  expect_true(e$g_available)
  p <- withr::local_tempfile(fileext = ".csv")
  expect_true(write_csv_if_changed(dplyr::tibble(a = 1), p))
  m1 <- file.mtime(p)
  Sys.sleep(1.1)
  expect_false(write_csv_if_changed(dplyr::tibble(a = 1), p))
  expect_equal(file.mtime(p), m1)
  expect_true(write_csv_if_changed(dplyr::tibble(a = 2), p))
})

test_that("bracketing check of part E on hand rows", {
  r <- dplyr::tibble(site_id = "XX-Tst", zl_class = c("stable", "stable", "neutral", "unstable", "stable"),
                     daynight = c("night", "night", "night", "day", "night"),
                     h__era5 = c(100, 100, 300, 1000, 100), h__pblhmod = c(200, 200, 600, 800, 200),
                     h__pblhmod_clim = c(150, 150, 500, 900, 150))
  for (cs in c("pblhmod", "pblhmod_clim", "era5")) {
    r[[paste("v", "etc", "valid_qc_used", cs, sep = "__")]] <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
    r[[paste0("ok__", cs)]] <- TRUE
  }
  h_ec <- c(50, 150, 1200, 900, 10) # below by ln 2, inside, above by ln 2, inside, not in the row set
  b <- gamma_ec_bracketing(r, h_ec)
  a <- b[b$regime == "all" & b$daynight == "all", ]
  expect_equal(a$n, 4)
  expect_equal(a$share_outside, 0.5)
  expect_equal(a$share_below, 0.25)
  expect_equal(a$share_above, 0.25)
  expect_equal(a$median_dist_ln, log(2))
  sn <- b[b$regime == "stable" & b$daynight == "night", ]
  expect_equal(sn$n, 2)
  expect_equal(sn$share_below, 0.5)
  expect_equal(b$share_outside[b$regime == "unstable" & b$daynight == "day"], 0)
  expect_true(is.na(b$median_dist_ln[b$regime == "unstable" & b$daynight == "day"]))
})

test_that("monthly gamma_EC vs ERA5 table and the site summary", {
  ec <- dplyr::tibble(month = 1:12, median = c(NA, rep(0.006, 10), NA), p25 = 0.005, p75 = 0.007, n_days = 6)
  era <- dplyr::tibble(month = 1:12, median = 0.004, n_days = 30)
  m <- gamma_ec_vs_era5_table(ec, era, "XX-Tst")
  expect_equal(nrow(m), 12)
  expect_equal(m$ratio_ec_era5[2], 1.5)
  expect_true(is.na(m$ratio_ec_era5[1]))
  daily <- dplyr::tibble(gamma_A = c(0.004, 0.006, 0.008, NA), gamma_B = c(0.004, 0.003, 0.008, NA),
                         used = c(TRUE, TRUE, TRUE, FALSE))
  s <- gamma_ec_site_summary(m, daily, dplyr::tibble(gamma = c(0.004, 0.005, 0.006)),
                             dplyr::tibble(ebr = 0.8, ebr_no_g = 0.75, n_rows = 100, g_available = TRUE))
  expect_equal(s$n_used_days, 3)
  expect_equal(s$months_with_value, 10)
  expect_equal(s$gamma_ec_site, 0.006)
  expect_equal(s$gamma_era5_site, 0.005)
  expect_equal(s$ratio_ec_era5, 1.5)
  expect_equal(s$ratio_B_A, 1)
  expect_equal(s$ebr, 0.8)
  expect_s3_class(ggplot2::ggplot_build(plot_gamma_ec_vs_era5(m)), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(plot_gamma_ec_seasonal(m)), "ggplot_built")
})

test_that("offset folders and heights of the spec-17 supplement", {
  expect_equal(gamma_ec_dir_name(0), "gamma_ec")
  expect_equal(gamma_ec_dir_name(2), "gamma_ec_off2h")
  expect_equal(gamma_ec_dir_name(2.5), "gamma_ec_off2.5h")
  sites <- dplyr::tibble(site_id = "XX-Tst", utc_offset_h = 1, zm_m = 30, canopy_height_m = 15, igbp = "ENF")
  tv <- dplyr::tibble(site_id = "XX-Tst", variable = c("zm_m", "zm_m"), value = c(30, 36),
                      date_start = as.Date(c("2022-01-01", "2023-01-01")), date_end = as.Date(c("2023-01-01", NA)),
                      notes = "")
  daily <- dplyr::tibble(window_start_utc = as.POSIXct(c("2022-06-01 08:00", "2023-06-01 08:00", "2023-06-02 08:00",
                                                         "2023-06-03 08:00"), tz = "UTC"),
                         used = c(TRUE, TRUE, TRUE, FALSE))
  h <- gamma_ec_heights(daily, "XX-Tst", sites, tv)
  expect_equal(h$ec_height_m, 36)  # median of 30, 36, 36 (unused day ignored)
  expect_equal(h$canopy_height_m, 15)
  expect_equal(h$zm_m, 26)         # 36 - 2/3 * 15
  expect_equal(h$igbp, "ENF")
  expect_true(is.na(gamma_ec_heights(daily[0, ], "XX-Tst", sites, tv)$zm_m))
  x <- dplyr::tibble(site_id = rep(c("AA-One", "BB-Two"), 2), offset_h = rep(c(0, 2), each = 2),
                     ratio_ec_era5 = c(2, 3, 1.2, 1.5), zm_m = c(5, 20, 5, 20), ec_height_m = c(6, 35, 6, 35),
                     igbp = c("CRO", "ENF", "CRO", "ENF"))
  expect_s3_class(ggplot2::ggplot_build(plot_gamma_ratio_vs_height(x)), "ggplot_built")
})

test_that("land-cover groups, ratio-height statistics and caption (17s_gamma_ratio_vs_height)", {
  expect_equal(land_cover_group(c("CRO", "gra", "CSH", "DBF", "ENF", "MF", "EBF", "WET", "SAV")),
               c("lc_crop_grass", "lc_crop_grass", "lc_crop_grass", "lc_dbf", "lc_evergreen_mixed",
                 "lc_evergreen_mixed", "lc_evergreen_mixed", "lc_wetland", NA))
  x <- dplyr::tibble(site_id = rep(c("AA-One", "BB-Two", "CC-Thr", "DD-For", "EE-Fiv"), 2),
                     offset_h = rep(c(0, 2), each = 5),
                     ratio_ec_era5 = c(3.5, 3, 2, 1, NA, 3, 2.5, 1.5, 0.8, NA), zm_m = rep(c(2, 3, 20, 25, 5), 2),
                     igbp = rep(c("CRO", "GRA", "ENF", "ENF", "WET"), 2))
  st <- gamma_ratio_height_stats(x)
  expect_equal(st$rho$n_sites, c(4L, 4L))
  expect_equal(st$rho$rho_zm, c(-1, -1))
  g0 <- st$groups[st$groups$offset_h == 0, ]
  expect_equal(g0$group, c("lc_crop_grass", "lc_evergreen_mixed"))  # style order; the NA ratio row dropped
  expect_equal(g0$median_ratio, c(3.25, 1.5))
  cap <- gamma_ratio_height_caption(x)
  expect_match(cap, "-1.00", fixed = TRUE)
  expect_match(cap, "crop / grass 3.25 (2)", fixed = TRUE)
  p <- plot_gamma_ratio_vs_height(x, cap)
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[2]]), 4L)                                # one segment per site with both fits
  expect_setequal(b$plot$layers[[4]]$data$site_id, c("AA-One", "DD-For"))  # labels: < 1.5 or > 3.2 at +0 h
  expect_warning(plot_gamma_ratio_vs_height(dplyr::mutate(x, igbp = "SAV")), "outside the land-cover groups")
})
