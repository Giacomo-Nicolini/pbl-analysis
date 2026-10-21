site <- test_site_meta()

test_that("H_F_MDS: QC masking and fallback only where missing in the file", {
  d <- make_icos_like(4,
    fn_H_F_MDS = c(100, 100, NA, NA),
    fn_H_F_MDS_QC = c(0, 2, NA, NA),
    H = c(90, 95, 80, NA)
  )
  f <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 0))
  expect_equal(f$H, c(100, NA, 80, NA))
  expect_equal(f$H_source, c("H_F_MDS", NA, "H", NA))
  expect_equal(f$flag_H_qc_masked, c(FALSE, TRUE, FALSE, FALSE))

  f <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 0, h_fallback = FALSE))
  expect_equal(f$H, c(100, NA, NA, NA))
  f <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 0, h_qc_max = 2))
  expect_equal(f$H[2], 100)
})

test_that("h_source H and H_CORR", {
  d <- make_icos_like(2, H = c(10, 20), fn_H_CORR = c(30, 40))
  expect_equal(prepare_pblh_forcing(d, site, pblh_options(h_source = "H"))$H, c(10, 20))
  f <- prepare_pblh_forcing(d, site, pblh_options(h_source = "H_CORR"))
  expect_equal(f$H, c(30, 40))
  expect_equal(f$H_source, c("H_CORR", "H_CORR"))
})

test_that("T and p fallbacks, barometric estimate", {
  d <- make_icos_like(3,
    fn_TA_F_K = c(290, NA, NA), met_TA_K = c(289, 288, NA), T_SONIC_K = 291,
    fn_PA_F_PA = c(100000, NA, NA), met_PA_PA = c(99000, 98000, NA)
  )
  f <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 0))
  expect_equal(f$T_K, c(290, 288, 291))
  expect_equal(f$T_source, c("TA_F", "TA_meteo", "T_SONIC"))
  # barometric: 101325 * exp(-16 / 8434) = 101132.96 Pa
  expect_equal(f$p_Pa, c(100000, 98000, 101325 * exp(-16 / 8434)))
  expect_equal(f$p_Pa[3], 101132.96, tolerance = 1e-6)
  expect_equal(f$p_source, c("PA_F", "PA_meteo", "barometric"))

  site_na <- site
  site_na$altitude_m <- NA
  expect_error(prepare_pblh_forcing(d, site_na, pblh_options(max_gap_interp = 0)), "altitude_m")
})

test_that("density, kinematic heat flux and Coriolis parameter (hand-computed)", {
  d <- make_icos_like(1, fn_TA_F_K = 300, fn_PA_F_PA = 100000, fn_H_F_MDS = 100)
  f <- prepare_pblh_forcing(d, site)
  # rho = 100000 / (287.05 * 300) = 100000 / 86115 = 1.161238 kg m-3
  expect_equal(f$rho, 1.161238, tolerance = 1e-6)
  # wtheta = 100 / (1.161238 * 1005) = 100 / 1167.044 = 0.0856866 K m s-1
  expect_equal(f$wtheta, 0.0856866, tolerance = 1e-6)
  # f = 2 * 7.2921159e-5 * sin(51.31 deg) = 1.138365e-4 s-1
  expect_equal(f$f_cor, 1.138365e-4, tolerance = 1e-6)
})

test_that("short gaps interpolated and flagged; long gaps and edges not", {
  d <- make_icos_like(10, USTAR = c(0.1, NA, NA, 0.4, NA, NA, NA, 0.8, 0.9, NA))
  f <- prepare_pblh_forcing(d, site)
  expect_equal(f$ustar[1:4], c(0.1, 0.2, 0.3, 0.4))
  expect_equal(f$flag_ustar_interp, c(FALSE, TRUE, TRUE, rep(FALSE, 7)))
  expect_true(all(is.na(f$ustar[5:7]))) # gap of 3 > 2
  expect_true(is.na(f$ustar[10]))       # trailing gap: no extrapolation
})

test_that("L is interpolated through 1/L", {
  # 1/L: 0.01, (0.015), 0.02 -> L = 66.67 m (linear in L would give 75 m)
  # l_source = "L2": synthetic L not consistent with the constant H of the fixture (spec 06 check off)
  o <- pblh_options(l_source = "L2")
  d <- make_icos_like(3, MO_LENGTH = c(100, NA, 50))
  f <- prepare_pblh_forcing(d, site, o)
  expect_equal(f$L[2], 1 / 0.015)
  expect_equal(f$flag_L_interp, c(FALSE, TRUE, FALSE))
  expect_equal(f$L_source, c("L2", "interpolated", "L2"))
  # across neutral: -1000, NA, 1000 -> 1/L = 0 -> L = Inf (neutral), not ~0
  d <- make_icos_like(3, MO_LENGTH = c(-1000, NA, 1000))
  expect_equal(prepare_pblh_forcing(d, site, o)$L[2], Inf)
})

test_that("excluded rows: all forcing NA, never interpolated across", {
  d <- make_icos_like(5, excl_flag = c(FALSE, TRUE, FALSE, FALSE, FALSE),
                      USTAR = c(0.1, 0.2, NA, 0.4, 0.5))
  f <- prepare_pblh_forcing(d, site)
  expect_true(all(is.na(unlist(f[2, c("ustar", "L", "H", "T_K", "p_Pa", "rho", "wtheta")]))))
  expect_true(is.na(f$ustar[3])) # gap adjacent to an excluded row is not filled
  expect_false(f$flag_ustar_interp[3])
  expect_equal(f$excl_flag, c(FALSE, TRUE, FALSE, FALSE, FALSE))
})

test_that("pblh_options validates", {
  expect_error(pblh_options(stable_method = "foo"), "stable_method")
  expect_error(pblh_options(gamma = -1), "gamma")
  expect_error(pblh_options(n_sub = 0), "n_sub")
  expect_error(pblh_options(h_source = "LE"), "h_source")
  expect_error(pblh_options(combine = "min"), "combine")
  expect_error(pblh_options(neutral_method = "B9"), "neutral_method")
  expect_equal(pblh_options()$combine, "bg_only")
  expect_equal(pblh_options()$neutral_method, "auto")
  expect_error(pblh_options(foo = 1), "Unknown")
  expect_equal(pblh_options()$gamma, 0.01)
  expect_equal(pblh_options(gamma = c(0.01, 0.02))$gamma, c(0.01, 0.02))
})

test_that("u* ratio test and cap: hand cases", {
  o <- pblh_options()
  # u* 2.5 / WS 4 = 0.625 > 0.5 -> rejected; 0.8 / 4 = 0.2 -> kept; 0.8 / 0.5 = 1.6 but WS < 1 -> kept;
  # WS NA -> only the cap (u* 2.5 < 3: kept); 3.5 with WS 20 (ratio 0.175) -> rejected by the cap
  r <- ustar_implausible(c(2.5, 0.8, 0.8, 2.5, 3.5, NA), c(4, 4, 0.5, NA, 20, 4), o)
  expect_equal(r$ratio, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(r$cap, c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  # WS exactly at the threshold is tested; ratio exactly 0.5 is kept (strict >)
  r2 <- ustar_implausible(c(0.6, 0.5), c(1, 1), o)
  expect_equal(r2$ratio, c(TRUE, FALSE))
  none <- ustar_implausible(c(2.5, 49), c(4, 4), pblh_options(ustar_ws_max = Inf, ustar_max = Inf))
  expect_false(any(none$ratio | none$cap))
})

test_that("rejected u*: u* and L set to NA, H untouched, flags; short gaps interpolated, long runs NA", {
  d <- make_icos_like(10, USTAR = c(0.3, 2.5, 0.5, 0.3, 4, 4, 4, 0.3, 0.3, 0.3),
                      WS = c(3, 4, 4, 3, 20, 20, 20, 3, 3, 3),
                      MO_LENGTH = c(-100, -50, -80, -100, -60, -60, -60, -100, -100, -100),
                      fn_H_F_MDS = 1:10 * 10)
  f <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 0, l_source = "L2"))
  expect_equal(f$flag_ustar_implausible, c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(f$ustar, c(0.3, NA, 0.5, 0.3, NA, NA, NA, 0.3, 0.3, 0.3))
  expect_true(all(is.na(f$L[c(2, 5:7)])))
  expect_equal(f$H, 1:10 * 10)
  expect_equal(f$WS, c(3, 4, 4, 3, 20, 20, 20, 3, 3, 3))

  f2 <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 2, l_source = "L2"))
  # row 2: single rejected row, interpolated (0.3 -> 0.5), keeps both flags; L via 1/L
  expect_equal(f2$ustar[2], 0.4)
  expect_true(f2$flag_ustar_implausible[2] && f2$flag_ustar_interp[2] && f2$flag_L_interp[2])
  expect_equal(f2$L[2], 1 / mean(c(1 / -100, 1 / -80)))
  # rows 5-7: run of 3 > max_gap_interp, stay NA
  expect_true(all(is.na(f2$ustar[5:7])) && all(is.na(f2$L[5:7])))
  expect_false(any(f2$flag_ustar_interp[5:7]))
})

test_that("u* filter options are validated", {
  expect_error(pblh_options(ustar_ws_max = 0), "ustar_ws_max")
  expect_error(pblh_options(ustar_max = -1), "ustar_max")
  expect_error(pblh_options(ws_min_ratio_test = NA_real_), "ws_min_ratio_test")
  expect_error(pblh_options(ustar_ws_max = c(0.4, 0.5)), "ustar_ws_max")
  expect_no_error(pblh_options(ustar_ws_max = Inf, ustar_max = Inf))
})

test_that("L_forcing: hand case", {
  # rho0 = 100000 / (287.05 * 293.15) = 1.188372 kg m-3
  # wtheta0 = 200 / (1.188372 * 1005) = 0.167460 K m s-1
  # L = -0.5^3 * 293.15 / (0.4 * 9.80665 * 0.167460) = -36.64375 / 0.656890 = -55.784 m
  d <- make_icos_like(1, start = "2024-06-01 10:00:00", USTAR = 0.5, fn_H_F_MDS = 200, fn_TA_F_K = 293.15, fn_PA_F_PA = 1e5, MO_LENGTH = -50)
  f <- prepare_pblh_forcing(d, site)
  expect_equal(f$L_forcing, -55.7838, tolerance = 1e-5)
  expect_equal(f$L_forcing, obukhov_length(0.5, 200 / (1e5 / (287.05 * 293.15) * 1005), 293.15))
  expect_equal(f$L_check, "consistent") # -50 vs -55.8: ratio 1.12
  expect_equal(f$L, -50)
  expect_equal(f$L_source, "L2")
  expect_true(is.na(obukhov_length(0.3, 0, 290)))
})

test_that("l_consistency: sign, magnitude (both directions), near-neutral not checked", {
  o <- pblh_options()
  lf <- rep(-50, 7)
  # L2: +50 (sign); -145 (2.9x, kept); -155 (3.1x); -50/2.9 (kept); -50/3.1; -50 near neutral; NA
  L <- c(50, -145, -155, -50 / 2.9, -50 / 3.1, -500, NA)
  w <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.005, 0.1)
  expect_equal(l_consistency(L, lf, w, o),
               c("sign", "consistent", "magnitude", "consistent", "magnitude", "not_checked", "not_checked"))
  expect_equal(l_consistency(-155, -50, 0.1, pblh_options(l_consistency_factor = 5)), "consistent")
})

test_that("inconsistent L2 L replaced by the forcing L and flagged; L_l2 kept", {
  # H = 200 W m-2 (wtheta0 ~ 0.167), u* 0.5: L_forcing ~ -55.8 m
  d <- make_icos_like(4, start = "2024-06-01 10:00:00", USTAR = 0.5, fn_H_F_MDS = 200, fn_TA_F_K = 293.15, fn_PA_F_PA = 1e5,
                      MO_LENGTH = c(-50, 40, -3e6, -500))
  f <- prepare_pblh_forcing(d, site)
  expect_equal(f$L_check, c("consistent", "sign", "magnitude", "magnitude"))
  expect_equal(f$flag_L_replaced, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(f$L_source, c("L2", "forcing", "forcing", "forcing"))
  expect_equal(f$L[2:4], f$L_forcing[2:4])
  expect_equal(f$L_l2, c(-50, 40, -3e6, -500))
  expect_false(any(f$flag_L_interp))
})

test_that("missing L2 L stays missing (never filled from the forcing), interpolated only in short gaps", {
  d <- make_icos_like(8, start = "2024-06-01 10:00:00", USTAR = 0.5, fn_H_F_MDS = 200, fn_TA_F_K = 293.15, fn_PA_F_PA = 1e5,
                      MO_LENGTH = c(-50, NA, -60, -50, NA, NA, NA, -50))
  f <- prepare_pblh_forcing(d, site)
  expect_equal(f$L_source, c("L2", "interpolated", "L2", "L2", NA, NA, NA, "L2"))
  expect_true(all(is.na(f$L[5:7])))
  expect_false(any(f$flag_L_replaced))
  expect_equal(f$L_check[c(2, 5)], c("not_checked", "not_checked"))
  expect_true(all(!is.na(f$L_forcing[5:7]))) # available, but not used
})

test_that("rows rejected by the u* filter: L NA, not checked", {
  d <- make_icos_like(3, start = "2024-06-01 10:00:00", USTAR = c(0.5, 2.5, 0.5), WS = 4, fn_H_F_MDS = 200, MO_LENGTH = 30)
  f <- prepare_pblh_forcing(d, site, pblh_options(max_gap_interp = 0))
  expect_true(f$flag_ustar_implausible[2])
  expect_true(is.na(f$L[2]) && is.na(f$L_l2[2]) && is.na(f$L_forcing[2]))
  expect_equal(f$L_check[2], "not_checked")
  expect_equal(f$L_check[c(1, 3)], c("sign", "sign"))
})

test_that("l_source L2 leaves L2 L as is; forcing uses L_forcing everywhere it exists", {
  d <- make_icos_like(4, start = "2024-06-01 10:00:00", USTAR = 0.5, fn_H_F_MDS = c(200, 200, 0, 200), fn_TA_F_K = 293.15, fn_PA_F_PA = 1e5,
                      MO_LENGTH = c(-50, 40, -3e6, NA))
  a <- prepare_pblh_forcing(d, site, pblh_options(l_source = "L2", max_gap_interp = 0))
  expect_equal(a$L, c(-50, 40, -3e6, NA))
  expect_false(any(a$flag_L_replaced))
  expect_equal(a$L_check, c("consistent", "sign", "not_checked", "not_checked")) # still reported
  b <- prepare_pblh_forcing(d, site, pblh_options(l_source = "forcing", max_gap_interp = 0))
  expect_equal(b$L[c(1, 2, 4)], b$L_forcing[c(1, 2, 4)]) # also where L2 L is missing
  expect_true(is.na(b$L[3]))                              # wtheta0 = 0: NA
  expect_equal(b$L_source, c("forcing", "forcing", NA, "forcing"))
  expect_equal(b$flag_L_replaced, c(TRUE, TRUE, FALSE, TRUE))
})

test_that("L consistency options are validated", {
  expect_error(pblh_options(l_source = "ERA5"), "l_source")
  expect_error(pblh_options(l_consistency_factor = 1), "l_consistency_factor")
  expect_error(pblh_options(l_check_wtheta_min = -0.1), "l_check_wtheta_min")
  expect_no_error(pblh_options(l_consistency_factor = 1.5, l_check_wtheta_min = 0))
})

test_that("stable rows are not checked: inconsistent L2 L kept (spec 06 amendment 1)", {
  # H = -20 W m-2, u* 0.5: L_forcing ~ +560 m; L2 L = +50 (11x) or -50 (sign): both kept
  d <- make_icos_like(3, start = "2024-06-01 10:00:00", USTAR = 0.5, fn_H_F_MDS = c(-20, -20, -5), fn_TA_F_K = 293.15, fn_PA_F_PA = 1e5,
                      MO_LENGTH = c(50, -50, 40))
  f <- prepare_pblh_forcing(d, site)
  expect_true(all(f$L_forcing[1:2] > 0))
  expect_equal(f$L_check, c("not_checked", "not_checked", "not_checked"))
  expect_false(any(f$flag_L_replaced))
  expect_equal(f$L, c(50, -50, 40))
  expect_equal(l_consistency(c(50, -50), c(560, 560), c(-0.02, -0.02)), c("not_checked", "not_checked"))
})

test_that("night rows are not checked (spec 06 amendment 2); same row in daytime replaced; FALSE = V1", {
  # positive H (+35 W m-2 gap-filled at night) with a stable L2 L: sign inconsistent
  mk <- function(start) make_icos_like(1, start = start, USTAR = 0.4, fn_H_F_MDS = 35, fn_TA_F_K = 275,
                                       fn_PA_F_PA = 1e5, MO_LENGTH = 80)
  night <- prepare_pblh_forcing(mk("2023-01-20 01:00:00"), site)
  expect_false(night$is_day)
  expect_equal(night$L_check, "not_checked")
  expect_false(night$flag_L_replaced)
  expect_equal(night$L, 80)
  day <- prepare_pblh_forcing(mk("2023-01-20 11:00:00"), site)
  expect_true(day$is_day)
  expect_equal(day$L_check, "sign")
  expect_true(day$flag_L_replaced)
  expect_equal(day$L, day$L_forcing)
  v1 <- prepare_pblh_forcing(mk("2023-01-20 01:00:00"), site, pblh_options(l_check_daytime_only = FALSE))
  expect_equal(v1$L_check, "sign")
  expect_true(v1$flag_L_replaced)
  expect_error(pblh_options(l_check_daytime_only = NA), "l_check_daytime_only")
  expect_error(pblh_options(l_check_daytime_only = "yes"), "l_check_daytime_only")
  # l_consistency without is_day: no daytime restriction
  expect_equal(l_consistency(80, -50, 0.03), "sign")
})
