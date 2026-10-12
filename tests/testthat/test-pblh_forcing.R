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
  d <- make_icos_like(3, MO_LENGTH = c(100, NA, 50))
  f <- prepare_pblh_forcing(d, site)
  expect_equal(f$L[2], 1 / 0.015)
  expect_equal(f$flag_L_interp, c(FALSE, TRUE, FALSE))
  # across neutral: -1000, NA, 1000 -> 1/L = 0 -> L = Inf (neutral), not ~0
  d <- make_icos_like(3, MO_LENGTH = c(-1000, NA, 1000))
  expect_equal(prepare_pblh_forcing(d, site)$L[2], Inf)
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
