utc <- function(x) as.POSIXct(x, tz = "UTC")

test_that("footprint_heights uses d = 2/3 hc, z0 = 0.15 hc, zm = EC height - d", {
  h <- footprint_heights(30, 15)
  expect_equal(c(h$d_m, h$z0_m, h$zm_m), c(10, 2.25, 20))
})

test_that("ERA5 BLH interpolated to half-hour midpoints, NA without both neighbours", {
  et <- utc(c("2022-09-01 00:00:00", "2022-09-01 01:00:00"))
  h <- era5_at_halfhour_mid(utc(c("2022-09-01 00:00:00", "2022-09-01 00:30:00", "2022-09-01 01:00:00")), et, c(100, 200))
  expect_equal(h, c(125, 175, NA))
})

test_that("footprint_h_cases builds the cases (scaled PBLHmod, measured exact / bound)", {
  t <- utc("2022-09-01 00:00:00") + 1800 * (0:2)
  std <- function(v, q = NA_integer_) dplyr::tibble(timestamp_start_utc = t, pblh_m = v, qc_flag = q)
  era <- dplyr::tibble(timestamp_start_utc = utc("2022-09-01 00:00:00") + 3600 * (0:2), pblh_m = c(100, 200, 300))
  meas <- std(c(500, NA, 800), c(0L, 1L, 3L))
  x <- footprint_h_cases(t, std(c(400, 500, NA)), era, era, measured = meas)
  expect_equal(x$pblhmod, c(400, 500, NA))
  expect_equal(x$`pblhmod_x0.8`, c(320, 400, NA))
  expect_equal(x$`pblhmod_x1.2`, c(480, 600, NA))
  expect_equal(x$era5, c(125, 175, 225))
  expect_equal(x$measured, c(500, NA, NA))
  expect_equal(x$bound, c(NA, NA, 800))
  expect_false("oper" %in% names(x))
})

test_that("etc_z0_multiplier reproduces the legacy loop", {
  expect_equal(etc_z0_multiplier(rep(10, 4), rep(1, 4)), 5) # all retained
  # ratios 10, 10, 10, 3.33, 3.33: rejected 2 > 0.5 * 3 until m <= 3.33
  expect_equal(etc_z0_multiplier(rep(10, 5), c(1, 1, 1, 3, 3)), 3.3)
  # all rows below any multiplier: 30 decrements to the minimum 2
  expect_equal(etc_z0_multiplier(rep(1, 3), rep(1, 3)), 2)
  # ratio exactly 2 with half retained: rejected 1 <= 0.5 * 2 at m = 5
  expect_equal(etc_z0_multiplier(c(10, 10, 2), c(1, 1, 1)), 5)
})

test_that("footprint_validity: conditions at the boundaries", {
  n <- 8
  v <- footprint_validity(
    zm = c(20, 20, 20, 20, 20, 20, 20, 20), z0 = c(0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), ec_height = 30,
    h = c(1000, 1000, 1000, 25, 9, 1000, 1000, 1000), umean = 3, L = c(100, 100, 20, 100, 100, 100, 19.8, -20 / 15.5),
    sigmav = 0.5, ustar = c(0.3, 0.3, 0.1, 0.3, 0.3, 0.2, 0.3, 0.3), wind_dir = 180, m_etc = 5
  )
  expect_equal(v$k_rsl, c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))   # 20 z0 < zm (zm = 20 z0: FALSE)
  expect_equal(v$e_el[4], TRUE)                                               # zm = 20 = 0.8 * 25 (ETC: invalid only if zm > 0.8 h)
  expect_equal(v$e_tower[4], FALSE)                                           # EC 30 > h 25
  expect_equal(v$code_h10[5], FALSE)                                          # h < 10
  expect_equal(v$domain[5], FALSE)                                            # zm >= h
  expect_equal(c(v$code_ustar01[3], v$e_ustar[3]), c(TRUE, FALSE))             # u* = 0.1
  expect_equal(v$e_ustar[6], TRUE)                                            # u* = 0.2
  expect_equal(c(v$e_stab[3], v$e_stab[7]), c(TRUE, FALSE))                   # zm/L = 1 vs 1.01
  expect_equal(c(v$k_stab[8], v$code_stab[8]), c(TRUE, FALSE))                # zm/L = -15.5
  expect_equal(v$valid_kljun, c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(v$valid_kljun_norsl, c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(v$valid_etc, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(v$valid_none, c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE))
  # missing input, implausible u*, exclusion
  w <- footprint_validity(20, 0.5, 30, c(NA, 1000, 1000), 3, 100, 0.5, 0.3, 180, 5,
                          ustar_implausible = c(FALSE, TRUE, FALSE), excluded = c(FALSE, FALSE, TRUE))
  expect_equal(w$valid_none, c(FALSE, FALSE, FALSE))
})

test_that("footprint_subsets are nested", {
  s <- footprint_subsets(c(TRUE, TRUE, TRUE, TRUE, FALSE), nee = c(1, 1, NA, 1, 1), nee_data_flag = c(0, 0, 0, 2, 0),
                         nee_vut_ref_qc = c(0, 1, 0, 0, 0))
  expect_equal(s$valid, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(s$valid_qc, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(s$valid_qc_used, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_true(all(s$valid >= s$valid_qc & s$valid_qc >= s$valid_qc_used))
})

test_that("ffp_ci_paper: distances scale with zm / (1 - zm/h) umean / u*", {
  a <- ffp_ci_paper(zm = 20, h = c(1000, 200, 20, NA), umean = 3, ustar = 0.3)
  expect_true(all(is.na(a[3:4, ])))
  expect_equal(a$x_80[2] / a$x_80[1], (1 - 20 / 1000) / (1 - 20 / 200))
  expect_true(all(a$x_peak[1:2] < a$x_50[1:2] & a$x_50[1:2] < a$x_70[1:2] & a$x_70[1:2] < a$x_80[1:2]))
  # peak: X*max = -c/b + d (Kljun et al. 2015, Eq. 20) on the discrete grid
  k <- ffp_constants()
  s <- 20 / (1 - 20 / 1000) * 3 / 0.3 * k$k
  expect_equal(a$x_peak[1] / s, -k$c / k$b + k$d, tolerance = 0.03)
  # x_R by the ETC convention is normalised by the integral over the code's X* grid (X* <= 30), which holds
  # exp(-c / (30 - d)) = 95 % of the footprint (Kljun et al. 2015, Eq. 23): Eq. 24 with R * F(30), not R
  f30 <- exp(-k$c / (k$xstar_end - k$d))
  expect_equal(a$x_80[1] / s, -k$c / log(0.8 * f30) + k$d, tolerance = 0.03)
  expect_true(a$x_80[1] / s < -k$c / log(0.8) + k$d) # shorter than the "true" 80 % distance of Eq. 24
})

test_that("ffp_smooth equals the code's padded EBImage::filter2 (FFT circular convolution)", {
  set.seed(1)
  f <- matrix(stats::runif(7 * 9), 7, 9)
  k <- matrix(c(0.05, 0.1, 0.05, 0.1, 0.4, 0.1, 0.05, 0.1, 0.05), 3)
  # filter2 emulation: circular convolution of the zero-padded field via FFT, kernel centred at (1, 1)
  filter2_circ <- function(img, ker) {
    K <- matrix(0, nrow(img), ncol(img))
    K[1:3, 1:3] <- ker
    K <- K[c(2:nrow(K), 1), c(2:ncol(K), 1)] # centre the kernel on element (1, 1)
    Re(stats::fft(stats::fft(img) * stats::fft(K), inverse = TRUE)) / length(img)
  }
  code_like <- f
  for (t in 1:2) {
    p <- matrix(0, nrow(f) + 2, ncol(f) + 2)
    p[2:(nrow(f) + 1), 2:(ncol(f) + 1)] <- code_like
    code_like <- filter2_circ(p, k)[2:(nrow(f) + 1), 2:(ncol(f) + 1)]
  }
  expect_equal(ffp_smooth(f), code_like, tolerance = 1e-12)
  expect_equal(sum(ffp_smooth(matrix(c(rep(0, 12), 1, rep(0, 12)), 5))), 1) # interior mass conserved
})

test_that("ffp_field_metrics: domain fraction, target area, source regions", {
  g <- ffp_grid(c(-2, 2, -2, 2), 1)
  f <- matrix(0, 5, 5)
  f[3, 3] <- 0.5
  f[3, 4] <- 0.3
  f[2, 3] <- 0.2
  mask <- matrix(FALSE, 5, 5)
  mask[3, 3:5] <- TRUE
  m <- ffp_field_metrics(f, g, mask, r = c(0.5, 0.8))
  expect_equal(m$metrics$f_domain, 1)
  expect_equal(m$metrics$f_target_area, 0.8)
  expect_equal(m$metrics$area_50, 1)
  expect_equal(m$metrics$area_80, 2)
  expect_equal(m$metrics$xmax_80, 1)
  expect_false(m$metrics$edge_80)
  expect_equal(m$region, sort(c(13, 18)))
  expect_equal(jaccard_index(c(1, 2, 3), c(2, 3, 4)), 0.5)
})

test_that("point_in_rings handles holes (even-odd rule)", {
  outer <- cbind(c(-2, 2, 2, -2, -2), c(-2, -2, 2, 2, -2))
  hole <- cbind(c(-1, 1, 1, -1, -1), c(-1, -1, 1, 1, -1))
  expect_equal(point_in_rings(c(0, 1.5, 3), c(0, 0, 0), list(outer, hole)), c(FALSE, TRUE, FALSE))
})

test_that("FFP code wrappers agree with the paper's equations (skipped without the Kljun code)", {
  dir <- here::here("legacy")
  skip_if_not(all(file.exists(file.path(dir, c("calc_footprint_FFP.R", "calc_footprint_FFP_climatology.R")))),
              "Kljun FFP code not available (local only, cfg$ffp_code_dir)")
  env <- load_ffp_code(dir)
  a <- ffp_ci_paper(zm = c(20, 12), h = c(800, 300), umean = c(3, 2), ustar = c(0.4, 0.3))
  for (i in 1:2) {
    b <- ffp_ci_code(env, zm = c(20, 12)[i], h = c(800, 300)[i], umean = c(3, 2)[i], L = c(-50, 80)[i],
                     sigmav = c(0.8, 0.5)[i], ustar = c(0.4, 0.3)[i], wind_dir = 200)
    expect_equal(b$flag_err, 0)
    expect_equal(unlist(a[i, ]), unlist(b[1, names(a)]), tolerance = 1e-9)
  }
  g <- ffp_grid(c(-200, 200, -200, 200), 2)
  f <- ffp_field_halfhour(env, g, zm = 20, h = 800, umean = 3, L = -50, sigmav = 0.8, ustar = 0.4, wind_dir = 200)
  expect_equal(dim(f), c(length(g$y), length(g$x)))
  expect_true(sum(f) * g$dx^2 > 0.5)
})
