# Spec 18: PBLH source and EC fluxes (R/pblh_ec_fluxes.R, R/pblh_ec_collect.R)

hh <- function(start, n) as.POSIXct(start, tz = "UTC") + 1800 * (seq_len(n) - 1)

test_that("night grouping: one id across midnight, polar night split per solar day, polar day none", {
  t <- hh("2024-06-01 12:00", 48)
  night <- as.integer(format(t, "%H", tz = "UTC")) >= 21 | as.integer(format(t, "%H", tz = "UTC")) < 4
  id <- night_groups(t, night, lon = 0)
  expect_true(all(is.na(id[!night])))
  expect_length(unique(stats::na.omit(id)), 1)                  # 21:00 - 04:00 is one night
  expect_equal(unique(stats::na.omit(id)), as.integer(as.Date("2024-06-01")))
  tp <- hh("2024-12-20 00:00", 96)                              # polar night: 48 h of night
  idp <- night_groups(tp, rep(TRUE, 96), lon = 0)
  expect_length(unique(idp), 3)                                 # cut at each solar noon
  expect_true(all(diff(idp) >= 0))
  expect_true(all(is.na(night_groups(tp, rep(FALSE, 96), lon = 0))))  # polar day
  # longitude shifts the cut: at 180 E local noon is 00 UTC
  expect_false(identical(night_groups(tp, rep(TRUE, 96), 180), idp))
})

test_that("u* bins and the u*-filter agreement (QC flag vs threshold) on a hand case", {
  expect_equal(ustar_ratio_bin(c(0.9, 1, 1.49, 1.5, 2, 5, NA)), c(NA, "1-1.5", "1-1.5", "1.5-2", ">2", ">2", NA))
  a <- ustar_filter_agreement(c(TRUE, TRUE, FALSE, FALSE, NA, TRUE), c(1, 0, 2, 0, 0, NA))
  expect_equal(a$n, 4L)
  expect_equal(c(a$n_both, a$n_thr_only, a$n_qc_only, a$n_neither), c(1L, 1L, 1L, 1L))
  expect_equal(a$agreement, 0.5)
})

test_that("SC = NEE - FC check and profile detection", {
  d <- dplyr::tibble(FC = c(1, 2, 3, NA), NEE = c(1.5, 2.01, 4, 1), SC = c(0.5, 0, 0.2, 0.1), SC_FLAG = c(0, 0, 1, 0))
  x <- storage_profile_info(c("FC", "CO2_DRY_1", "CO2_DRY_2", "CO2_DRY", "SC"), d)
  expect_equal(x$n_profile_levels, 2L)                          # CO2_DRY alone is not a level
  expect_equal(x$sc_check_n, 3L)
  expect_equal(x$sc_check_share, 2 / 3)                         # 1.5-1-0.5 = 0; 2.01-2-0 = 0.01 <= tol; 4-3-0.2 off
  expect_equal(c(x$n_sc_profile, x$n_sc_flag1), c(3L, 1L))
  expect_true(x$storage)
  expect_false(storage_profile_info(c("FC", "SC"), d)$storage)  # no profile levels
  expect_false(storage_profile_info(character(), d[c("FC", "NEE")])$has_sc)
})

mk_rows <- function(n = 8) {
  dplyr::tibble(timestamp_start_utc = hh("2024-07-01 00:00", n), night = TRUE, zl = "stable", season = "JJA",
                zm_m = 10, h_era5_blh_center = c(50, 50, 200, 200, 400, 400, 80, NA)[seq_len(n)],
                h_pblhmod_default = c(50, 20, 20, 300, 300, 300, 300, 300)[seq_len(n)],
                h_pblhmod_clim = 300, h_meas = NA_real_, meas_kind = NA_character_)
}

test_that("section 1: zm/h shares on a synthetic series", {
  r <- mk_rows()
  s <- sl_share_table(r, "XX-Tst")
  g <- s[s$daynight == "night" & s$zl == "stable" & s$season == "all", ]
  e <- g[g$source == "era5_blh_center" & g$row_basis == "own", ]
  expect_equal(e$n, 7L)                                         # NA h dropped
  expect_equal(e$share_gt_0.1, 3 / 7)                           # zm/h = 0.2, 0.2, 0.125 > 0.1
  expect_equal(e$share_gt_0.05, 3 / 7)                          # 0.05 is not > 0.05
  expect_equal(e$share_gt_0.2, 0)
  c1 <- g[g$source == "pblhmod_default" & g$row_basis == "common", ]
  expect_equal(c1$n, 7L)                                        # common = rows with all three h
  expect_equal(c1$share_gt_0.1, 3 / 7)                          # 0.2, 0.5, 0.5
  expect_equal(nrow(s[s$daynight == "day" & s$zl == "all" & s$season == "all" & s$n > 0, ]), 0)
  # measured: censored rows counted apart, exact values only in the shares
  r$h_meas <- c(NA, 50, NA, 200, NA, NA, NA, NA)
  r$meas_kind <- c("below", "exact", "lower_bound", "exact", NA, NA, NA, NA)
  m <- sl_share_table(r, "XX-Tst", meas_src = "cloudnet_eps_mlh")
  m <- m[m$source == "cloudnet_eps_mlh" & m$daynight == "all" & m$zl == "all" & m$season == "all", ]
  expect_equal(c(m$n, m$n_censored_below, m$n_lower_bound), c(2L, 1L, 1L))
  expect_equal(m$share_gt_0.1, 0.5)
})

test_that("section 2: 2 x 2 table of u*-filtered x above the surface layer", {
  r <- mk_rows()
  r$below_thr <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE)
  r$nee_qc <- c(1, 1, 0, 0, 1, 0, 0, 1)
  u <- ustar_cross_table(r, "XX-Tst")
  e <- u[u$source == "era5_blh_center" & u$filter_def == "threshold", ]
  # rows 1-7 (row 8 has no ERA5 h); above (zm/h > 0.1) = rows 1, 2, 7
  expect_equal(c(e$n, e$n_filt_above, e$n_filt_below, e$n_pass_above, e$n_pass_below), c(7L, 2L, 1L, 1L, 3L))
  expect_equal(e$share_above_removed, 2 / 3)
  expect_equal(e$share_pass_above, 1 / 4)
  q <- u[u$source == "era5_blh_center" & u$filter_def == "qc", ]
  expect_equal(q$n_filt_above, 2L)
})

test_that("night-block bootstrap resamples whole nights", {
  night <- rep(1:20, each = 5)
  v <- rep(1:20, each = 5)                                      # constant within a night
  b <- night_block_boot(night, function(i) mean(v[i]), n_boot = 300, seed = 1)
  expect_equal(b$est, 10.5)
  expect_true(b$lo < 10.5 && b$hi > 10.5)
  # every replicate mean is a mean of whole nights: 5 rows each, so 100 rows per replicate
  b2 <- night_block_boot(night, function(i) length(i), n_boot = 50, seed = 1)
  expect_true(all(b2$reps == 100))
  expect_equal(night_block_boot(night, function(i) mean(v[i]), 300, 1)$lo, b$lo)  # seeded
  b1 <- night_block_boot(rep(1, 10), function(i) mean(i), n_boot = 10)
  expect_true(is.na(b1$lo))                                     # one night: no interval
})

test_that("section 3: residual difference above minus within the surface layer, per cell", {
  set.seed(3)
  n <- 4000
  r <- dplyr::tibble(night = TRUE, night_id = rep(1:200, each = 20), below_thr = FALSE, nee_qc = 0,
                     reco_nt = 3, ta_f_k = stats::runif(n, 280, 290), zm_m = 10,
                     ustar_ratio = stats::runif(n, 1, 3), sc_flag = 0,
                     h_era5_blh_center = ifelse(stats::runif(n) < 0.5, 50, 500), h_pblhmod_default = 500)
  r$nee_vut <- r$reco_nt + stats::rnorm(n, 0, 0.1) + ifelse(r$h_era5_blh_center == 50, 1, 0)  # +1 above the SL
  r$fc <- r$nee_vut
  x <- nee_residual_table(r, "XX-Tst", n_boot = 100)
  a <- x[x$source == "era5_blh_center" & x$flux == "NEE" & x$measure == "abs" & x$stratum == "all", ]
  expect_equal(a$diff, 1, tolerance = 0.05)
  expect_true(a$lo90 < a$diff && a$hi90 > a$diff)
  rel <- x[x$source == "era5_blh_center" & x$flux == "NEE" & x$measure == "rel" & x$stratum == "all", ]
  expect_equal(rel$diff, 1 / 3, tolerance = 0.03)
  p <- x[x$source == "pblhmod_default" & x$stratum == "all", ]
  expect_true(all(!p$valid))                                    # nothing above the SL with PBLHmod
  expect_equal(nrow(x[x$stratum == "cell" & x$source == "era5_blh_center" & x$flux == "NEE" & x$measure == "abs", ]), 9L)
  small <- nee_residual_table(r[r$night_id <= 3, ], "XX-Tst", n_boot = 20)
  expect_true(all(!small$valid[small$stratum == "cell"]))       # < 30 rows per side
  expect_equal(nrow(nee_residual_table(dplyr::mutate(r, reco_nt = NA_real_), "XX-Tst", n_boot = 5)), 0L)
})

test_that("section 4: night accumulation, morning release and growth of h", {
  t <- hh("2024-07-01 12:00", 48)
  hr <- as.integer(format(t, "%H", tz = "UTC"))
  night <- hr >= 20 | hr < 4
  r <- dplyr::tibble(timestamp_start_utc = t, night = night, night_id = night_groups(t, night, 0), ref_year = TRUE,
                     sc = ifelse(night, 1, -2), sc_flag = 0, zm_m = 10,
                     h_era5_blh_center = ifelse(night, 100, 100 + 50 * (seq_along(t) - 33)),
                     h_pblhmod_default = 200)
  ns <- storage_nights(r, "XX-Tst")
  expect_equal(nrow(ns), 1L)
  expect_equal(ns$accum, 16 * 1800)                             # 16 night half-hours of SC = 1
  expect_equal(ns$release, -2 * 8 * 1800)
  expect_equal(ns$dh_era5_blh_center, 50 * 7)                   # from the first day row (index 33) to index 40
  expect_equal(ns$h_night_pblhmod_default, 200)
  r$sc[which(night)[1:5]] <- NA                                 # coverage 11/16 < 0.8
  expect_true(is.na(storage_nights(r, "XX-Tst")$accum))
})

test_that("section 5: Panofsky prediction on hand values and the within-site ranking", {
  expect_equal(panofsky_sv_ustar(1000, -100), 17^(1 / 3))
  expect_equal(panofsky_sv_ustar(c(1000, 1000, -5), c(50, 0, -10)), rep(NA_real_, 3))
  set.seed(5)
  n <- 400
  h_true <- stats::runif(n, 500, 2000)
  L <- -stats::runif(n, 10, 150)
  r <- dplyr::tibble(sw_in_pot = 500, zm_m = 5, L = L, ustar = 0.4,
                     h_era5_blh_center = h_true, h_pblhmod_default = sample(h_true),
                     h_pblhmod_clim = h_true * 1.5, h_meas = NA_real_)
  r$v_sigma <- r$ustar * panofsky_sv_ustar(h_true, L) * exp(stats::rnorm(n, 0, 0.01))
  p <- panofsky_table(r, "XX-Tst")$stats
  cm <- p[p$row_basis == "common_models", ]
  expect_equal(cm$rank[cm$source == "era5_blh_center"], 1L)      # the h that generated sigma_v
  expect_equal(cm$median_log_ratio[cm$source == "era5_blh_center"], 0, tolerance = 0.01)
  expect_equal(cm$slope[cm$source == "era5_blh_center"], 1 / 3, tolerance = 0.02)  # Panofsky exponent
  expect_lt(cm$median_abs_dev[cm$source == "era5_blh_center"], cm$median_abs_dev[cm$source == "pblhmod_default"])
  expect_true(all(p$row_basis != "common_with_measured"))
})

test_that("section 6 and missing variables: a part is skipped, a missing SCF column is dropped", {
  expect_equal(pec_missing(list(a = TRUE, b = FALSE, c = NA)), c("b", "c"))
  expect_length(pec_missing(list(a = TRUE)), 0)
  r <- dplyr::tibble(zl = "unstable", L = -50, zm_m = 5, h_era5_blh_center = c(500, 1000), h_pblhmod_default = 800,
                     h_scf = 1.02, le_scf = 1.1, fc_scf = 1.05, tau_scf = NA_real_)
  s <- scf_h_table(r, "XX-Tst")
  expect_false("TAU_SCF_STAT" %in% s$scf_var)
  expect_true(pec_has(r, "h_scf"))
  expect_false(pec_has(r, "tau_scf"))
  expect_false(pec_has(r, "absent"))
})

test_that("collection: network medians, stable-night flags and Panofsky ranks", {
  sl <- dplyr::tibble(site_id = c("AA-One", "BB-Two"), source = "era5_blh_center", row_basis = "common",
                      daynight = "night", zl = "stable", season = "all", share_gt_0.1 = c(0.1, 0.5))
  net <- pec_network_median(sl, c("source", "row_basis", "daynight", "zl", "season"), "share_gt_0.1")
  expect_equal(c(net$site_id, net$n_sites, net$share_gt_0.1), c("network", "2", "0.3"))
  f <- pec_stable_night_flags(sl)
  expect_equal(f$site_id, "BB-Two")
  pan <- dplyr::tibble(site_id = c("AA-One", "AA-One", "BB-Two", "BB-Two"), row_basis = "common_models",
                       source = rep(c("era5_blh_center", "pblhmod_default"), 2), rank = c(1L, 2L, 2L, 1L))
  rk <- pec_panofsky_rank(pan)
  expect_equal(rk$n_first, c(1L, 1L))
  expect_equal(rk$share_first, c(0.5, 0.5))
})
