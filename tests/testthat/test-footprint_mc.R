# Spec 16: footprint uncertainty from the FFP inputs (draws, sampling, validity, attribution, runner).

base_row <- function(n = 1) {
  dplyr::tibble(row_id = seq_len(n), ec_height_m = 30, canopy_height_m = 21, umean = 3, L = -50, sigmav = 0.8,
                ustar = 0.5, wind_dir = 200, h = 1200, day = TRUE, ustar_implausible = FALSE, excluded = FALSE)
}

test_that("lognormal factors: median 1 and the requested CV; +-1 sigma quantiles; truncated hc", {
  set.seed(1)
  for (cv in c(0.05, 0.15, 0.5)) {
    f <- mc_factor(stats::rnorm(4e5), "lognormal", cv)
    expect_equal(stats::median(f), 1, tolerance = 0.01)
    expect_equal(stats::sd(f) / mean(f), cv, tolerance = 0.02)
  }
  oat <- mc_oat_runs()
  expect_equal(nrow(oat), 21)
  expect_equal(sum(oat$validity == "etc"), 17)
  expect_setequal(oat$run[oat$validity == "rsl_off"], c("hc_minus_rsl_off", "hc_plus_rsl_off", "z0f_minus_rsl_off",
                                                        "z0f_plus_rsl_off"))
  expect_equal(oat$z[oat$run == "hc_plus_rsl_off"], oat$z[oat$run == "hc_plus"])
  expect_equal(mc_oat_ffp_runs(oat), 19L)
  expect_equal(oat$z[oat$run == "sigmav_minus"], stats::qnorm(0.1587))
  expect_equal(mc_factor(oat$z[oat$run == "sigmav_plus"], "lognormal", 0.15), exp(mc_sigma_ln(0.15) * stats::qnorm(0.8413)))
  # hc: normal truncated at +-3 sigma
  set.seed(2)
  z <- mc_draw_z(2e4)
  expect_true(all(abs(z[, "hc"]) <= 3))
  expect_lt(abs(mc_qz(0.8413, trunc = TRUE) - 1), 0.01)
  expect_equal(mc_factor(1, "normal_trunc", 0.1), 1.1)
})

test_that("wind direction: additive and wrapped to [0, 360); heights follow hc", {
  b <- base_row(3)
  b$wind_dir <- c(355, 2, 180)
  p <- mc_perturb(b, list(wd = c(10, -7.5, 0)))
  expect_equal(p$wind_dir, c(5, 354.5, 180))
  p2 <- mc_perturb(base_row(), list(hc = 1.1, z0f = 2))
  expect_equal(p2$d_m, 2 / 3 * 23.1)
  expect_equal(p2$zm_m, 30 - 2 / 3 * 23.1)
  expect_equal(p2$z0_m, 0.15 * 23.1 * 2)
})

test_that("L: alone via 1/L in design A; linked to u* in the joint draws (L'/L0 = (u*'/u*0)^3 when eps_L = 1)", {
  b <- base_row()
  pa <- mc_perturb(b, list(ustar = 1.2), link_L = FALSE)
  expect_equal(pa$L, b$L)
  pb <- mc_perturb(b, list(ustar = 1.2, L = 1), link_L = TRUE)
  expect_equal(pb$L / b$L, 1.2^3)
  pc <- mc_perturb(b, list(L = 1.25))
  expect_equal(pc$L, b$L / 1.25) # 1/L times 1.25, sign kept
  expect_true(sign(pc$L) == sign(b$L))
})

test_that("stratified sample and weights: n per stratum, small strata kept whole, weights sum to 1", {
  st <- c(rep("unstable_day", 300), rep("neutral_day", 5), rep("stable_night", 50), NA)
  set.seed(3)
  s <- mc_sample(st, 20)
  expect_equal(as.integer(table(s$stratum)[c("unstable_day", "neutral_day", "stable_night")]), c(20L, 5L, 20L))
  expect_equal(sum(s$w_row), 1)
  expect_equal(unique(s$w_stratum[s$stratum == "neutral_day"]), 5 / 355)
  expect_false(anyNA(st[s$i]))
  s$in_b <- mc_subsample(s, 3)
  expect_equal(as.integer(table(s$stratum[s$in_b])[c("unstable_day", "neutral_day", "stable_night")]), c(3L, 3L, 3L))
  w <- mc_b_weights(s)
  expect_equal(sum(w, na.rm = TRUE), 1)
  expect_equal(mc_stratum(c(10, 10, 10, NA), c(-50, 200, 50, 1), c(TRUE, FALSE, NA, TRUE)),
               c("unstable_day", "neutral_night", NA, NA))
  expect_equal(mc_site_seed("BE-Bra"), mc_site_seed("BE-Bra"))
  expect_false(mc_site_seed("BE-Bra") == mc_site_seed("FI-Hyy"))
})

test_that("validity under perturbation: zm/L, u* and roughness boundaries of the etc set", {
  b <- base_row(4)
  m_etc <- 5
  # zm = 30 - 14 = 16 m; z0 = 3.15 m -> m z0 = 15.75 < 16 (valid); z0 factor 1.1 -> 17.3 > 16 (fails e_rsl)
  v <- mc_validity(mc_perturb(b, list(z0f = c(1, 1.1, 1, 1), ustar = c(1, 1, 0.39, 1), L = c(1, 1, 1, -1000))), m_etc)
  expect_equal(v$valid_etc, c(TRUE, FALSE, FALSE, FALSE))
  # u* 0.5 x 0.39 = 0.195 < 0.2 (e_ustar); a factor -1000 on 1/L flips the sign: L = 0.05 m, zm/L = 320 > 1 (e_stab)
  expect_equal(v$fail, c(NA, "e_rsl", "e_ustar", "e_stab"))
  # without e_rsl (rsl_off): only the roughness failure becomes valid; the other conditions are kept
  expect_equal(v$valid_norsl, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(v$rsl_fail, c(FALSE, TRUE, FALSE, FALSE))
})

test_that("rsl_off runs: hc and z0 valid without e_rsl, same FFP inputs as the etc runs, own effects and flip share", {
  b <- base_row(1) # zm 16 m, z0 3.15 m: zm / z0 = 5.08, just above m = 5
  a <- mc_run_inputs(b, "A", m_etc = 5)
  expect_equal(nrow(a), 21)
  g <- function(r, k) a[[k]][a$run == r]
  expect_false(g("hc_plus", "valid_run")) # +1 sigma hc: zm / z0 = 4.2 < 5
  expect_true(g("hc_plus", "rsl_fail"))
  expect_true(g("hc_plus_rsl_off", "valid_run"))
  expect_false(g("hc_plus_rsl_off", "valid_etc"))
  expect_true(g("hc_minus", "valid_run"))
  expect_false(g("z0f_plus", "valid_run"))
  expect_true(g("z0f_plus_rsl_off", "valid_run"))
  k <- mc_run_key(a)
  expect_equal(k[a$run == "hc_plus_rsl_off"], k[a$run == "hc_plus"])
  expect_equal(k[a$run == "z0f_plus_rsl_off"], k[a$run == "base"]) # z0 enters validity only
  # synthetic results: TA +1 pp per sigma of hc where the run is valid under its set
  s <- b
  s$stratum <- "unstable_day"; s$w_stratum <- 1; s$w_row <- 1; s$in_b <- FALSE; s$valid_qc <- TRUE
  s$valid_qc_used <- TRUE; s$timestamp_start_utc <- as.POSIXct("2025-01-01", tz = "UTC")
  s$ta__pblhmod <- 0.9; s$x80__pblhmod <- 200; s$ta__era5 <- 0.9; s$x80__era5 <- 200
  s$ta__pblhmod_clim <- NA_real_; s$x80__pblhmod_clim <- NA_real_
  a$ffp_ok <- a$valid_run
  a$f_target_area <- ifelse(a$valid_run, 0.9 + ifelse(a$input %in% "hc", 0.01 * a$side, 0), NA)
  a$x_80 <- ifelse(a$valid_run, 200, NA)
  oat <- mc_oat_table(mc_join_runs(list(sample = s), a))
  eff <- mc_oat_effects(oat)
  expect_true(is.na(eff$eff_ta_pp[eff$input == "hc" & eff$validity == "etc"])) # one side flipped
  expect_equal(eff$eff_ta_pp[eff$input == "hc" & eff$validity == "rsl_off"], 1)
  expect_equal(eff$eff_ta_pp[eff$input == "z0f" & eff$validity == "rsl_off"], 0)
  expect_equal(eff$rsl_flip[eff$input == "hc"], c(0.5, 0.5))
  expect_equal(eff$rsl_flip[eff$input == "umean"], 0)
  src <- mc_source_effects(s)
  sm <- mc_site_summary("XX-Tst", eff, src, dplyr::tibble())
  r <- sm[sm$item == "hc|rsl_off" & sm$quantity == "ta_pp" & sm$stratum == "all", ]
  expect_equal(r$median, 1)
  expect_equal(r$share_rsl_flip, 0.5)
  net <- mc_network_summary(sm)
  expect_equal(net$share_rsl_flip[net$item == "z0f|rsl_off" & net$quantity == "ta_pp" & net$stratum == "all"], 0.5)
  p <- plot_mc_input_ranking(dplyr::bind_rows(sm, net), "ta_pp")
  expect_s3_class(p, "ggplot")
  lab <- ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$get_labels()
  expect_true(any(grepl("rsl test off", lab)))
})

test_that("validity: zm/L at the stable (1) and unstable (-15.5) limits", {
  b <- base_row(2)
  b$L <- c(16 / 1.0, -16 / 15.4) # zm/L = 1 (valid: <= 1), -15.4 (valid: > -15.5)
  v <- mc_validity(mc_perturb(b, list()), 5)
  expect_equal(v$valid_etc, c(TRUE, TRUE))
  v2 <- mc_validity(mc_perturb(b, list(L = c(1.01, 1.01))), 5) # 1/L x 1.01: zm/L = 1.01, -15.55
  expect_equal(v2$valid_etc, c(FALSE, FALSE))
  expect_equal(v2$fail, c("e_stab", "code_stab"))
})

test_that("run tables: 21 runs per row in A (17 etc + 4 rsl_off), draws x h options in B, identical inputs keyed once", {
  b <- base_row(2)
  a <- mc_run_inputs(b, "A", m_etc = 2)
  expect_equal(nrow(a), 42)
  expect_equal(sum(a$run == "base"), 2)
  z <- list(`1` = mc_draw_z(4), `2` = mc_draw_z(4))
  b$day <- c(TRUE, FALSE)
  bb <- mc_run_inputs(b, "B", draws = list(z = z, h_options = list(h_pct30 = list(day = 0.3, night = 0.3),
                                                                    h_empirical = list(day = 0.5, night = 0.3))), m_etc = 2)
  expect_equal(nrow(bb), 2 * 4 * 2)
  k <- mc_run_key(bb)
  # night row: h_empirical equals h_pct30 (same deviates, same sd) -> same key; day row differs
  n1 <- bb[bb$row_id == 2, ]
  expect_equal(k[bb$row_id == 2 & bb$h_option == "h_pct30"], k[bb$row_id == 2 & bb$h_option == "h_empirical"])
  expect_false(any(k[bb$row_id == 1 & bb$h_option == "h_pct30"] == k[bb$row_id == 1 & bb$h_option == "h_empirical"]))
  expect_true(all(c("z_sigmav", "z_h") %in% names(bb)))
})

test_that("attribution recovers the shares of a synthetic linear response", {
  set.seed(4)
  n_row <- 30
  n_draw <- 200
  z <- matrix(stats::rnorm(n_row * n_draw * 3), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
  row <- rep(seq_len(n_row), each = n_draw)
  y <- 10 * row + 3 * z[, "a"] + 1 * z[, "b"] + 0 * z[, "c"] # row offsets removed by the method
  at <- mc_attribution(y, z, row)
  expect_equal(at$share, c(0.9, 0.1, 0), tolerance = 0.02)
  expect_equal(at$r2[1], 1, tolerance = 1e-8)
  y2 <- y + stats::rnorm(length(y), sd = sqrt(10)) # half the within-row variance unexplained
  at2 <- mc_attribution(y2, z, row)
  expect_equal(at2$r2[1], 0.5, tolerance = 0.03)
  expect_equal(at2$share, c(0.45, 0.05, 0), tolerance = 0.03)
})

test_that("half-range, weighted quantiles, chunk list", {
  expect_equal(mc_half_range(c(-1, NA), c(3, 1)), c(2, NA))
  expect_equal(mc_wquantile(c(1, 2, 3), c(0.1, 0.1, 0.8), 0.5), 3)
  expect_equal(mc_wquantile(c(1, 2, 3), c(1, 1, 1), 0.5), 2)
  cl <- mc_chunk_list("XX-Tst", 1020, 300, runs_a = 17, runs_b = 100, sec_per_run = 3, cpus = 2, target_min = 120)
  expect_equal(sum(cl$n_rows[cl$design == "A"]), 1020)
  expect_equal(sum(cl$n_rows[cl$design == "B"]), 300)
  expect_true(all(cl$est_min <= 120))
  expect_equal(cl$first_row[2], cl$last_row[1] + 1L)
})

test_that("design A and B tables from synthetic runs: base difference, half-range, quantiles, ICOS probability", {
  s <- base_row(2)
  s$stratum <- c("unstable_day", "stable_night")
  s$w_stratum <- c(0.6, 0.4); s$w_row <- c(0.6, 0.4); s$in_b <- TRUE
  s$valid_qc <- TRUE; s$valid_qc_used <- TRUE; s$timestamp_start_utc <- as.POSIXct("2025-01-01", tz = "UTC")
  s$ta__pblhmod <- c(0.95, 0.80); s$ta__era5 <- c(0.94, 0.70); s$ta__pblhmod_clim <- c(0.95, 0.79)
  s$x80__pblhmod <- c(100, 300); s$x80__era5 <- c(100, 330); s$x80__pblhmod_clim <- c(100, 300)
  a <- mc_run_inputs(s, "A", m_etc = 2)
  a$ffp_ok <- TRUE
  a$f_target_area <- ifelse(a$run == "base", s$ta__pblhmod[a$row_id], s$ta__pblhmod[a$row_id] + 0.01 * a$side)
  a$x_80 <- ifelse(a$run == "base", s$x80__pblhmod[a$row_id], s$x80__pblhmod[a$row_id] * exp(0.1 * a$side))
  draws <- list(sample = s)
  x <- mc_join_runs(draws, a)
  oat <- mc_oat_table(x)
  expect_equal(max(abs(oat$base_diff_pp)), 0)
  eff <- mc_oat_effects(oat)
  expect_equal(unique(round(eff$eff_ta_pp, 10)), 1)
  expect_equal(unique(round(eff$eff_lnx80, 10)), 0.1)
  src <- mc_source_effects(s)
  expect_equal(src$src_ta_pp__era5_vs_pblhmod, c(-1, -10))
  # B: 10 draws per row, TA = base + draw deviations
  zb <- list(`1` = mc_draw_z(10), `2` = mc_draw_z(10))
  bb <- mc_run_inputs(s, "B", draws = list(z = zb, h_options = list(h_pct30 = list(day = 0.3, night = 0.3))), m_etc = 2)
  bb$valid_etc <- TRUE; bb$ffp_ok <- TRUE
  bb$f_target_area <- s$ta__pblhmod[bb$row_id] + 0.001 * (bb$draw - 5.5)
  bb$x_80 <- s$x80__pblhmod[bb$row_id]
  x2 <- mc_join_runs(draws, dplyr::bind_rows(a, bb))
  mc <- mc_mc_table(x2)
  expect_equal(nrow(mc), 2)
  expect_equal(mc$n_valid, c(10L, 10L))
  expect_gt(min(mc$width_ta_pp), 0)
  ic <- mc_icos_table(x2, src, "XX-Tst")
  r <- ic[ic$case == "pblhmod" & ic$criterion == "rebmann_ta90", ]
  expect_equal(r$share_base, 0.6) # only row 1 (weight 0.6) has TA >= 90 %
  expect_equal(r$prob_pass, 1)
  e <- ic[ic$case == "era5" & ic$criterion == "icos_lab_ta70", ]
  expect_equal(e$share_base, 1) # 94 % and 70 % both >= 70 %
  sm <- mc_site_summary("XX-Tst", eff, src, mc)
  expect_true(all(c("oat_1sigma", "source_effect", "joint_width", "ratio") %in% sm$kind))
  net <- mc_network_summary(sm)
  expect_true(all(net$scope == "network"))
})

test_that("base run equals the 08b computation of the same half-hour (small grid; FFP code local only)", {
  cfg <- get_cfg()
  ffp_dir <- tryCatch(resolve_path(cfg$ffp_code_dir), error = function(e) NA_character_)
  skip_if(is.na(ffp_dir) || !file.exists(file.path(ffp_dir, "calc_footprint_FFP.R")), "Kljun FFP code not on this machine")
  grid <- ffp_grid(c(-200, 200, -200, 200), 2)
  mask <- grid$rho < 150 & outer(grid$y, grid$x, function(y, x) x > -50)
  b <- base_row(1)
  runs <- mc_run_inputs(b, "A", m_etc = 2)
  runs <- runs[runs$run %in% c("base", "sigmav_plus", "wd_minus"), ]
  res <- mc_run_ffp(runs, ffp_dir, grid, mask, c(0.5, 0.7, 0.8, 0.9), workers = 1)
  env <- load_ffp_code(ffp_dir)
  ref <- ffp_halfhour_run(env, grid, mask, c(0.5, 0.7, 0.8, 0.9), runs$zm_m[1], runs$h[1], runs$umean[1], runs$L[1],
                          runs$sigmav[1], runs$ustar[1], runs$wind_dir[1])$metrics
  expect_equal(res$f_target_area[1], ref$f_target_area)
  expect_equal(res$x_80[1], ffp_ci_paper(runs$zm_m[1], runs$h[1], runs$umean[1], runs$ustar[1])$x_80)
  expect_false(isTRUE(all.equal(res$f_target_area[1], res$f_target_area[3]))) # wind direction moves the footprint
  expect_equal(res$x_80[2], res$x_80[1]) # sigma_v does not change the crosswind-integrated distances
  expect_equal(attr(res, "n_ffp"), 3L)
  # an invalid etc run listed before its valid rsl_off twin: the twin is computed (m = 5: hc_plus fails e_rsl)
  r2 <- mc_run_inputs(b, "A", m_etc = 5)
  r2 <- r2[r2$run %in% c("hc_plus", "hc_plus_rsl_off", "z0f_plus", "z0f_plus_rsl_off"), ]
  res2 <- mc_run_ffp(r2, ffp_dir, grid, mask, c(0.5, 0.7, 0.8, 0.9), workers = 1)
  expect_equal(res2$ffp_ok, c(FALSE, FALSE, TRUE, TRUE)) # order of mc_oat_runs(): etc runs first
  expect_equal(attr(res2, "n_ffp"), 2L)
})
