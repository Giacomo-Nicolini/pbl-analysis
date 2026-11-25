# Spec 19: circular statistics, row selection (= spec-15 common_pe), descriptor layout, row-level R2, plot helpers.

test_that("circular mean, difference and reversal handle the 0/360 wrap", {
  expect_equal(circ_mean_deg(c(350, 10)), 0, tolerance = 1e-9)
  expect_equal(circ_mean_deg(c(80, 100)), 90, tolerance = 1e-9)
  expect_equal(circ_mean_deg(c(340, 350, 360, 10, 20)), 0, tolerance = 1e-9)
  expect_equal(circ_mean_deg(c(170, 190)), 180, tolerance = 1e-9)
  expect_true(is.na(circ_mean_deg(c(0, 180))))
  expect_true(is.na(circ_mean_deg(c(NA, NaN))))
  expect_equal(circ_diff_deg(10, 350), 20)
  expect_equal(circ_diff_deg(350, 10), -20)
  expect_equal(circ_diff_deg(180, 0), 180)
  expect_equal(circ_diff_deg(0, 180), 180)
  expect_equal(circ_diff_deg(c(270, 90), c(90, 270)), c(180, 180))
  # night from 10 deg, day from 190 deg: reversal 180 deg
  expect_equal(abs(circ_diff_deg(circ_mean_deg(c(0, 20)), circ_mean_deg(c(180, 200)))), 180, tolerance = 1e-9)
})

test_that("steadiness and dominant sector on hand cases", {
  expect_equal(wind_steadiness(c(90, 90, 90), c(1, 2, 3)), 1)
  expect_equal(wind_steadiness(c(0, 180), c(2, 2)), 0, tolerance = 1e-12)
  # two equal speeds 90 deg apart: |mean vector| = sqrt(2)/2 of the speed
  expect_equal(wind_steadiness(c(0, 90), c(1, 1)), sqrt(2) / 2, tolerance = 1e-12)
  expect_equal(wind_steadiness(c(355, 5), c(1, 1)), cos(5 * pi / 180), tolerance = 1e-12)
  expect_true(is.na(wind_steadiness(numeric(), numeric())))
  ds <- dominant_sector(c(350, 355, 5, 10, 180))
  expect_equal(ds$share, 0.8)
  expect_true(abs(circ_diff_deg(ds$centre, 0)) <= 35)
  ds2 <- dominant_sector(c(rep(270, 6), rep(90, 4)))
  expect_equal(ds2$share, 0.6)
  expect_equal(abs(circ_diff_deg(ds2$centre, 270)) <= 45, TRUE)
})

# Synthetic site: input table, subsets and chunks of one month
make_alpine_site <- function(dir, site = "XX-Tst", n = 400, seed = 1) {
  set.seed(seed)
  t <- as.POSIXct("2024-03-01", tz = "UTC") + 1800 * (seq_len(n) - 1)
  tab <- dplyr::tibble(timestamp_start_utc = t, ec_height_m = 30, canopy_height_m = 15, zm_m = 20,
                       umean = runif(n, 0.5, 4), L = sample(c(-50, -500, 500, 40, 15), n, TRUE),
                       ustar = runif(n, 0.05, 0.8), wind_dir = runif(n, 0, 360),
                       day = rep(c(TRUE, FALSE), each = 12, length.out = n), regime = "stable",
                       h_pblhmod = runif(n, 80, 1500), h_pblhmod_clim = runif(n, 80, 1500), h_era5 = runif(n, 50, 1500))
  sub <- list()
  for (cs in c("pblhmod", "pblhmod_clim", "era5")) for (u in c("valid", "valid_qc", "valid_qc_used")) {
    sub[[paste(u, "etc", cs, sep = "__")]] <- runif(n) > (if (u == "valid") 0.1 else 0.3)
  }
  sd <- file.path(dir, site)
  dir.create(file.path(sd, "chunks"), recursive = TRUE)
  saveRDS(list(site_id = site, tab = tab, subsets = dplyr::as_tibble(sub)),
          file.path(sd, paste0(site, "_footprint_inputs.rds")))
  for (cs in c("pblhmod", "pblhmod_clim", "era5")) {
    ok <- runif(n) > 0.05
    readr::write_csv(dplyr::tibble(timestamp_start_utc = t, h = tab[[paste0("h_", cs)]], x_80 = runif(n, 100, 2000),
                                   ffp_ok = ok, f_target_area = ifelse(ok, runif(n, 0.5, 1), NA)),
                     file.path(sd, "chunks", paste0(site, "_", cs, "_2024-03_hh.csv.gz")))
    saveRDS(list(), file.path(sd, "chunks", paste0(site, "_", cs, "_2024-03_meta.rds")))
  }
  site
}

test_that("row selection equals the spec-15 common_pe rows and its target-area statistics", {
  d <- withr::local_tempdir()
  s <- make_alpine_site(d)
  r <- fp_site_rows(s, 2024, d, sets = "etc", subsets = c("valid", "valid_qc_used"), tab_cols = alpine_tab_cols())
  expect_true(all(alpine_tab_cols() %in% names(r)))
  r <- alpine_row_vars(r, w_sigma = rep(0.4, nrow(r)), ustar_thr = rep(0.2, nrow(r)))
  b <- alpine_basis(r, "etc", "valid_qc_used", require_ffp = TRUE)
  expect_identical(b, fp_common(r, "etc", "valid_qc_used"))
  expect_gt(sum(b), 50)
  # flags only: ignores the FFP result, so it is a superset of common_pe
  bf <- alpine_basis(r, "etc", "valid_qc_used", require_ffp = FALSE)
  expect_true(all(bf[b]))
  expect_gt(sum(bf), sum(b))
  desc <- alpine_site_descriptors(r, list(valid_qc_used = b, valid = alpine_basis(r, "etc", "valid")), "common_pe",
                                  list(cells = NULL, igbp = "ENF", ustar_thr_source = "VUT"))
  ta <- fp_target_area_table(r, "etc", "valid_qc_used")
  ta <- ta[ta$scope == "site" & ta$pair == "era5_vs_pblhmod", ]
  for (k in seq_len(nrow(ta))) {
    g <- function(v) desc$value[desc$variable == v & desc$subset %in% "valid_qc_used" &
                                  desc$regime == as.character(ta$regime[k]) & desc$daynight == as.character(ta$daynight[k])]
    expect_equal(g("ta_diff_n"), ta$n[k])
    expect_equal(g("ta_diff_p10"), ta$p10_pp[k], tolerance = 1e-9)
    expect_equal(g("ta_diff_p90"), ta$p90_pp[k], tolerance = 1e-9)
  }
  # spec-15 exception status from the same rule
  ex <- fp_exceptions_table(fp_target_area_table(r, "etc", "valid_qc_used"), r)
  ex <- ex[ex$pair == "era5_vs_pblhmod", ]
  st <- desc[desc$variable == "spec15_status" & desc$subset %in% "valid_qc_used" &
               desc$regime %in% c("unstable", "neutral", "stable") & desc$daynight %in% c("day", "night"), ]
  for (k in seq_len(nrow(ex))) {
    expect_equal(st$value_chr[st$regime == ex$regime[k] & st$daynight == ex$daynight[k]], ex$status[k])
  }
  expect_true(all(st$value_chr[!paste(st$regime, st$daynight) %in% paste(ex$regime, ex$daynight)] %in% "none"))
  # u* threshold share and sigma_w / u*
  i <- b & r$daynight == "night"
  g1 <- function(v, reg = "all", dn = "night") desc$value[desc$variable == v & desc$subset %in% "valid_qc_used" &
                                                            desc$regime == reg & desc$daynight == dn]
  expect_equal(g1("share_below_ustar_thr"), mean(r$ustar[i] < 0.2))
  expect_equal(g1("sw_ustar_median"), stats::median(0.4 / r$ustar[i]))
  expect_equal(g1("n"), sum(i))
})

test_that("descriptor table layout", {
  d <- withr::local_tempdir()
  s <- make_alpine_site(d)
  r <- alpine_row_vars(fp_site_rows(s, 2024, d, sets = "etc", subsets = c("valid", "valid_qc_used"),
                                    tab_cols = alpine_tab_cols()), rep(0.3, 400), rep(NA_real_, 400))
  cells <- dplyr::tibble(site_id = s, site_altitude_m = 1600, orog_m = 1960, dalt_m = -360, sdor_max_m = 380,
                         orog_min_m = 1200, orog_max_m = 2500)
  desc <- alpine_site_descriptors(r, list(valid_qc_used = alpine_basis(r, "etc", "valid_qc_used"),
                                          valid = alpine_basis(r, "etc", "valid")), "common_pe",
                                  list(cells = cells, igbp = "ENF", ustar_thr_source = "none (test)"))
  expect_named(desc, c("site_id", "block", "variable", "unit", "subset", "regime", "daynight", "value", "value_chr",
                       "row_basis"))
  expect_setequal(unique(desc$block), c("rows", "pblh", "turbulence", "flows", "footprint", "terrain", "tower"))
  site_rows <- desc[is.na(desc$subset), ]
  expect_true(all(is.na(site_rows$regime) & is.na(site_rows$daynight)))
  expect_equal(site_rows$value[site_rows$variable == "orog_range_9cells_m"], 1300)
  expect_equal(site_rows$value[site_rows$variable == "site_minus_orog_m"], -360)
  expect_equal(site_rows$value_chr[site_rows$variable == "igbp"], "ENF")
  expect_equal(site_rows$value[site_rows$variable == "zm_m"], 20)
  # 5 regimes x 3 day/night x 2 subsets, each with the same variables; no duplicated keys
  key <- paste(desc$variable, desc$subset, desc$regime, desc$daynight)
  expect_false(anyDuplicated(key) > 0)
  expect_setequal(unique(desc$regime[!is.na(desc$regime)]), c("all", "unstable", "neutral", "stable", "neutral_stable"))
  expect_true(all(is.na(desc$value[desc$variable == "share_below_ustar_thr"])))
  # flags_only: footprint values NA
  r2 <- r
  for (cs in c("pblhmod", "era5")) r2[[paste0("ta__", cs)]] <- NA_real_
  r2 <- alpine_row_vars(r2, rep(0.3, 400), rep(NA_real_, 400))
  d2 <- alpine_site_descriptors(r2, list(valid_qc_used = alpine_basis(r2, "etc", "valid_qc_used", FALSE)), "flags_only",
                                list(cells = NULL, igbp = NA, ustar_thr_source = "VUT"))
  expect_true(all(d2$value[d2$variable == "ta_diff_n"] == 0))
  expect_true(all(is.na(d2$value_chr[d2$variable == "spec15_status"])))
  expect_true(all(d2$row_basis == "flags_only"))
})

test_that("row-level R2: exact linear relation gives 1, strata and models listed", {
  set.seed(3)
  n <- 60
  x1 <- rnorm(n)
  x2 <- runif(n, 0.01, 0.5)
  nr <- dplyr::tibble(site_id = "XX-Tst", zl_class = rep(c("neutral", "stable"), length.out = n),
                      ta_diff_pp = 2 - 3 * x1, lnx80_era5_pblhmod = 0.5 * x1 + 4 * x2 + rnorm(n, 0, 0.01),
                      ln_h_era5_pblhmod = x1, zm_h_era5 = x2)
  rl <- alpine_row_level(nr)
  expect_equal(nrow(rl), 3 * 2 * 3)
  g <- rl[rl$stratum == "neutral_stable_night" & rl$response == "ta_diff_pp" & rl$predictors == "ln_h_era5_pblhmod", ]
  expect_equal(g$r2, 1, tolerance = 1e-9)
  expect_equal(g$coef_1, -3, tolerance = 1e-9)
  expect_equal(g$n, n)
  b <- rl[rl$stratum == "neutral_stable_night" & rl$response == "lnx80_era5_pblhmod" & rl$predictors == "both", ]
  expect_gt(b$r2, 0.99)
  expect_equal(c(b$coef_1, b$coef_2), c(0.5, 4), tolerance = 0.01)
  # too few rows -> NA
  rl2 <- alpine_row_level(nr[1:2, ])
  expect_true(all(is.na(rl2$r2)))
})

test_that("wind-direction histogram bins wrap at 0/360", {
  r <- dplyr::tibble(site_id = "XX-Tst", wind_dir = c(355, 5, 14, 16, 180, 340), daynight = "night")
  h <- alpine_wd_hist(r, rep(TRUE, 6))
  hn <- h[h$daynight == "night", ]
  expect_equal(hn$n[hn$dir_bin == 0], 3)   # 355, 5, 14 (bin 345-15)
  expect_equal(hn$n[hn$dir_bin == 30], 1)  # 16
  expect_equal(hn$n[hn$dir_bin == 330], 1) # 340
  expect_equal(hn$n[hn$dir_bin == 180], 1)
  expect_equal(sum(hn$share), 1)
  expect_true(all(is.na(h$share[h$daynight == "day"])))
})

test_that("plot helpers draw on a synthetic descriptor table", {
  sites <- c("AA-One", "BB-Two", "CZ-BK1")
  mk <- function(s, var, reg, dn, val, sub = "valid_qc_used") dplyr::tibble(
    site_id = s, block = "x", variable = var, unit = "-", subset = sub, regime = reg, daynight = dn, value = val,
    value_chr = NA_character_, row_basis = "common_pe")
  desc <- dplyr::bind_rows(lapply(seq_along(sites), function(k) dplyr::bind_rows(
    mk(sites[k], "ln_h_pblhmod_era5_median", "neutral_stable", "night", k * 0.5),
    mk(sites[k], "ta_diff_p10", "stable", "night", -k * 3), mk(sites[k], "ta_diff_p10", "neutral", "night", -k),
    mk(sites[k], "ta_diff_n", "stable", "night", c(400, 30, 100)[k]), mk(sites[k], "ta_diff_n", "neutral", "night", 200),
    mk(sites[k], "ec_height_m", NA, NA, 10 * k, sub = NA), mk(sites[k], "sdor_max_m", NA, NA, 100 * k, sub = NA))))
  ft <- alpine_drivers_table(desc)
  expect_equal(nrow(ft), 6)
  expect_equal(ft$below_min_n[ft$site_id == "BB-Two" & ft$panel == "stable night"], TRUE)
  expect_equal(ft$label[ft$site_id == "AA-One" & ft$panel == "stable night"], "AA-One")   # |p10| = 3 > 2
  expect_equal(ft$label[ft$site_id == "AA-One" & ft$panel == "neutral night"], "")        # |p10| = 1
  expect_equal(ft$label[ft$site_id == "CZ-BK1" & ft$panel == "neutral night"], "CZ-BK1")  # always labelled
  expect_s3_class(ggplot2::ggplot_build(plot_alpine_drivers(ft)), "ggplot_built")
  nr <- dplyr::tibble(site_id = rep(sites, each = 5), zl_class = "stable", ta_diff_pp = rnorm(15),
                      lnx80_era5_pblhmod = rnorm(15), ln_h_era5_pblhmod = rnorm(15), zm_h_era5 = runif(15, 0.01, 1))
  expect_s3_class(ggplot2::ggplot_build(plot_alpine_row_level(nr)), "ggplot_built")
  wh <- dplyr::bind_rows(lapply(sites, function(s) alpine_wd_hist(
    dplyr::tibble(site_id = s, wind_dir = runif(50, 0, 360), daynight = rep(c("day", "night"), 25)), rep(TRUE, 50))))
  expect_s3_class(ggplot2::ggplot_build(plot_alpine_flows(wh)), "ggplot_built")
})
