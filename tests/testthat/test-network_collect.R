# Spec 09 A4/A5: reference-year rule, year counts, gamma network summary, site flags.

test_that("footprint_year_counts: usable = both main cases valid under etc in valid_qc_used", {
  t <- as.POSIXct("2021-12-31 23:00", tz = "UTC") + 1800 * (0:5)
  tab <- dplyr::tibble(timestamp_start_utc = t, h_measured = c(NA, 500, NA, 600, NA, NA))
  sub <- dplyr::tibble(valid_qc_used__etc__pblhmod = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
                       valid_qc_used__etc__era5 = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
  y <- footprint_year_counts(tab, sub)
  expect_equal(y$year, c(2021L, 2022L))
  expect_equal(y$n_usable, c(2L, 2L))
  expect_equal(y$n_months, c(1L, 1L))
  expect_equal(y$n_measured, c(1L, 1L))
})

test_that("propose_reference_year: >= 9 months, most usable, 2 % tolerance, measured overlap, then most recent", {
  yc <- dplyr::tibble(year = 2020:2024, n_usable = c(10000, 9900, 9850, 9000, 12000),
                      n_months = c(12, 12, 12, 12, 8), n_measured = c(0, 0, 0, 0, 0))
  p <- propose_reference_year(yc)
  expect_equal(p$year, 2022L) # 2024 has 8 months; 2020-2022 within 2 % of 10000: most recent
  expect_equal(p$n_candidates, 3L)
  expect_match(p$reason, "most recent")
  yc$n_measured <- c(50, 0, 0, 0, 0)
  expect_equal(propose_reference_year(yc)$year, 2020L) # measured overlap preferred over the most recent
  expect_match(propose_reference_year(yc)$reason, "measured")
  yc$n_usable[2] <- 20000
  expect_equal(propose_reference_year(yc)$year, 2021L) # a clear maximum wins over measured overlap
  expect_true(is.na(propose_reference_year(dplyr::mutate(yc, n_months = 3))$year))
})

test_that("network_gamma_summary and network_flags", {
  m <- dplyr::bind_rows(lapply(c("A", "B", "C"), function(s) {
    dplyr::tibble(site_id = s, month = 1:12, median = c(A = 0.004, B = 0.005, C = 0.006)[[s]])
  }))
  su <- dplyr::tibble(site_id = c("A", "B", "C"), layer = "default", season = "all", share_2_levels = c(1, 0.3, 0))
  g <- network_gamma_summary(m, su)
  expect_equal(nrow(g), 13)
  expect_equal(g$median[g$month == "1"], 0.005)
  expect_equal(g$n_sites[g$month == "all"], 3L)
  expect_equal(g$n_sites_2_levels[1], 1L)

  st <- dplyr::tibble(site_id = c(rep("A", 5), "B", "B"), step = c("01", "11a", "02", "04", "08a", "01", "11a"),
                      state = c(rep("ok", 6), "failed"))
  ca <- dplyr::tibble(site_id = c("A", "A"), case = c("default", "max_B2"), gt5km = c(2, 40))
  cv <- dplyr::tibble(site_id = "A", case = "pblhmod", set = "etc", stratum = "all", n_rows = 100, n_valid_qc_used = 10)
  ce <- dplyr::tibble(site_id = c("A", "B"), cell_dist_km = c(5, 10), dalt_m = c(10, -362))
  f <- network_flags(c("A", "B"), st, ca, cv, ce)
  expect_equal(f$flag[f$site_id == "A"], c("pblh_gt_5km", "low_coverage"))
  expect_equal(f$flag[f$site_id == "B"], c("step_failed", "era5_orography"))
  expect_match(f$detail[f$flag == "step_failed"], "11a")
  expect_match(f$detail[f$flag == "era5_orography"], "\\+362 m")
})
