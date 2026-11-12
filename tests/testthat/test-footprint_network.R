# Spec 15: network footprint summary and gamma supplementary tables (synthetic rows only).

cases3 <- c("pblhmod", "pblhmod_clim", "era5")

# Synthetic rows in the layout of fp_site_rows(): every row valid in every set / subset unless changed.
make_rows <- function(site_id, ta, zl = "stable", dn = "night", h = 500, x80 = 300, regime = "stable",
                      season = "JJA") {
  n <- length(ta$pblhmod)
  r <- dplyr::tibble(site_id = site_id, timestamp_start_utc = as.POSIXct("2025-06-01", tz = "UTC") + 1800 * seq_len(n),
                     zl_class = rep_len(zl, n), daynight = rep_len(dn, n), season = rep_len(season, n),
                     regime = rep_len(regime, n))
  for (cs in cases3) {
    r[[paste0("h__", cs)]] <- if (is.list(h)) h[[cs]] else h
    r[[paste0("ta__", cs)]] <- ta[[cs]]
    r[[paste0("x80__", cs)]] <- if (is.list(x80)) x80[[cs]] else x80
    r[[paste0("ok__", cs)]] <- TRUE
    for (s in c("etc", "kljun_norsl")) for (u in c("valid", "valid_qc", "valid_qc_used")) {
      r[[paste("v", s, u, cs, sep = "__")]] <- TRUE
    }
  }
  r
}

test_that("pairs: difference = first minus second, in percentage points", {
  p <- fp_pairs()
  expect_equal(p$case_a, c("era5", "pblhmod_clim", "era5"))
  expect_equal(p$case_b, c("pblhmod", "pblhmod", "pblhmod_clim"))
  r <- make_rows("XX-Aaa", list(pblhmod = rep(0.60, 20), pblhmod_clim = rep(0.65, 20), era5 = rep(0.50, 20)))
  ta <- fp_target_area_table(r, "etc", "valid_qc_used")
  a <- ta[ta$scope == "site" & ta$regime == "all" & ta$daynight == "all", ]
  expect_equal(a$median_pp[a$pair == "era5_vs_pblhmod"], -10)
  expect_equal(a$median_pp[a$pair == "pblhmod_clim_vs_pblhmod"], 5)
  expect_equal(a$median_pp[a$pair == "era5_vs_pblhmod_clim"], -15)
  expect_equal(a$median_ta_a_pct[a$pair == "era5_vs_pblhmod"], 50)
  # strata: the crossed cell and the margins
  expect_setequal(paste(ta$regime, ta$daynight)[ta$scope == "site" & ta$pair == "era5_vs_pblhmod"],
                  c("stable night", "all night", "stable all", "all all"))
  f <- fp_fetch_table(make_rows("XX-Aaa", list(pblhmod = rep(0.6, 5), pblhmod_clim = rep(0.6, 5), era5 = rep(0.6, 5)),
                                x80 = list(pblhmod = 100, pblhmod_clim = 200, era5 = 400)), "etc", "valid")
  f1 <- f[f$scope == "site" & f$regime == "all" & f$daynight == "all", ]
  expect_equal(f1$median_lnratio[f1$pair == "era5_vs_pblhmod"], log(4))
  expect_equal(f1$median_lnratio[f1$pair == "era5_vs_pblhmod_clim"], log(2))
  expect_equal(f1$median_x80_a_m[f1$pair == "era5_vs_pblhmod"], 400)
})

test_that("common_pe: a row missing in one case is left out of every pair; network scopes", {
  r1 <- make_rows("XX-Aaa", list(pblhmod = rep(0.6, 10), pblhmod_clim = rep(0.6, 10), era5 = c(rep(0.5, 9), 0.0)))
  r1$v__etc__valid__pblhmod_clim[10] <- FALSE # the row with the big ERA5 difference is not valid in clim
  r2 <- make_rows("XX-Bbb", list(pblhmod = rep(0.6, 30), pblhmod_clim = rep(0.6, 30), era5 = rep(0.4, 30)))
  r <- dplyr::bind_rows(r1, r2)
  expect_equal(sum(fp_common(r, "etc", "valid")), 39)
  ta <- fp_target_area_table(r, "etc", "valid")
  e <- ta[ta$pair == "era5_vs_pblhmod" & ta$regime == "all" & ta$daynight == "all", ]
  expect_equal(e$p10_pp[e$scope == "site" & e$site_id == "XX-Aaa"], -10)
  expect_equal(e$n[e$scope == "network_pooled"], 39)
  expect_equal(e$n_sites[e$scope == "network_pooled"], 2)
  expect_equal(e$median_pp[e$scope == "network_pooled"], -20)  # 30 of 39 rows at -20
  expect_equal(e$median_pp[e$scope == "network_site_median"], -15) # median of -10 and -20
})

test_that("top 10 %: the rows with the largest |TA difference|, shares and removals", {
  n <- 100
  ta_e <- 0.6 - (1:n) / 1000 # |difference| 0.1 ... 10 pp
  r <- make_rows("XX-Aaa", list(pblhmod = rep(0.6, n), pblhmod_clim = rep(0.6, n), era5 = ta_e),
                 zl = c(rep("neutral", 90), rep("stable", 10)), dn = c(rep("day", 95), rep("night", 5)))
  r$v__etc__valid_qc__pblhmod[c(91:96)] <- FALSE      # 6 of the top 10 fail the NEE QC
  r$v__etc__valid_qc_used__era5[c(97, 98)] <- FALSE   # 2 more fail "used" only
  r$v__etc__valid_qc_used__pblhmod[c(91:96)] <- FALSE # (QC failures also fail "used")
  sel <- fp_top_rows(abs(100 * (ta_e - 0.6)), rep("s", n))
  expect_equal(which(sel$top), 91:100)
  t <- fp_topaffected_table(r, "etc", fp_pairs()[1, ])
  s <- t[t$scope == "site", ]
  expect_equal(s$n_top, 10)
  expect_equal(s$share_stable, 1)
  expect_equal(s$share_night, 0.5)
  expect_equal(s$share_removed_by_qc, 0.6)
  expect_equal(s$share_removed_by_qc_or_used, 0.8)
  expect_equal(s$threshold_pp, unname(stats::quantile((1:n) / 10, 0.9)))
  expect_equal(t$n_top[t$scope == "network_pooled"], 10)
})

test_that("ICOS criteria: Rebmann 90 % for >= 50 % of the time, labelling rule 70 % for >= 70 % of the good data", {
  expect_true(icos_ta_criteria(c(rep(0.95, 5), rep(0.5, 5)))$pass_ta90_50pct)  # exactly 50 %: pass
  expect_false(icos_ta_criteria(c(rep(0.95, 4), rep(0.5, 6)))$pass_ta90_50pct)
  expect_true(icos_ta_criteria(c(rep(0.70, 7), rep(0.1, 3)))$pass_lab_70_70)  # exactly 70 % at exactly 0.70: pass
  expect_false(icos_ta_criteria(c(rep(0.69, 7), rep(0.9, 3)))$pass_lab_70_70)
  expect_true(is.na(icos_ta_criteria(numeric())$pass_ta90_50pct))
  # a site: all rows pass (60 % >= 0.9), stable rows fail (40 %); ERA5 differs in the labelling rule
  zl <- c(rep("neutral", 5), rep("stable", 5))
  ta <- list(pblhmod = c(rep(0.95, 4), 0.75, rep(0.95, 2), rep(0.75, 3)), pblhmod_clim = c(rep(0.95, 4), 0.75, rep(0.95, 2), rep(0.75, 3)),
             era5 = c(rep(0.95, 4), 0.5, rep(0.95, 2), rep(0.5, 3)))
  r <- make_rows("XX-Aaa", ta, zl = zl)
  x <- fp_icos_criteria_table(r, "etc", "valid_qc_used")
  expect_equal(x$share_ta90, rep(0.6, 3))
  expect_true(all(x$pass_ta90))
  expect_equal(x$share_ta90_stable, rep(0.4, 3))
  expect_false(any(x$pass_ta90_stable))
  expect_equal(x$icos_lab_share, c(1, 1, 0.6))
  expect_equal(x$icos_lab_pass, c(TRUE, TRUE, FALSE))
  expect_true(all(x$icos_lab_differs))
  expect_false(any(x$ta90_differs))
})

test_that("exceptions: |p10| or |p90| strictly above the threshold, in crossed strata only", {
  r1 <- make_rows("XX-Aaa", list(pblhmod = rep(0.60, 20), pblhmod_clim = rep(0.60, 20), era5 = rep(0.55, 20)),
                  h = list(pblhmod = 300, pblhmod_clim = 300, era5 = 900))      # exactly -5 pp
  r2 <- make_rows("XX-Bbb", list(pblhmod = rep(0.60, 20), pblhmod_clim = rep(0.60, 20), era5 = rep(0.5499, 20)),
                  h = list(pblhmod = 300, pblhmod_clim = 300, era5 = 900))      # -5.01 pp
  r <- dplyr::bind_rows(r1, r2)
  ta <- fp_target_area_table(r, "etc", "valid_qc_used")
  ex <- fp_exceptions_table(ta, r, 5)
  expect_equal(unique(ex$site_id), "XX-Bbb")
  e <- ex[ex$pair == "era5_vs_pblhmod", ]
  expect_equal(nrow(e), 1) # crossed stratum only (stable x night), not the margins
  expect_equal(c(e$regime, e$daynight), c("stable", "night"))
  expect_equal(e$smaller_case, "era5")
  expect_equal(c(e$median_h_a_m, e$median_h_b_m), c(900, 300))
  expect_equal(nrow(fp_exceptions_table(ta, r, 10)), 0)
})

test_that("fp_site_rows reads the input table and the hh chunks of the reference year", {
  d <- withr::local_tempdir()
  sd <- file.path(d, "XX-Tst", "chunks")
  dir.create(sd, recursive = TRUE)
  t <- as.POSIXct("2025-01-31 23:00", tz = "UTC") + 1800 * 0:5 # spans January and February
  tab <- dplyr::tibble(timestamp_start_utc = t, zm_m = 20, L = c(-100, 1000, 50, NA, 10, -10),
                       day = rep(c(TRUE, FALSE), 3), regime = "stable")
  sub <- dplyr::tibble(valid__etc__pblhmod = TRUE, valid__etc__era5 = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))
  saveRDS(list(tab = tab, subsets = sub), file.path(d, "XX-Tst", "XX-Tst_footprint_inputs.rds"))
  for (cs in c("pblhmod", "era5")) for (mo in c("2025-01", "2025-02")) {
    i <- format(t, "%Y-%m") == mo
    saveRDS(list(), file.path(sd, paste0("XX-Tst_", cs, "_", mo, "_meta.rds")))
    readr::write_csv(dplyr::tibble(timestamp_start_utc = t[i], h = 500, x_80 = 300, ffp_ok = TRUE,
                                   f_target_area = if (cs == "era5") 0.5 else 0.6),
                     file.path(sd, paste0("XX-Tst_", cs, "_", mo, "_hh.csv.gz")))
  }
  r <- fp_site_rows("XX-Tst", 2025, d, cases = c("pblhmod", "pblhmod_clim", "era5"), sets = "etc", subsets = "valid")
  expect_equal(nrow(r), 6)
  expect_equal(r$zl_class, c("unstable", "neutral", "stable", NA, "stable", "unstable")) # zm/L with zm = 20 m
  expect_equal(r$daynight, rep(c("day", "night"), 3))
  expect_equal(r$ta__era5, rep(0.5, 6))
  expect_true(all(is.na(r$ta__pblhmod_clim)))
  expect_false(any(r$v__etc__valid__pblhmod_clim))
  expect_equal(r$v__etc__valid__era5, c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_setequal(attr(r, "missing"), c("pblhmod_clim 2025-01", "pblhmod_clim 2025-02"))
  expect_false(any(fp_common(r, "etc", "valid")))                       # clim missing: no common_pe row
  expect_equal(sum(fp_common(r, "etc", "valid", c("pblhmod", "era5"))), 5)
})

test_that("write_csv_header writes a commented header readable with comment = '#'", {
  p <- withr::local_tempfile(fileext = ".csv")
  write_csv_header(dplyr::tibble(a = 1:2, b = c("x", "y")), p, c("definitions", "units"))
  expect_equal(readLines(p)[1:2], c("# definitions", "# units"))
  expect_equal(readr::read_csv(p, comment = "#", show_col_types = FALSE)$a, 1:2)
})

test_that("gamma supplementary tables: site quantiles and a reproducible site bootstrap", {
  mk <- function(g, start = "2024-01-01") dplyr::tibble(date = as.Date(start) + seq_along(g) - 1, gamma = g)
  daily <- list(`XX-Aaa` = mk(c(0.004, 0.005, 0.006, NA)), `XX-Bbb` = mk(c(0.002, 0.003)),
                `XX-Ccc` = mk(seq(0.006, 0.008, length.out = 400)))
  st <- gamma_site_table(daily, c(`XX-Aaa` = "ENF", `XX-Bbb` = "GRA", `XX-Ccc` = "CRO"))
  expect_equal(st$gamma_site, c(0.005, 0.0025, 0.007))
  expect_equal(st$n_days, c(3L, 2L, 400L))
  expect_equal(st$months_used[st$site_id == "XX-Ccc"], "1-12")
  expect_equal(st$months_used[st$site_id == "XX-Aaa"], "1")
  expect_equal(st$igbp, c("ENF", "GRA", "CRO"))
  n1 <- gamma_network_table(daily, n_boot = 500, seed = 1)
  n2 <- gamma_network_table(daily, n_boot = 500, seed = 1)
  expect_equal(n1, n2)
  expect_equal(n1$gamma_med, 0.005)
  expect_equal(c(n1$site_median_min, n1$site_median_max), c(0.0025, 0.007))
  expect_true(n1$gamma_med_lo95 <= n1$gamma_med && n1$gamma_med <= n1$gamma_med_hi95)
  expect_equal(n1$n_days, 405L)
})
