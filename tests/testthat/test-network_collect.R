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

test_that("propose_reference_year: >= 9 months, measured overlap within 80 %, most usable (2 % tolerance), most recent", {
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
  expect_equal(propose_reference_year(yc)$year, 2021L) # measured year below 80 % of the best: standard rule
  expect_equal(propose_reference_year(yc)$reason, "most usable half-hours") # one candidate: no tolerance wording
  expect_true(is.na(propose_reference_year(dplyr::mutate(yc, n_months = 3))$year))

  # spec 09 amendment 3: FI-Hyy year counts of the cluster (Cloudnet 2022-01 to 2024-04): 2021 has the most usable
  # half-hours, 2022 has 88 % of them and the most measured half-hours -> 2022
  fh <- dplyr::tibble(year = 2017:2026, n_usable = c(0, 0, 0, 3648, 3986, 3510, 3574, 3450, 3574, 2772),
                      n_months = c(0, 0, 0, 12, 12, 12, 12, 12, 12, 8), n_measured = c(0, 0, 0, 0, 0, 1200, 657, 231, 0, 0))
  p <- propose_reference_year(fh)
  expect_equal(p$year, 2022L)
  expect_equal(p$n_measured, 1200L)
  expect_equal(p$n_candidates, 3L) # 2022, 2023, 2024 all >= 80 % of 3986
  expect_match(p$reason, "most measured half-hours [(]1200[)], usable 88 % of the best year [(]>= 80 %[)]")
  # threshold above 2022's share: 2023 (90 %) is the only measured candidate left
  expect_equal(propose_reference_year(fh, measured_min_share = 0.89)$year, 2023L)
  # without measured overlap: unchanged standard rule (2021)
  expect_equal(propose_reference_year(dplyr::mutate(fh, n_measured = 0))$year, 2021L)
  expect_equal(propose_reference_year(fh, measured_min_share = 0.95)$year, 2021L)
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

  # NEE QC (spec 09 amendment 2): whole-period row; CUT sites flagged, empty valid_qc_used flagged, VUT ok not flagged
  nq <- dplyr::tibble(site_id = c("A", "A", "B", "B", "C"), year = c(NA, 2022L, NA, 2023L, NA),
                      nee_qc_source = c("VUT", "VUT", "CUT", "CUT", "none"), n_fluxnet = c(10, 0, 64272, 64272, 5),
                      n_qc = c(10, 0, 50197, 50197, 0), n_qc0 = c(6, 0, 11524, 11524, 0), n_used = c(5, 0, 9000, 9000, 0),
                      status = c("ok", "no_fluxnet", "ok", "ok", "all_missing"))
  st3 <- dplyr::tibble(site_id = rep(c("A", "B", "C"), each = 5), step = rep(c("01", "11a", "02", "04", "08a"), 3), state = "ok")
  empty_ca <- dplyr::tibble(site_id = character(), case = character(), gt5km = integer())
  empty_cv <- dplyr::tibble(site_id = character(), case = character(), set = character(), stratum = character(),
                            n_rows = integer(), n_valid_qc_used = integer())
  empty_ce <- dplyr::tibble(site_id = character(), cell_dist_km = numeric(), dalt_m = numeric())
  f <- network_flags(c("A", "B", "C"), st3, empty_ca, empty_cv, empty_ce, nee_qc = nq)
  expect_equal(f$site_id, c("B", "C"))
  expect_equal(f$flag, c("nee_qc_cut", "nee_qc"))
  expect_match(f$detail[1], "NEE_CUT_REF_QC \\(64272 FLUXNET rows, 11524 QC 0")
  expect_match(f$detail[2], "NEE_VUT/CUT_REF_QC all_missing")
  # old 08_nee_qc tables (no year, no source column): VUT assumed
  old <- dplyr::tibble(site_id = "A", n_fluxnet = 5, n_qc = 5, n_qc0 = 0, n_used = 0, status = "none_measured")
  expect_match(network_flags("A", st3, empty_ca, empty_cv, empty_ce, nee_qc = old)$detail, "NEE_VUT_REF_QC none_measured")
})

test_that("stations not run: read, listed as not_run (not as a missing step); config lists consistent", {
  p <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("site_id,reason,detail,decided", "XX-Nor,no L2 data,no data in the L2I 2026 release,G 2026-11-03"), p)
  nr <- read_sites_not_run(p)
  expect_equal(nr$site_id, "XX-Nor")
  expect_equal(nrow(read_sites_not_run(file.path(tempdir(), "none.csv"))), 0)
  bad <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("site_id,reason,detail,decided", "XX-Nor,,x,y"), bad)
  expect_error(read_sites_not_run(bad), "reason")
  st <- dplyr::tibble(site_id = rep("A", 5), step = c("01", "11a", "02", "04", "08a"), state = "ok")
  empty_cv <- dplyr::tibble(site_id = character(), case = character(), set = character(), stratum = character(),
                            n_rows = integer(), n_valid_qc_used = integer())
  empty_ce <- dplyr::tibble(site_id = character(), cell_dist_km = numeric(), dalt_m = numeric())
  empty_ca <- dplyr::tibble(site_id = character(), case = character(), gt5km = integer())
  # also in the site list and without any run status: still "not_run", never "step_missing"
  f <- network_flags(c("A", "XX-Nor"), st, empty_ca, empty_cv, empty_ce, not_run = nr)
  expect_equal(f$site_id, "XX-Nor")
  expect_equal(f$flag, "not_run")
  expect_equal(f$detail, "not run (no L2 data): no data in the L2I 2026 release [G 2026-11-03]")
  expect_equal(network_flags(c("A", "XX-Nor"), st, empty_ca, empty_cv, empty_ce)$flag, "step_missing")

  # the committed lists: run_sites.txt + sites_not_run.csv = sites.csv, disjoint (IT-SR2 not run; 39 sites)
  ids <- trimws(readLines(here::here("config", "run_sites.txt"), warn = FALSE))
  ids <- ids[nzchar(ids) & !startsWith(ids, "#")]
  nr_cfg <- read_sites_not_run(here::here("config", "sites_not_run.csv"))
  all_sites <- read_sites(here::here("config", "sites.csv"))$site_id
  expect_length(intersect(ids, nr_cfg$site_id), 0)
  expect_setequal(c(ids, nr_cfg$site_id), all_sites)
  expect_true("IT-SR2" %in% nr_cfg$site_id && "IT-SR2" %in% all_sites)
  expect_length(ids, 39)
})

test_that("reference-year exception 'coverage' (spec 09 amendment 2): only when no year qualifies; config list", {
  yc <- dplyr::tibble(year = 2022:2025, n_usable = 0L, n_months = 0L, n_measured = 0L,
                      n_valid_qc = c(1200L, 4100L, 4100L, 3000L), n_months_valid_qc = c(4L, 6L, 5L, 5L))
  std <- propose_reference_year(yc)
  expect_true(is.na(std$year))
  expect_equal(std$rule, "standard")
  ex <- propose_reference_year(yc, rule = "coverage", rule_reason = "remote site")
  expect_equal(ex$year, 2024L) # tie on n_valid_qc between 2023 and 2024: most recent
  expect_equal(ex$rule, "coverage")
  expect_match(ex$reason, "no year with >= 9 months.*most valid_qc half-hours [(]4100, 5 months[)].*remote site")
  # a year qualifying under the standard rule wins: the exception is not used
  yc2 <- dplyr::mutate(yc, n_usable = c(0L, 0L, 0L, 900L), n_months = c(0L, 0L, 0L, 10L))
  expect_equal(propose_reference_year(yc2, rule = "coverage")$year, 2025L)
  expect_equal(propose_reference_year(yc2, rule = "coverage")$rule, "standard")
  # old year counts without n_valid_qc: NA with the reason
  expect_match(propose_reference_year(yc[1:4], rule = "coverage")$reason, "re-run step 08a")
  expect_error(propose_reference_year(yc, rule = "other"))

  # footprint_year_counts gives n_valid_qc (both main cases, etc, valid_qc)
  t <- as.POSIXct("2023-06-01", tz = "UTC") + 1800 * (0:3)
  sub <- dplyr::tibble(valid_qc_used__etc__pblhmod = FALSE, valid_qc_used__etc__era5 = FALSE,
                       valid_qc__etc__pblhmod = c(TRUE, TRUE, FALSE, TRUE), valid_qc__etc__era5 = c(TRUE, FALSE, TRUE, TRUE))
  y <- footprint_year_counts(dplyr::tibble(timestamp_start_utc = t), sub)
  expect_equal(y$n_usable, 0L)
  expect_equal(y$n_valid_qc, 2L)
  expect_equal(y$n_months_valid_qc, 1L)

  # the committed exception list: GL-ZaF, rule coverage, a run site
  ex_cfg <- read_reference_year_exceptions(here::here("config", "reference_year_exceptions.csv"))
  expect_equal(ex_cfg$site_id, "GL-ZaF")
  expect_equal(ex_cfg$rule, "coverage")
  ids <- trimws(readLines(here::here("config", "run_sites.txt"), warn = FALSE))
  expect_true(all(ex_cfg$site_id %in% ids))
  expect_equal(nrow(read_reference_year_exceptions(file.path(tempdir(), "none.csv"))), 0)
  bad <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("site_id,rule,reason,decided", "XX-Aaa,months6,x,y"), bad)
  expect_error(read_reference_year_exceptions(bad), "months6")
})

test_that("network figure (spec 09 amendment 2): IQR bar for PBLHmod only, PBLHmod-clim open circles on top", {
  met <- dplyr::bind_rows(lapply(c("AA-One", "BB-Two"), function(s) {
    dplyr::bind_rows(lapply(c("1_pblhmod_default_vs_era5_center", "6_pblhmod_clim_vs_era5_center"), function(p) {
      dplyr::tibble(site_id = s, pair = p, stratum = rep(c("all", "regime=stable", "regime=convective", "regime=mixed"), 4),
                    metric = rep(c("median_log_ratio", "p25_log_ratio", "p75_log_ratio", "n_log"), each = 4),
                    value = rep(c(0.1, -0.3, 0.5, 100), each = 4) + if (startsWith(p, "6_")) 0.05 else 0)
    }))
  }))
  d <- network_logratio_table(met)
  expect_setequal(unique(d$regime), c("all", "stable", "convective")) # mixed dropped
  expect_equal(d$p25[d$site_id == "AA-One" & d$source == "pblhmod_default" & d$regime == "all"], -0.3)
  p <- plot_network_logratio(d)
  b <- ggplot2::ggplot_build(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], "")
  i_bar <- which(geoms == "GeomLinerange")
  i_pts <- which(geoms == "GeomPoint")
  expect_length(i_bar, 1)
  expect_length(i_pts, 3) # PBLHmod, side-result constants (empty here), PBLHmod-clim
  # bar from p25 to p75, PBLHmod only
  expect_equal(unique(p$layers[[i_bar]]$data$source), "pblhmod_default")
  expect_equal(range(b$data[[i_bar]]$xmin), c(-0.3, -0.3))
  expect_equal(range(b$data[[i_bar]]$xmax), c(0.5, 0.5))
  # clim drawn last (on top), open circles; default filled dots
  expect_equal(unique(p$layers[[max(i_pts)]]$data$source), "pblhmod_clim")
  expect_equal(max(i_pts), length(p$layers))
  expect_equal(unique(b$data[[max(i_pts)]]$shape), 1)
  expect_equal(unique(b$data[[min(i_pts)]]$shape), 16)
  expect_match(p$labels$caption, "Stable: independent of γ")
  # metrics older than amendment 2 (no p25/p75): no bar rows, figure still builds
  d0 <- network_logratio_table(met[!met$metric %in% c("p25_log_ratio", "p75_log_ratio"), ])
  expect_true(all(is.na(d0$p25)))
  p0 <- plot_network_logratio(d0)
  expect_equal(nrow(p0$layers[[which(vapply(p0$layers, function(l) class(l$geom)[1], "") == "GeomLinerange")]]$data), 0)
  expect_s3_class(ggplot2::ggplot_build(p0), "ggplot_built")
})

test_that("supplementary network figure (spec 09 amendment 3): gamma_med and gamma_site, clim still on top", {
  prs <- network_logratio_pairs(supp = TRUE)
  expect_equal(names(prs), c("pblhmod_default", "pblhmod_clim", "pblhmod_gamma_med", "pblhmod_gamma_site"))
  expect_equal(network_logratio_pairs(), prs[1:2])
  met <- dplyr::bind_rows(lapply(c("AA-One", "BB-Two"), function(s) {
    dplyr::bind_rows(lapply(seq_along(prs), function(j) {
      dplyr::tibble(site_id = s, pair = unname(prs[j]), stratum = rep(c("all", "regime=stable"), 4),
                    metric = rep(c("median_log_ratio", "p25_log_ratio", "p75_log_ratio", "n_log"), each = 2),
                    value = rep(c(0.1 * j, -0.3, 0.5, 100), each = 2))
    }))
  }))
  expect_setequal(unique(network_logratio_table(met)$source), c("pblhmod_default", "pblhmod_clim")) # main figure
  d <- network_logratio_table(met, pairs = prs)
  expect_equal(d$median[d$site_id == "AA-One" & d$source == "pblhmod_gamma_site" & d$regime == "all"], 0.4)
  p <- plot_network_logratio(d)
  b <- ggplot2::ggplot_build(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], "")
  i_pts <- which(geoms == "GeomPoint")
  expect_setequal(unique(p$layers[[i_pts[2]]]$data$source), c("pblhmod_gamma_med", "pblhmod_gamma_site"))
  expect_setequal(unique(b$data[[i_pts[2]]]$shape), c(17, 2))
  expect_equal(unique(p$layers[[max(i_pts)]]$data$source), "pblhmod_clim") # clim still last
  expect_equal(unique(p$layers[[which(geoms == "GeomLinerange")]]$data$source), "pblhmod_default") # bar: PBLHmod only
  expect_true(all(c("pblhmod_gamma_med", "pblhmod_gamma_site") %in% style_ids("source")))
})

test_that("network and site gamma constants (spec 09 amendment 3)", {
  a <- dplyr::tibble(gamma = c(0.004, 0.005, NA, 0.006))
  b <- dplyr::tibble(gamma = c(0.003, 0.003))
  cc <- dplyr::tibble(gamma = c(NA_real_, NA_real_))
  expect_equal(gamma_site_constant(a), dplyr::tibble(gamma = 0.005, n_days = 3L))
  expect_true(is.na(gamma_site_constant(cc)$gamma))
  g <- gamma_network_median(list(A = a, B = b, C = cc))
  expect_equal(g$gamma_med, 0.004) # median of the site medians 0.005 (A), 0.003 (B); C has no valid day
  expect_equal(g$gamma_med_pooled, 0.004) # reference: pooled 0.003 0.003 0.004 0.005 0.006
  # the two differ when the sites hold different numbers of valid days: A weighs 3 days, B 1 day in the pooled value
  g2 <- gamma_network_median(list(A = a, B = dplyr::tibble(gamma = 0.001)))
  expect_equal(c(g2$gamma_med, g2$gamma_med_pooled), c(0.003, 0.0045))
  expect_equal(c(g$n_sites, g$n_days), c(2L, 5L))
  expect_equal(g$sites, "A;B")
  expect_error(gamma_network_median(list()), "non-empty")
  expect_error(gamma_network_median(list(C = cc)), "No valid")
  # written only if changed (step-02 inputs hash uses size and mtime); read back
  cfg <- list(derived_dir = withr::local_tempdir())
  p <- gamma_network_path(cfg)
  expect_null(read_gamma_network(p))
  expect_true(write_gamma_network(g, p))
  expect_false(write_gamma_network(g, p))
  expect_equal(read_gamma_network(p)$gamma_med, 0.004)
  expect_true(write_gamma_network(dplyr::mutate(g, gamma_med = 0.0045), p))
  bad <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("gamma_med", "-1"), bad)
  expect_error(read_gamma_network(bad), "positive")
  # step 02 hashes the file: input group gamma_network
  expect_true("gamma_network" %in% pipeline_steps()[["02"]]$inputs)
  expect_equal(site_input_files("XX-Tst", cfg, "gamma_network")$path, p)
  expect_equal(nrow(site_input_files("XX-Tst", list(derived_dir = withr::local_tempdir()), "gamma_network")), 0)
})
