# Spec 11 amendment 1: reference per site, empty stratum sets, decision and season tables of the gamma evaluation.

ms_row <- function(...) {
  utils::modifyList(list(icos_site = "XX-Aer", provider = "MeteoSwiss", provider_site = "AAA", distance_km = 3,
                         method_class = "aerosol", primary_source_id = "meteoswiss_cabam_mlh", status = "used",
                         borderline = FALSE), list(...))
}
ms_tab <- function(...) dplyr::bind_rows(lapply(list(...), dplyr::as_tibble))

test_that("measured_reference() takes the primary source, its sets and the pairing flags", {
  ms <- ms_tab(ms_row(),
               ms_row(icos_site = "XX-Aer", provider = "cloudnet", provider_site = "aaa", primary_source_id = NA,
                      status = "inventory", method_class = "turbulence"),
               ms_row(icos_site = "XX-Trb", provider = "cloudnet", provider_site = "bbb", distance_km = 0.5,
                      method_class = "turbulence", primary_source_id = "cloudnet_eps_mlh"),
               ms_row(icos_site = "XX-Far", provider_site = "CCC", distance_km = 24.9, borderline = TRUE),
               ms_row(icos_site = "XX-Tfr", provider = "cloudnet", provider_site = "ddd", distance_km = 10,
                      method_class = "turbulence", primary_source_id = "cloudnet_eps_mlh"))
  a <- measured_reference(ms, "XX-Aer")
  expect_equal(a$src, "meteoswiss_cabam_mlh")
  expect_equal(a$type, "aerosol")
  expect_equal(a$sets, c("window_convective", "convective"))
  expect_equal(a$primary_set, "window_convective")
  expect_equal(nrow(a$row), 1)
  t <- measured_reference(ms, "XX-Trb")
  expect_equal(t$sets, c("convective", "all_regimes"))
  expect_equal(t$distance_class, "all_regimes")
  # turbulence beyond 5 km: the distance rule keeps the convective set only
  expect_equal(measured_reference(ms, "XX-Tfr")$sets, "convective")
  f <- measured_reference(ms, "XX-Far")
  expect_true(f$borderline)
  expect_equal(f$distance_class, "daytime_convective")
  expect_null(measured_reference(ms, "XX-Non"))
})

test_that("measured_reference() refuses what the comparison cannot use", {
  expect_error(measured_reference(ms_tab(ms_row(), ms_row(provider_site = "BBB")), "XX-Aer"), "primary")
  expect_error(measured_reference(ms_tab(ms_row(method_class = "turbulence")), "XX-Aer"), "Cloudnet")
  expect_error(measured_reference(ms_tab(ms_row(method_class = "radiosonde")), "XX-Aer"), "not supported")
  expect_error(measured_reference(ms_tab(ms_row(distance_km = 60)), "XX-Aer"), "not used")
})

test_that("measured_reference() on config/measured_sites.csv gives the four-site references of spec 11 D", {
  ms <- read_measured_sites(here::here("config", "measured_sites.csv"))
  old <- list(`FR-Gri` = c("sirta_mlh_chm15k", "aerosol"), `FI-Hyy` = c("cloudnet_eps_mlh", "turbulence"),
              `CH-Dav` = c("meteoswiss_stratfinder_mlh", "aerosol"))
  for (s in names(old)) {
    r <- measured_reference(ms, s)
    expect_equal(c(r$src, r$type), old[[s]], info = s)
  }
  expect_equal(measured_reference(ms, "FI-Hyy")$sets, c("convective", "all_regimes"))
  expect_equal(measured_reference(ms, "FR-Gri")$sets, c("window_convective", "convective"))
  expect_null(measured_reference(ms, "BE-Bra"))
  expect_equal(measured_reference(ms, "FI-Lom")$src, "meteoswiss_cabam_mlh")
  expect_true(measured_reference(ms, "DE-HoH")$borderline)
  expect_true(measured_reference(ms, "FR-Lus")$borderline)
  expect_false(measured_reference(ms, "DE-RuS")$borderline)
})

std_rows <- function(from, n, step = 1800, pblh = 500) {
  t <- as.POSIXct(from, tz = "UTC") + step * (seq_len(n) - 1)
  dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t + step, pblh_m = rep(pblh, n))
}

test_that("comparison_period() reports an inverted overlap instead of returning it (FI-Lom)", {
  # FI-Lom on the cluster: PBLHmod valid from 2024-11-04 22:00, CABAM valid until 2024-05-15 19:30
  mod <- std_rows("2023-12-31 22:00", 48 * 900)
  mod$pblh_m[mod$timestamp_start_utc < as.POSIXct("2024-11-04 22:00", tz = "UTC")] <- NA
  cab <- std_rows("2024-01-01 00:00", 48 * 135)
  cp <- comparison_period(pblhmod_default = mod, meteoswiss_cabam_mlh = cab)
  expect_false(cp$ok)
  expect_true(is.na(cp$start) && is.na(cp$end))
  expect_match(cp$reason, "no common period")
  expect_match(cp$reason, "pblhmod_default 2024-11-04 22:00 to")
  expect_equal(cp$ranges$first_valid_utc[1], as.POSIXct("2024-11-04 22:00", tz = "UTC"))
  expect_equal(cp$ranges$n_valid[2], 48 * 135)
  # an input without any valid value
  cab$pblh_m <- NA
  cp0 <- comparison_period(pblhmod_default = mod, meteoswiss_cabam_mlh = cab)
  expect_false(cp0$ok)
  expect_match(cp0$reason, "no valid value in meteoswiss_cabam_mlh")
  # a normal overlap: the range of the old scripts (max of the first, min of the last valid time)
  ok <- comparison_period(a = std_rows("2024-01-01", 100), b = std_rows("2024-01-01 10:00", 100))
  expect_true(ok$ok)
  expect_equal(ok$start, as.POSIXct("2024-01-01 10:00", tz = "UTC"))
  expect_equal(ok$end, as.POSIXct("2024-01-01", tz = "UTC") + 1800 * 99)
  expect_true(is.na(ok$reason))
  expect_error(comparison_period(std_rows("2024-01-01", 2)), "named")
})

test_that("measured_step_s() detects the step and refuses an empty table with a clear error", {
  expect_equal(measured_step_s(std_rows("2024-01-01", 10)), 1800)
  expect_equal(measured_step_s(std_rows("2024-01-01", 10, step = 3600)), 3600)
  expect_error(measured_step_s(std_rows("2024-01-01", 10)[0, ]), "cannot be detected")
  # the old failure: "missing value where TRUE/FALSE needed" in measured_time_base() on an empty cut
  expect_error(measured_time_base(std_rows("2024-01-01", 10)[0, ], NULL, NULL), "cannot be detected")
  expect_error(measured_time_base(std_rows("2024-01-01", 10)[0, ], NULL, NULL, step_s = 900), "not supported")
})

pair_rows <- function(n, regime = "convective") {
  t <- as.POSIXct("2024-04-01 10:00", tz = "UTC") + 1800 * seq_len(n)
  pair_table(dplyr::tibble(timestamp_utc = t, pblh_m = seq(300, by = 50, length.out = n), x_regime = regime),
             dplyr::tibble(timestamp_utc = t, pblh_m = seq(400, by = 40, length.out = n)), "1_x")
}

test_that("an empty stratum set writes the 'all' row with n = 0 and the reason (FI-Lom)", {
  e <- pair_rows(4, "stable")
  e <- e[e$regime %in% "convective", ]
  m <- pblh_metrics_long_set(e, "1_x", strata = c("season", "hour_utc"), empty_reason = "no window interval")
  expect_equal(unique(m$pair), "1_x")
  expect_equal(unique(m$stratum), "all")
  expect_equal(m$value[m$metric == "n"], 0)
  expect_true(all(is.na(m$value[m$metric == "median_log_ratio"])))
  expect_equal(unique(m$no_data_reason), "no window interval")
  expect_setequal(m$metric, names(pblh_metrics(1, 1)))
  # the caller's step of scripts/10_measured_compare.R works on it
  expect_no_error(dplyr::mutate(m, stratum_set = "window_convective", .after = "pair"))
})

test_that("a non-empty set gives the metrics of pblh_metrics_long() and no reason", {
  pt <- pair_rows(6)
  m <- pblh_metrics_long_set(pt, "1_x", strata = c("season", "hour_utc"))
  ref <- pblh_metrics_long(pt, strata = c("season", "hour_utc"))
  expect_equal(m[names(ref)], ref)
  expect_true(all(is.na(m$no_data_reason)))
})

# Wide metrics as built by scripts/11_gamma_evaluation.R (only the columns the tables use)
wide_row <- function(site_id, case, reference, set, n, n_days, mlr) {
  dplyr::tibble(site_id = site_id, case = case, reference = reference, stratum_set = set, n = n, n_days = n_days,
                median_log_ratio = mlr, mlr_lo = mlr - 0.1, mlr_hi = mlr + 0.1, frac_within_factor2 = 0.8, fb = 0.1,
                nrmsd = 0.3, r_spearman = 0.5, slope_sym = 0.9)
}

test_that("gamma_decision_table() gives one row per site with the per-case columns", {
  src <- "meteoswiss_cabam_mlh"
  w <- dplyr::bind_rows(
    wide_row("XX-Aer", "default", src, "window_convective", 500, 60, -0.4),
    wide_row("XX-Aer", "default", src, "window_convective:JJA", 300, 30, -0.6),
    wide_row("XX-Aer", "default", src, "window_convective:DJF", 5, 3, -1.5), # low n: not the season maximum
    wide_row("XX-Aer", "default", src, "window_convective:MAM", 195, 27, -0.2),
    wide_row("XX-Aer", "default", src, "convective", 900, 80, -0.3),
    wide_row("XX-Aer", "default", "era5_blh_center", "convective", 2000, 300, -0.25),
    wide_row("XX-Aer", "gamma_med", src, "window_convective", 500, 60, -0.05),
    wide_row("XX-Aer", "gamma_med", "era5_blh_center", "convective", 2000, 300, 0.1),
    wide_row("XX-Aer", "era5_blh_center", src, "window_convective", 250, 55, 0.02),
    # FI-Lom-like: empty window set
    wide_row("XX-Low", "default", src, "window_convective", 0, 0, NA_real_),
    wide_row("XX-Low", "gamma_med", src, "window_convective", 0, 0, NA_real_),
    wide_row("XX-Era", "default", "era5_blh_center", "convective", 1500, 200, -0.5),
    wide_row("XX-Era", "gamma_med", "era5_blh_center", "convective", 1500, 200, -0.1)
  )
  ms <- ms_tab(ms_row(), ms_row(icos_site = "XX-Low", provider_site = "LLL", borderline = TRUE))
  refs <- list(`XX-Aer` = measured_reference(ms, "XX-Aer"), `XX-Low` = measured_reference(ms, "XX-Low"),
               `XX-Era` = NULL)
  d <- gamma_decision_table(w, c("XX-Era", "XX-Aer", "XX-Low"), refs, c("default", "gamma_med", "era5_blh_center"),
                            min_days = 10)
  expect_equal(d$site_id, c("XX-Era", "XX-Aer", "XX-Low"))
  expect_equal(nrow(d), 3)
  cols <- c("site_id", "reference", "method_class", "stratum_set", "distance_km", "distance_class", "borderline",
            paste0(rep(c("n_", "n_days_", "mlr_", "mlr_lo_", "mlr_hi_", "low_n_", "season_max_", "season_max_mlr_",
                         "season_max_n_days_"), 3), rep(c("default", "gamma_med", "era5_blh_center"), each = 9)),
            "era5_mlr_default", "era5_n_default", "era5_mlr_gamma_med", "era5_n_gamma_med")
  expect_setequal(names(d), cols)
  a <- d[d$site_id == "XX-Aer", ]
  expect_equal(a$mlr_default, -0.4)
  expect_equal(c(a$mlr_lo_default, a$mlr_hi_default), c(-0.5, -0.3))
  expect_equal(a$season_max_default, "JJA")
  expect_equal(a$season_max_mlr_default, -0.6)
  expect_false(a$low_n_default)
  expect_equal(a$era5_mlr_default, -0.25)
  expect_equal(a$mlr_era5_blh_center, 0.02)
  expect_equal(a$stratum_set, "window_convective")
  l <- d[d$site_id == "XX-Low", ]
  expect_equal(l$n_default, 0L)
  expect_true(is.na(l$mlr_default))
  expect_true(l$low_n_default)
  expect_true(l$borderline)
  expect_true(is.na(l$season_max_default))
  e <- d[d$site_id == "XX-Era", ]
  expect_true(is.na(e$reference))
  expect_true(is.na(e$mlr_default))
  expect_equal(e$era5_mlr_gamma_med, -0.1)
})

test_that("gamma_season_table() splits the season sets and flags low n", {
  w <- dplyr::bind_rows(wide_row("XX-Aer", "default", "s", "window_convective", 500, 60, -0.4),
                        wide_row("XX-Aer", "default", "s", "window_convective:JJA", 300, 30, -0.6),
                        wide_row("XX-Aer", "default", "s", "window_convective:DJF", 5, 3, -1.5))
  s <- gamma_season_table(w, min_days = 10)
  expect_equal(nrow(s), 2)
  expect_equal(s$stratum_set, c("window_convective", "window_convective"))
  expect_equal(s$season, c("JJA", "DJF"))
  expect_equal(s$low_n, c(FALSE, TRUE))
})
