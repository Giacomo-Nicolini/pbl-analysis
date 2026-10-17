utc <- function(x) as.POSIXct(x, tz = "UTC")

test_that("standardize_pblhmod maps flags to qc (highest wins) and drops excluded rows", {
  s <- utc("2020-06-01 00:00:00") + 1800 * 0:5
  out <- dplyr::tibble(
    site_id = "XX-Tst", timestamp_start_utc = s, timestamp_end_utc = s + 1800,
    pblh_m = c(100, 200, 20, 400, 500, 600),
    regime = c("stable", "stable", "stable", "convective", "convective", "excluded"),
    method = c("B1", "B1", "B1", "B5", "B5", NA), H_source = "H_F_MDS",
    flag_ustar_interp = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
    flag_L_interp = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    flag_H_interp = FALSE,
    flag_sign_inconsistent = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
    flag_below_zm = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
    flag_no_h0 = TRUE, # not mapped
    excl_flag = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_message(r <- standardize_pblhmod(out, "default"), "1 of 6 rows in exclusion periods dropped")
  expect_equal(r$n_excluded, 1)
  expect_equal(r$data$qc_flag, c(0L, 1L, 3L, 2L, 2L))
  expect_equal(r$data$x_regime, c("stable", "stable", "stable", "convective", "convective"))
  expect_equal(unique(r$data$source_id), "pblhmod_default")
  expect_equal(r$meta$method_class, "parameterisation")
  expect_equal(r$meta$time_reference, "interval_mean")
  expect_true(validate_standard_pblh(r$data, r$meta))
  expect_match(r$meta$processing_history, "options: unknown")
  r2 <- suppressMessages(standardize_pblhmod(out, "gamma0025", opts = pblh_options(gamma = 0.0025)))
  expect_match(r2$meta$processing_history, "gamma = 0.0025", fixed = TRUE)
})

test_that("standardize_pblhmod accepts the derived csv (character ISO times)", {
  s <- utc("2020-06-01 00:00:00") + 1800 * 0:1
  out <- dplyr::tibble(
    site_id = "XX-Tst", timestamp_start_utc = format(s, "%Y-%m-%dT%H:%M:%SZ"),
    timestamp_end_utc = format(s + 1800, "%Y-%m-%dT%H:%M:%SZ"), pblh_m = c(1, 2),
    regime = "stable", method = "B1", H_source = "H", excl_flag = FALSE
  )
  r <- suppressMessages(standardize_pblhmod(out, "default"))
  expect_equal(r$data$timestamp_start_utc, s)
  expect_equal(r$data$qc_flag, c(0L, 0L))
})

test_that("standardize_era5_blh: instantaneous rows, spread extras, qc NA, exclusions dropped", {
  ex <- read_exclusions(fixture_path("exclusion_periods.csv")) # XX-Nor 2022-10-01..2023-04-20 local (UTC+1)
  t <- utc("2022-09-30 21:00:00") + 3600 * 0:4                # 22:00..02:00 local
  center <- dplyr::tibble(site_id = "XX-Nor", cell_row = 0L, cell_col = 0L, timestamp_utc = t,
                          blh_m = c(300, 310, 320, 330, 340))
  bilin <- dplyr::tibble(timestamp_utc = t, blh_m_bilin = c(301, 311, 321, 331, 341))
  spread <- dplyr::tibble(timestamp_utc = t, mean9 = 1, sd9 = 10, min9 = 1, max9 = 2, range9 = 30,
                          cv9 = 0.1, n9 = 9L)
  expect_message(r <- standardize_era5_blh(center, bilin, spread, "center", ex, utc_offset_h = 1),
                 "3 of 5 rows in exclusion periods dropped")
  # 2022-10-01 00:00 local = 2022-09-30 23:00 UTC: first excluded hour
  expect_equal(r$data$timestamp_start_utc, t[1:2])
  expect_equal(r$data$timestamp_end_utc, r$data$timestamp_start_utc)
  expect_equal(r$data$pblh_m, c(300, 310))
  expect_true(all(is.na(r$data$qc_flag)))
  expect_equal(r$data$x_blh_sd9_m, c(10, 10))
  expect_equal(r$data$x_n9, c(9L, 9L))
  expect_equal(r$meta$source_id, "era5_blh_center")
  expect_equal(r$meta$time_reference, "instantaneous")
  expect_equal(r$meta$licence, "Licence to use Copernicus Products")
  expect_match(r$meta$citation, "10.1002/qj.3803", fixed = TRUE)
  expect_match(r$meta$citation, "10.24381/cds.bd0915c6", fixed = TRUE)
  expect_true(validate_standard_pblh(r$data, r$meta))

  rb <- suppressMessages(standardize_era5_blh(center, bilin, spread, "bilinear", ex, 1))
  expect_equal(rb$data$pblh_m, c(301, 311))
  expect_equal(unique(rb$data$source_id), "era5_blh_bilinear")
})

test_that("ERA5 provenance comes from the file header", {
  f <- suppressMessages(find_era5_files("XX-Tst", list(era5_sl_dir = file.path(normalizePath(fixture_path("era5_2026", "sl")), "{site}")),
                                        "single_level", "boundary_layer_height"))
  h <- f$header[[which(f$cell_row == 0 & f$cell_col == 0)]]
  txt <- era5_provenance(h)
  expect_match(txt, "era5_download_arco.py", fixed = TRUE)
  expect_match(txt, "3x3 cells", fixed = TRUE)
  expect_match(txt, "cell centre 50 N, 10 E", fixed = TRUE)
})

# Synthetic operational series: floor rows, two fill constants, a frequent but
# too rare value, realistic noise.
oper_series <- function() {
  set.seed(42)
  n <- 20000
  x <- round(exp(stats::rnorm(n, log(600), 0.8)), 2)
  x[x %in% c(50, 1234.56, 777.77, 999.99)] <- 600.01 # noise must not hit the test constants
  x[1:200] <- 50          # floor
  x[201:500] <- 1234.56   # fill 1: 300 rows
  x[501:650] <- 777.77    # fill 2: 150 rows
  x[651:710] <- 999.99    # 60 rows: >= 50 but < 0.5 % of non-missing -> not a fill
  x[711:810] <- NA
  x
}

test_that("detect_fill_values finds exactly the constants", {
  x <- oper_series()
  fv <- detect_fill_values(x)
  expect_equal(fv$value, c(1234.56, 777.77))
  expect_equal(fv$n, c(300L, 150L))
  expect_equal(fv$frac, c(300, 150) / sum(!is.na(x)))
  expect_equal(nrow(detect_fill_values(c(NA, NA))), 0)
})

test_that("standardize_icos_oper: qc 0/1/2/NA, excluded rows dropped", {
  x <- oper_series()
  n <- length(x)
  s <- utc("2020-01-01 00:00:00") + 1800 * (seq_len(n) - 1)
  d <- dplyr::tibble(site_id = "XX-Tst", timestamp_start_utc = s, timestamp_end_utc = s + 1800,
                     PBLH = x, src_fluxes = TRUE, excl_flag = FALSE)
  d$excl_flag[n] <- TRUE
  d$src_fluxes[n - 1] <- FALSE
  expect_message(r <- standardize_icos_oper(d), "1 of 19999 rows in exclusion periods dropped")
  expect_equal(nrow(r$data), n - 2)
  expect_equal(r$n_floor, 200)
  q <- r$data$qc_flag
  expect_equal(sum(q == 1, na.rm = TRUE), 200)
  expect_equal(sum(q == 2, na.rm = TRUE), 450)
  expect_equal(sum(is.na(q)), 100)
  expect_true(all(q[651:710] == 0))
  expect_equal(r$fill_values$value, c(1234.56, 777.77))
  expect_match(r$meta$processing_history, "1234.56 m (300 rows)", fixed = TRUE)
  expect_equal(r$meta$licence, "CC BY 4.0")
  expect_true(validate_standard_pblh(r$data, r$meta))
})
