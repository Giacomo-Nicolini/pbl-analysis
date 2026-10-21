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
  # rejected u* alone (left NA, not interpolated) gives qc 1
  out$flag_ustar_implausible <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  r3 <- suppressMessages(standardize_pblhmod(out, "default"))
  expect_equal(r3$data$qc_flag, c(1L, 1L, 3L, 2L, 2L))
  expect_match(r3$meta$qc_flag_meaning$`1`, "flag_ustar_implausible", fixed = TRUE)
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

# Synthetic operational series, 2020-01-01 to 2023-06 (local standard time):
# a long run 2020-2022 and a short run in 2023, with noise, floor rows, NA and
# - 1234.56 m: long-run fill, 200 rows in each of 2020-2022 (> 1 % per year);
# - 682.86 m: short-run fill, 80 rows in 2023 (1.08 %) + 2 coincidences in 2021;
#   whole-record share 0.14 %, missed by a whole-record 0.5 % rule;
# - 999.99 m: 20 rows per year 2020-2022 (total 60 >= 50, every year < 0.5 %): not a fill;
# - 333.33 m: 40 rows in 2023 (0.54 %, but total < 50): not a fill.
oper_series <- function() {
  set.seed(42)
  s <- utc("2020-01-01 00:00:00") + 1800 * (seq_len(60000) - 1) # read as local standard time
  year <- as.integer(format(s, "%Y"))
  x <- round(exp(stats::rnorm(length(s), log(600), 0.8)), 2)
  x[x %in% c(50, 1234.56, 682.86, 999.99, 333.33)] <- 600.01 # noise must not hit the test constants
  pick <- function(y, k, used) sample(setdiff(which(year == y), used), k)
  used <- integer(0)
  put <- function(y, k, v) {
    i <- pick(y, k, used)
    used <<- c(used, i)
    x[i] <<- v
  }
  for (y in 2020:2022) put(y, 200, 1234.56)
  put(2023, 80, 682.86)
  put(2021, 2, 682.86)
  for (y in 2020:2022) put(y, 20, 999.99)
  put(2023, 40, 333.33)
  put(2020, 200, 50)
  put(2022, 100, NA)
  list(s = s, year = year, x = x)
}

test_that("detect_fill_values: per-year rule finds long- and short-run fills only", {
  o <- oper_series()
  fv <- detect_fill_values(o$x, o$year)
  expect_equal(fv$value, c(1234.56, 682.86))
  expect_equal(fv$n_total, c(600L, 82L))
  expect_equal(fv$years, list(2020:2022, 2023L))
  expect_equal(fv$n_flagged, c(600L, 80L))
  expect_equal(fv$n_outside, c(0L, 2L))
  n2023 <- sum(o$year == 2023 & !is.na(o$x))
  expect_equal(fv$share[[2]], 80 / n2023)
  # the whole-record share of the short-run fill is below 0.5 %
  expect_lt(82 / sum(!is.na(o$x)), 0.005)
  expect_equal(nrow(detect_fill_values(c(NA, NA), c(2020L, 2020L))), 0)
  expect_error(detect_fill_values(1:3, 2020L), "length")
})

test_that("standardize_icos_oper: qc 0/1/2/NA by value and year, excluded rows dropped", {
  o <- oper_series()
  n <- length(o$x)
  s_utc <- o$s - 3600 # site UTC+1
  d <- dplyr::tibble(site_id = "XX-Tst", timestamp_start_utc = s_utc, timestamp_end_utc = s_utc + 1800,
                     TIMESTAMP_START = format(o$s, "%Y%m%d%H%M"), PBLH = o$x, src_fluxes = TRUE,
                     excl_flag = FALSE)
  d$excl_flag[n] <- TRUE
  d$src_fluxes[n - 1] <- FALSE
  expect_message(r <- standardize_icos_oper(d), "1 of 59999 rows in exclusion periods dropped")
  expect_equal(nrow(r$data), n - 2)
  expect_equal(r$n_floor, 200)
  q <- r$data$qc_flag
  x <- r$data$pblh_m
  yr <- o$year[seq_len(n - 2)]
  expect_equal(sum(q == 1, na.rm = TRUE), 200)
  expect_equal(sum(q == 2, na.rm = TRUE), 680)
  expect_equal(sum(is.na(q)), 100)
  expect_true(all(q[x %in% 682.86 & yr == 2021] == 0)) # coincidences keep qc 0
  expect_true(all(q[x %in% c(999.99, 333.33)] == 0))
  expect_match(r$meta$processing_history,
               "682.86 m: 82 rows in total; years 2023 (1.08 %): 80 rows flagged; 2 rows outside those years",
               fixed = TRUE)
  expect_equal(r$meta$licence, "CC BY 4.0")
  expect_true(validate_standard_pblh(r$data, r$meta))
  d$TIMESTAMP_START <- NULL
  expect_error(suppressMessages(standardize_icos_oper(d)), "TIMESTAMP_START")
})
