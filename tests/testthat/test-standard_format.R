utc <- function(x) as.POSIXct(x, tz = "UTC")

hh_table <- function(n = 4) {
  s <- utc("2020-06-01 00:00:00") + 1800 * (seq_len(n) - 1)
  dplyr::tibble(
    timestamp_start_utc = s,
    timestamp_end_utc = s + 1800,
    pblh_m = c(250.5, NA, 1200.25, 80)[seq_len(n)],
    source_id = "pblhmod_default",
    qc_flag = c(0L, 1L, 0L, 3L)[seq_len(n)],
    x_regime = c("stable", NA, "convective", "stable")[seq_len(n)]
  )
}
hh_meta <- function() {
  standard_meta("XX-Tst", "pblhmod_default", "parameterisation",
                time_reference = "interval_mean", qc_flag_meaning = list(`0` = "no flag", `1` = "interpolated"))
}

test_that("write/read round trip keeps values, UTC times, sidecar", {
  dir <- withr::local_tempdir()
  df <- hh_table()
  paths <- write_standard_pblh(df, hh_meta(), dir)
  expect_equal(unname(basename(paths)), c("XX-Tst__pblhmod_default.csv.gz", "XX-Tst__pblhmod_default.yml"))
  first <- readLines(gzfile(paths[["data"]]), n = 2)
  expect_match(first[2], "^2020-06-01T00:00:00Z,2020-06-01T00:30:00Z,250.5,")

  back <- read_standard_pblh(paths[["data"]])
  expect_equal(back$timestamp_start_utc, df$timestamp_start_utc)
  expect_equal(attr(back$timestamp_end_utc, "tzone"), "UTC")
  expect_equal(back$pblh_m, df$pblh_m)
  expect_identical(back$qc_flag, df$qc_flag)
  expect_equal(back$x_regime, df$x_regime)
  expect_true(all(is.na(back$uncertainty_m)) && all(is.na(back$n_samples)))
  expect_equal(names(back), c(standard_required_cols, standard_recommended_cols, "x_regime"))
  m <- attr(back, "meta")
  expect_equal(m$method_class, "parameterisation")
  expect_equal(m$qc_flag_meaning$`1`, "interpolated")
  expect_equal(m$instrument, "unknown")
  expect_match(m$created_utc, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$")
  expect_true(nzchar(m$code_version))
  expect_equal(m$metadata_completeness, round(2 / 12, 4))
})

test_that("metadata_completeness: 12 optional fields, automatic fields excluded", {
  m <- standard_meta("XX-Tst", "s", "model")
  expect_equal(metadata_completeness(m), 0)
  m$licence <- "CC BY 4.0"
  m$citation <- "Someone et al. 2020"
  m$instrument_altitude_m <- 120
  m$created_utc <- "2026-10-17T00:00:00Z" # not counted
  m$code_version <- "abc1234"             # not counted
  expect_equal(metadata_completeness(m), 3 / 12)
  m$contact <- "UNKNOWN"
  expect_equal(metadata_completeness(m), 3 / 12)
  expect_length(standard_optional_fields(), 12)
  expect_error(standard_meta("XX-Tst", "s", "model", colour = "red"), "Unknown sidecar")
})

test_that("instantaneous rows with start == end are accepted", {
  t <- utc("2020-06-01 00:00:00") + 3600 * 0:2
  df <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t, pblh_m = c(300, 400, 500),
                      source_id = "era5_blh_center", qc_flag = NA_integer_)
  meta <- standard_meta("XX-Tst", "era5_blh_center", "model", time_reference = "instantaneous")
  expect_true(validate_standard_pblh(df, meta))
  dir <- withr::local_tempdir()
  back <- read_standard_pblh(write_standard_pblh(df, meta, dir)[["data"]])
  expect_equal(back$timestamp_start_utc, back$timestamp_end_utc)
  expect_true(all(is.na(back$qc_flag)))
})

test_that("validation lists all problems", {
  df <- hh_table()
  df$source_id <- NULL
  expect_error(validate_standard_pblh(df, hh_meta()), "missing required column")

  df <- hh_table()
  df$timestamp_end_utc[2] <- df$timestamp_start_utc[2] - 60
  df$timestamp_start_utc[4] <- df$timestamp_start_utc[3]
  df$timestamp_end_utc[4] <- df$timestamp_end_utc[3]
  df$pblh_m[1] <- Inf
  df$regime <- "x"
  meta <- hh_meta()
  meta$method_class <- "lidar"
  err <- tryCatch(validate_standard_pblh(df, meta), error = conditionMessage)
  for (pat in c("start > end", "duplicated start", "infinite", "without 'x_' prefix", "invalid method_class")) {
    expect_match(err, pat, fixed = TRUE)
  }

  t <- utc("2020-06-01 00:00:00")
  inst <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t + 1800, pblh_m = 1, source_id = "e")
  expect_error(validate_standard_pblh(inst, standard_meta("XX-Tst", "e", "model", time_reference = "instantaneous")),
               "start != end")
  expect_error(validate_standard_pblh(hh_table(), standard_meta("XX-Tst", "other", "model")), "differs from the sidecar")
  expect_error(validate_standard_pblh(hh_table(), list(site_id = "XX-Tst", source_id = "pblhmod_default")),
               "method_class")
})

test_that("read aborts without sidecar", {
  dir <- withr::local_tempdir()
  p <- write_standard_pblh(hh_table(), hh_meta(), dir)
  file.remove(p[["meta"]])
  expect_error(read_standard_pblh(p[["data"]]), "Sidecar")
})
