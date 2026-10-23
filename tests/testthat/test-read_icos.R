utc <- function(x) as.POSIXct(x, tz = "UTC")

l2_dir <- normalizePath(fixture_path("icos_l2"))
site_dir <- file.path(l2_dir, "XX-Tst_output_L2I_2026")
cfg <- list(icos_l2_dir = l2_dir)
f_fluxes <- file.path(site_dir, "ICOSETC_XX-Tst_FLUXES_INTERIM_L2.zip")
f_fluxnet <- file.path(site_dir, "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2020_v1.0_r1.zip")
f_meteo <- file.path(site_dir, "ICOSETC_XX-Tst_METEO_INTERIM_L2.zip")

# Plain CSV in a temp file, for edge cases (readr reads csv and zip alike).
write_l2 <- function(start, end) {
  f <- withr::local_tempfile(fileext = ".csv", .local_envir = parent.frame())
  df <- data.frame(TIMESTAMP_START = start, TIMESTAMP_END = end, H = seq_along(start))
  utils::write.csv(df, f, row.names = FALSE, quote = FALSE)
  f
}
q <- function(expr) suppressMessages(expr)

test_that("local standard time -> UTC at different offsets", {
  expect_equal(icos_time_to_utc("202001010000", 1), utc("2019-12-31 23:00:00"))
  expect_equal(icos_time_to_utc("202001010000", 2), utc("2019-12-31 22:00:00"))
  expect_equal(icos_time_to_utc("202001010000", 0), utc("2020-01-01 00:00:00"))
  expect_equal(attr(icos_time_to_utc("202001010000", 1), "tzone"), "UTC")
})

test_that("FLUXES: -9999 -> NA, degC -> K, grid completion adds exactly the missing row", {
  expect_message(
    expect_warning(d <- read_icos_fluxes(f_fluxes, utc_offset_h = 1), "1 step"),
    "1 row added"
  )
  expect_equal(nrow(d), 7)
  expect_equal(d$timestamp_start_utc[1], utc("2019-12-31 23:00:00"))
  expect_equal(d$timestamp_end_utc[1], utc("2019-12-31 23:30:00"))
  expect_true(all(diff(as.numeric(d$timestamp_start_utc)) == 1800))

  # the missing 01:00 local record (= 00:00 UTC) is an added row of NA
  i <- which(d$timestamp_start_utc == utc("2020-01-01 00:00:00"))
  expect_false(d$src_fluxes[i])
  expect_true(is.na(d$TIMESTAMP_START[i]))
  expect_true(is.na(d$H[i]))
  expect_equal(d$timestamp_end_utc[i], utc("2020-01-01 00:30:00"))
  expect_equal(sum(d$src_fluxes), 6)

  expect_true(is.na(d$H[1]))     # -9999
  expect_true(is.na(d$PBLH[2]))  # -9999.0
  expect_equal(d$TIMESTAMP_START[1], "202001010000")
  expect_equal(d$T_SONIC_K[1], 283.15)
  expect_false("T_SONIC" %in% names(d))
  expect_false("CO2" %in% names(d)) # not requested
  info <- attr(d, "read_info")
  expect_equal(info$n_rows_added, 1)
  expect_equal(info$n_rows_file, 6)
})

test_that("FLUXNET: kPa -> Pa, prefix fn_, QC columns untouched", {
  d <- q(read_icos_fluxnet(f_fluxnet, 1))
  expect_equal(d$fn_TA_F_K[1], 278.15)
  expect_equal(d$fn_PA_F_PA[1], 100500)
  expect_true("fn_TA_F_QC" %in% names(d))
  expect_true(is.na(d$fn_H_CORR[4]))
  expect_false(any(c("TA_F", "fn_TA_F", "fn_PA_F", "fn_TA_F_MDS") %in% names(d)))
  expect_true(all(c("TIMESTAMP_START", "timestamp_start_utc", "src_fluxnet") %in% names(d)))
})

test_that("METEO: TA profile levels detected from the header, no warning for absent levels", {
  expect_no_warning(d <- q(read_icos_meteo(f_meteo, 1)))
  expect_equal(attr(d, "read_info")$optional_cols, paste0("TA_", 1:9, collapse = ", "))
  expect_equal(attr(d, "read_info")$missing_cols, "")
  expect_equal(d$met_TA_9_K[1], 1.9 + 273.15)
  expect_equal(nrow(d), 8)
  expect_true(is.na(d$met_TA_K[3])) # -9999.0
  expect_equal(d$met_TA_K[1], 275.15)
  expect_equal(d$met_PA_PA[1], 99500)
  expect_equal(d$met_TA_1_K[1], 1.1 + 273.15)
  expect_false("met_TA_10_K" %in% names(d))
})

test_that("duplicated timestamps abort", {
  f <- write_l2(c("202001010000", "202001010000"), c("202001010030", "202001010030"))
  expect_error(q(read_icos_csv(f, c("TIMESTAMP_START", "TIMESTAMP_END", "H"), 1)), "duplicated")
})

test_that("non-monotonic order is sorted with a warning", {
  f <- write_l2(c("202001010030", "202001010000"), c("202001010100", "202001010030"))
  expect_warning(d <- q(read_icos_csv(f, c("TIMESTAMP_START", "TIMESTAMP_END", "H"), 1)), "order")
  expect_equal(d$H, c(2, 1))
})

test_that("step and duration != 30 min are reported", {
  f <- write_l2(c("202001010000", "202001010100"), c("202001010030", "202001010130"))
  expect_warning(
    d <- q(read_icos_csv(f, c("TIMESTAMP_START", "TIMESTAMP_END", "H"), 1)),
    "1 step"
  )
  expect_equal(nrow(d), 3)
  f <- write_l2(c("202001010000", "202001010030"), c("202001010030", "202001010200"))
  expect_warning(q(read_icos_csv(f, c("TIMESTAMP_START", "TIMESTAMP_END", "H"), 1)), "TIMESTAMP_END")
  f <- write_l2(c("202001010000", "202001010045"), c("202001010030", "202001010115"))
  expect_error(suppressWarnings(q(
    read_icos_csv(f, c("TIMESTAMP_START", "TIMESTAMP_END", "H"), 1)
  )), "grid")
})

test_that("find_icos_files: required, optional and ambiguous files", {
  f <- find_icos_files("XX-Tst", cfg)
  expect_equal(basename(f$fluxnet), basename(f_fluxnet))
  expect_false(is.null(f$meteo))

  tmp <- withr::local_tempdir()
  d <- file.path(tmp, "XX-Tst_output_L2I_2026")
  dir.create(d)
  file.copy(f_fluxes, d)
  w <- testthat::capture_warnings(f2 <- find_icos_files("XX-Tst", list(icos_l2_dir = tmp)))
  expect_length(w, 2)
  expect_match(w[1], "FLUXNET")
  expect_match(w[2], "METEO")
  expect_null(f2$fluxnet)

  file.copy(f_fluxnet, file.path(d, "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2021_v1.0_r1.zip"))
  file.copy(f_fluxnet, file.path(d, "ICOSETC_XX-Tst_FLUXNET_FLUXMET_HH_2020-2022_v1.0_r1.zip"))
  expect_error(suppressWarnings(find_icos_files("XX-Tst", list(icos_l2_dir = tmp))), "More than one")

  unlink(file.path(d, basename(f_fluxes)))
  expect_error(suppressWarnings(find_icos_files("XX-Tst", list(icos_l2_dir = tmp))), "FLUXES")
})

test_that("read_icos_site merges the three files on a complete grid", {
  sites <- read_sites(fixture_path("sites.csv"))
  sites_tv <- read_sites_timevarying(fixture_path("sites_timevarying.csv"))
  ex <- read_exclusions(fixture_path("exclusion_periods.csv"))
  d <- suppressWarnings(q(read_icos_site("XX-Tst", cfg, sites, sites_tv, ex)))

  # FLUXES range 00:00-03:00 local (reference period): METEO 03:30 dropped
  expect_equal(nrow(d), 7)
  expect_equal(attr(d, "reference_period_utc"), utc(c("2019-12-31 23:00:00", "2020-01-01 02:00:00")))
  expect_false(anyDuplicated(d$timestamp_start_utc) > 0)
  expect_true(all(diff(as.numeric(d$timestamp_start_utc)) == 1800))
  expect_equal(d$site_id[1], "XX-Tst")
  expect_equal(sum(d$src_fluxes), 6)
  expect_equal(sum(d$src_fluxnet), 4)
  expect_equal(sum(d$src_meteo), 7)
  # 01:00 local is missing in FLUXES: TIMESTAMP_START coalesced from FLUXNET
  i <- which(d$timestamp_start_utc == utc("2020-01-01 00:00:00"))
  expect_equal(d$TIMESTAMP_START[i], "202001010100")
  expect_false(d$src_fluxes[i])
  expect_true(all(c("T_SONIC_K", "fn_TA_F_K", "met_TA_K", "zm_m", "zm_source",
                    "canopy_height_m", "hc_source", "excl_flag", "excl_reason") %in% names(d)))
  expect_false(any(grepl("\\.(fx|fn|met)$", names(d))))
  expect_equal(unique(d$zm_source), "timevarying")
  expect_false(any(d$excl_flag))
  expect_equal(nrow(attr(d, "read_info")), 3)
  expect_equal(attr(d, "read_info")$n_rows_outside_ref, c(0L, 0L, 1L))
})

test_that("icos_reference_period returns the FLUXES range", {
  p <- icos_reference_period("XX-Tst", cfg, 1)
  expect_equal(p, utc(c("2019-12-31 23:00:00", "2020-01-01 02:00:00")))
})

test_that("read_icos_site drops FLUXNET / METEO rows outside the FLUXES range", {
  sites <- read_sites(fixture_path("sites.csv"))
  sites_tv <- read_sites_timevarying(fixture_path("sites_timevarying.csv"))
  ex <- read_exclusions(fixture_path("exclusion_periods.csv"))
  tmp <- withr::local_tempdir()
  d_site <- file.path(tmp, "XX-Tst_output_L2I_2026")
  dir.create(d_site)
  file.copy(c(f_fluxnet, f_meteo), d_site)
  # FLUXES starting at 01:30 local (first two rows removed)
  fl <- readr::read_csv(f_fluxes, col_types = readr::cols(.default = "c"))
  fl <- fl[-(1:2), ]
  csv <- file.path(tmp, "ICOSETC_XX-Tst_FLUXES_INTERIM_L2.csv")
  readr::write_csv(fl, csv)
  utils::zip(file.path(d_site, "ICOSETC_XX-Tst_FLUXES_INTERIM_L2.zip"), csv, flags = "-j9Xq")

  d <- suppressWarnings(q(read_icos_site("XX-Tst", list(icos_l2_dir = tmp), sites, sites_tv, ex)))
  expect_equal(nrow(d), 4) # 01:30-03:00 local
  expect_equal(min(d$timestamp_start_utc), utc("2020-01-01 00:30:00"))
  expect_equal(max(d$timestamp_start_utc), utc("2020-01-01 02:00:00"))
  expect_equal(sum(d$src_fluxnet), 1) # FLUXNET 01:30 only
  # FLUXNET 00:00-01:00 (3 rows) and METEO 00:00-01:00 + 03:30 (4 rows) dropped
  expect_equal(attr(d, "read_info")$n_rows_outside_ref, c(0L, 3L, 4L))
  # inside the period FLUXNET still has priority over METEO for T
  expect_equal(d$fn_TA_F_K[1], 6.5 + 273.15)
})

test_that("icos_range_report counts but never modifies", {
  d <- suppressWarnings(q(read_icos_fluxes(f_fluxes, 1)))
  r <- icos_range_report(d)
  expect_named(r, c("variable", "n_total", "n_excluded", "n_valid", "pct_valid", "n_out_of_range", "min", "max"))
  h <- r[r$variable == "H", ]
  expect_equal(h$n_total, 7)
  expect_equal(h$n_valid, 5)
  expect_equal(h$n_out_of_range, 1) # 900
  expect_equal(r$n_out_of_range[r$variable == "MO_LENGTH"], 1) # L = 0
  expect_equal(r$n_out_of_range[r$variable == "PBLH"], 1)      # 5 m
  expect_equal(r$n_valid[r$variable == "fn_PA_F_PA"], 0)       # absent in FLUXES alone
  expect_equal(max(d$H, na.rm = TRUE), 900)
})

test_that("T_SONIC is already in K: renamed, not converted", {
  d <- suppressWarnings(q(read_icos_fluxes(f_fluxes, 1)))
  expect_equal(d$T_SONIC_K[1], 283.15)
})

test_that("unit sanity check aborts on implausible medians of declared units", {
  write_col <- function(name, values) {
    f <- withr::local_tempfile(fileext = ".csv", .local_envir = parent.frame())
    n <- length(values)
    s <- format(as.POSIXct("2020-01-01", tz = "UTC") + (seq_len(n) - 1) * 1800, "%Y%m%d%H%M")
    e <- format(as.POSIXct("2020-01-01", tz = "UTC") + seq_len(n) * 1800, "%Y%m%d%H%M")
    df <- data.frame(TIMESTAMP_START = s, TIMESTAMP_END = e)
    df[[name]] <- values
    utils::write.csv(df, f, row.names = FALSE, quote = FALSE)
    f
  }
  rd <- function(f, col) q(read_icos_csv(f, c("TIMESTAMP_START", "TIMESTAMP_END", col), 1))

  # T_SONIC declared K: degC values abort, K values pass
  expect_error(rd(write_col("T_SONIC", c(10, 11, 12)), "T_SONIC"), "unit check failed")
  expect_no_error(rd(write_col("T_SONIC", c(283, 284, 285)), "T_SONIC"))
  # TA_F declared degC: K values abort; boundaries [-40, 50] inclusive
  expect_error(rd(write_col("TA_F", c(283, 284, 285)), "TA_F"), "TA_F")
  expect_no_error(rd(write_col("TA_F", c(-40, -40, -40)), "TA_F"))
  expect_error(rd(write_col("TA_F", c(-41, -41, -41)), "TA_F"), "unit check failed")
  # PA declared kPa: Pa values abort
  expect_error(rd(write_col("PA", c(100000, 100100, 100200)), "PA"), "kPa")
  expect_no_error(rd(write_col("PA", c(99.5, 99.6, 99.7)), "PA"))
  # all missing: no check
  expect_no_error(rd(write_col("PA", c(-9999, -9999, -9999)), "PA"))
})

test_that("icos_range_report skips excluded rows by default", {
  d <- suppressWarnings(q(read_icos_fluxes(f_fluxes, 1)))
  d$excl_flag <- FALSE
  d$excl_flag[d$H %in% 900] <- TRUE
  expect_message(r <- icos_range_report(d), "1 excluded")
  h <- r[r$variable == "H", ]
  expect_equal(h$n_total, 6)
  expect_equal(h$n_excluded, 1)
  expect_equal(h$n_out_of_range, 0)
  r2 <- suppressMessages(icos_range_report(d, include_excluded = TRUE))
  expect_equal(r2$n_total[r2$variable == "H"], 7)
  expect_equal(r2$n_out_of_range[r2$variable == "H"], 1)
  expect_equal(r2$n_excluded[r2$variable == "H"], 1)
})
