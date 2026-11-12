# Spec 14: secondary measured sources (DE-Tha SkyVUE PRO built-in MLH, CZ-Lnz BL-View L3) and the explicit column
# types of standard files. Synthetic data only; checks on the real samples are skipped when the samples are absent.

# --- SkyVUE ---------------------------------------------------------------------------------------------------------

# Bitwise CRC-16 exactly as the C code of the SkyVUE PRO manual (Sect. 6.2), for one record.
crc16_manual <- function(bytes) {
  crc <- 0xFFFF
  for (x in bytes) {
    crc <- bitwXor(crc, bitwShiftL(x, 8L))
    for (j in 1:8) {
      m <- if (bitwAnd(crc, 0x8000) != 0) 0x1021 else 0
      crc <- bitwXor(bitwAnd(bitwShiftL(crc, 1L), 0xFFFF), m)
    }
  }
  bitwXor(crc, 0xFFFF)
}

# Synthetic SkyVUE day: n records of message 006 (20 gates) and the matching CSV export.
csv_header <- c("UTC Timestamp", "Cloud Base 1( Meters )", "Cloud Base 2( Meters )", "Cloud Base 3( Meters )",
                "Cloud Base 4( Meters )", "Vertical Visibility( Meters )", "Ceil Status( Status )",
                "Backscatter Profile( (100000*srad*km)^-1 )", paste0("Sky Condition ", 1:5, "( Oktas )"),
                paste0("Sky Condition Height ", 1:5, "( Distance )"), paste0("Mixing Layer ", 1:3, "( Meters )"),
                paste0("Confidence ", 1:3, "( Confidence 1-3 )"), "Highest Signal( Meters )",
                "Ceil Resolution( Meters )")
fmt_e <- function(v) {
  # CSV style: decimal comma, mantissa without trailing zeros, "0" for zero
  s <- ifelse(v == 0, "0", sub("\\.?0+E", "E", formatC(v, format = "E", digits = 3)))
  gsub(".", ",", s, fixed = TRUE)
}
make_skyvue_day <- function(dir, n = 4, t0 = as.POSIXct("2026-07-15 00:00:09", tz = "UTC"), step = 9,
                            mlh = matrix(c(375, 650, NA), n, 3, byrow = TRUE), conf = matrix(c(2L, 1L, 0L), n, 3, byrow = TRUE),
                            cbh = rep(2802, n), corrupt = integer(), status = rep(1L, n), vv = rep(NA_real_, n),
                            layout = "full", vv_csv = vv) {
  set.seed(1)
  ints <- matrix(sample(-200:6000, n * 20, replace = TRUE), n, 20)
  ints[, 3] <- 0L # a zero gate: the CSV parse stops there (ambiguous)
  hex <- apply(ints, 1, function(r) paste(sprintf("%05x", ifelse(r < 0, r + 1048576L, r)), collapse = ""))
  m5 <- function(x, w = 5) ifelse(is.na(x), strrep("/", w), formatC(x, width = w, flag = "0", format = "d"))
  txt <- character(n)
  rows <- character(n)
  t <- t0 + step * (seq_len(n) - 1)
  for (i in seq_len(n)) {
    h1 <- if (status[i] == 5L) vv[i] else if (status[i] %in% 1:4) cbh[i] else NA # line 2: h1 = VV at status 5
    body <- paste0("CS0016006\002\r\n",
                   status[i], "0 096 ", m5(h1), " ///// ///// ///// c00000000000\r\n",
                   "  3 0141  8 0241  0 ////  0 ////  0 ////\r\n",
                   "00100 05 0020 100 +39 00 0096 0020 30 131\r\n",
                   paste(m5(mlh[i, 1]), m5(conf[i, 1]), m5(mlh[i, 2]), m5(conf[i, 2]), m5(mlh[i, 3]), m5(conf[i, 3])),
                   "\r\n", hex[i], "\r\n\003")
    crc <- crc16_manual(as.integer(charToRaw(body)))
    if (i %in% corrupt) crc <- bitwXor(crc, 1L)
    txt[i] <- paste0("\001", body, sprintf("%04x", crc), "\004\r\n")
    prof <- paste(fmt_e(ints[i, ] * 1e-8), collapse = ",")
    q <- function(x) paste0("\"", x, "\"")
    num <- function(x) ifelse(is.na(x), "/////", gsub(".", ",", format(x, nsmall = 0), fixed = TRUE))
    cb_csv <- if (status[i] %in% 1:4) q(paste0(cbh[i], ",00000425904")) else "/////"
    if (layout == "full") {
      rows[i] <- paste(format(t[i], "%d.%m.%Y %H:%M:%S"), cb_csv, "/////", "/////", "/////",
                       num(vv_csv[i]), "c00000000000", q(prof), "3", "8", "0", "0", "0", "1410", "2410", "/////", "/////",
                       "/////", num(mlh[i, 1]), num(mlh[i, 2]), num(mlh[i, 3]), conf[i, 1], conf[i, 2], conf[i, 3],
                       "/////", "5", sep = ",")
    } else { # 7-column layout before the MLH option (DE-Tha until about 2025-06-24)
      rows[i] <- paste(format(t[i], "%d.%m.%Y %H:%M:%S"), cb_csv, "/////", "/////", "/////", num(vv_csv[i]), "/////",
                       sep = ",")
    }
  }
  csv <- file.path(dir, "E1227_150726_Midnight.csv")
  txtf <- file.path(dir, "E1227_150726_Midnight.txt")
  hdr <- if (layout == "full") csv_header else csv_header[c(1:6, 25)] # UTC Timestamp, Cloud Base 1-4, VV, Highest Signal
  writeLines(c(paste(hdr, collapse = ","), rows), csv)
  writeBin(charToRaw(paste(txt, collapse = "")), txtf)
  list(csv = csv, txt = txtf, ints = ints, t = t)
}

test_that("skyvue_num parses decimal commas and the missing code", {
  expect_equal(skyvue_num(c("2802,00000425904", "/////", "375", "", "-2,4E-07", "0")),
               c(2802.00000425904, NA, 375, NA, -2.4e-07, 0))
})

test_that("SkyVUE CRC-16: manual algorithm, table form and the CRC-16/GENIBUS check value", {
  b <- as.integer(charToRaw("123456789"))
  expect_equal(crc16_manual(b), 0xD64E) # published check value of CRC-16/GENIBUS (poly 0x1021, init/xorout 0xFFFF)
  expect_equal(skyvue_crc16(matrix(b, 1)), 0xD64E)
  m <- rbind(c(b, 0L, 0L), c(as.integer(charToRaw("CS0016006")), 1L, 2L))
  expect_equal(skyvue_crc16(m, len = c(9L, 11L)), c(0xD64E, crc16_manual(m[2, ])))
})

test_that("SkyVUE hex profile: 20-bit two's complement times 1e-8, scale undone", {
  p <- skyvue_decode_profile("0027efffe8000000000a", scale = 100)
  expect_equal(p[1, ], c(638e-8, -24e-8, 0, 10e-8))
  expect_equal(skyvue_decode_profile("7ffff80000", 100)[1, ], c(524287e-8, -524288e-8))
  expect_equal(skyvue_decode_profile("0027e", 50)[1, 1], 638e-8 * 2)
  expect_error(skyvue_decode_profile("0027", 100), "multiple of 5")
})

test_that("SkyVUE CSV profile head: decimal-comma values, stops at the ambiguous zero", {
  x <- c("6,38E-06,1,895E-05,2,674E-05", "3,91E-06,-2E-07,0,6E-08,1,8E-07", "0,5E-06")
  h <- skyvue_csv_profile_head(x, 3)
  expect_equal(h[1, ], c(6.38e-06, 1.895e-05, 2.674e-05))
  expect_equal(h[2, ], c(3.91e-06, -2e-07, NA)) # "0,6E-08" is "0" then "6E-08": not parsed
  expect_true(all(is.na(h[3, ])))
})

test_that("read_skyvue_txt splits the messages, parses the fields and checks the CRC", {
  d <- withr::local_tempdir()
  s <- make_skyvue_day(d, n = 3, mlh = rbind(c(375, 650, 1009), c(NA, NA, NA), c(400, NA, NA)),
                       conf = rbind(c(2L, 1L, 1L), c(0L, 0L, 0L), c(3L, 0L, 0L)), corrupt = 3L)
  r <- read_skyvue_txt(s$txt)
  expect_equal(nrow(r), 3)
  expect_equal(r$mlh1_m, c(375, NA, 400))
  expect_equal(r$mlh3_m, c(1009, NA, NA))
  expect_equal(r$mlh1_q, c(2L, 0L, 3L))
  expect_equal(r$h1_m, rep(2802, 3))
  expect_equal(r$detection_status, rep(1L, 3))
  expect_equal(r$window_tr_pct, rep(96, 3))
  expect_equal(r$n_gates, rep(20L, 3))
  expect_equal(r$scale_pct, rep(100, 3))
  expect_equal(r$laser_temp_c, rep(39, 3))
  expect_equal(r$crc_ok, c(TRUE, TRUE, FALSE))
  expect_equal(attr(r, "skyvue_txt_meta")$n_crc_fail, 1L)
  bad <- withr::local_tempfile()
  writeBin(charToRaw("no message here"), bad)
  expect_error(read_skyvue_txt(bad), "SOH")
})

test_that("read_skyvue_day pairs CSV and TXT and fails on a shifted record", {
  d <- withr::local_tempdir()
  mlh <- cbind(c(300, 350, 400, 450), c(900, NA, 950, NA), NA)
  conf <- cbind(c(1L, 2L, 3L, 1L), c(1L, 0L, 1L, 0L), 0L)
  s <- make_skyvue_day(d, n = 4, mlh = mlh, conf = conf)
  r <- read_skyvue_day(s$csv, s$txt)
  expect_equal(r$time_end_utc, s$t)
  expect_equal(r$mlh1_m, mlh[, 1])
  expect_equal(r$mlh2_q, conf[, 2])
  m <- attr(r, "skyvue_meta")
  expect_equal(m$n_records, 4)
  expect_equal(m$n_profile_records_checked, 4)
  expect_lt(m$max_rel_diff_profile, 0.01)
  # a deliberately shifted record: drop the first CSV row and repeat the last one
  l <- readLines(s$csv)
  writeLines(c(l[1], l[3:5], l[5]), s$csv)
  expect_error(read_skyvue_day(s$csv, s$txt), "differ")
  # unequal record counts
  writeLines(l[1:4], s$csv)
  expect_error(read_skyvue_day(s$csv, s$txt), "3 records")
  # fix-14a item 10: the count mismatch has its own class with both counts (skipped day, not a failed check),
  # while a shifted record with equal counts stays an ordinary error
  e <- tryCatch(read_skyvue_day(s$csv, s$txt), error = function(e) e)
  expect_s3_class(e, "skyvue_count_mismatch")
  expect_equal(c(e$n_csv, e$n_txt), c(3L, 4L))
  writeLines(c(l[1], l[3:5], l[5]), s$csv)
  e2 <- tryCatch(read_skyvue_day(s$csv, s$txt), error = function(e) e)
  expect_false(inherits(e2, "skyvue_count_mismatch"))
})

test_that("a record with a CRC mismatch keeps its time but loses its layers", {
  d <- withr::local_tempdir()
  s <- make_skyvue_day(d, n = 3, corrupt = 2L)
  r <- read_skyvue_day(s$csv, s$txt)
  expect_equal(r$crc_ok, c(TRUE, FALSE, TRUE))
  expect_true(is.na(r$mlh1_m[2]) && is.na(r$mlh1_q[2]))
  expect_equal(attr(r, "skyvue_meta")$n_crc_fail, 1L)
})

test_that("skyvue_halfhour: half-hour of the stamp minus 1 s, medians per layer, no layer chosen", {
  t <- as.POSIXct("2025-07-01 00:00:00", tz = "UTC") + c(9, 600, 1800, 1809, 2000)
  r <- dplyr::tibble(time_end_utc = t, crc_ok = c(TRUE, TRUE, TRUE, TRUE, FALSE), detection_status = c(1L, 0L, 0L, 0L, 1L),
                     window_tr_pct = 96, mlh1_m = c(300, 400, 500, 700, 900), mlh1_q = c(1L, 3L, 3L, 2L, 3L),
                     mlh2_m = c(NA, 800, NA, NA, NA), mlh2_q = c(0L, 1L, 0L, 0L, 0L), mlh3_m = NA_real_, mlh3_q = 0L)
  h <- skyvue_halfhour(r)
  expect_equal(nrow(h), 2)
  expect_equal(h$n_samples, c(3L, 1L)) # the record ending at 00:30:00 belongs to 00:00; the CRC failure is not used
  expect_equal(h$x_n_records, c(3L, 2L))
  expect_equal(h$x_layer1_m, c(400, 700))
  expect_equal(h$x_layer1_conf, c(3, 2))
  expect_equal(h$x_layer1_share_conf3, c(2 / 3, 0))
  expect_equal(h$x_layer2_m, c(800, NA))
  expect_equal(h$x_layer3_n, c(0L, 0L))
  s <- suppressMessages(standardize_skyvue_mlh(h, "XX-Tst", read_exclusions(fixture_path("exclusion_periods.csv")), 1))
  expect_true(all(is.na(s$data$pblh_m)))
  expect_equal(s$meta$source_id, "tud_skyvue_mlh")
  expect_equal(s$meta$method_class, "aerosol")
  expect_match(s$meta$licence, "internal use only")
  expect_no_error(validate_standard_pblh(s$data, s$meta))
})

test_that("find_skyvue_days pairs the day files and ignores the .nc", {
  d <- withr::local_tempdir()
  for (f in c("E1227_150726_Midnight.csv", "E1227_150726_Midnight.txt", "E1227_150726_Midnight.nc",
              "E1227_160726_Midnight.csv")) writeLines("x", file.path(d, f))
  x <- suppressMessages(find_skyvue_days(d))
  expect_equal(x$day, as.Date(c("2026-07-15", "2026-07-16")))
  expect_match(x$csv[1], "150726_Midnight.csv$")
  # fix-14c item 9: a day with only one of the two files is kept with a reason, not dropped
  expect_equal(x$skip_reason, c(NA, "no TXT file"))
  expect_true(is.na(x$txt[2]))
})

test_that("short CSV layout (7 columns, before the MLH option): paired by order, cloud base and VV checked", {
  d <- withr::local_tempdir()
  st <- c(1L, 0L, 5L, 1L)
  s <- make_skyvue_day(d, n = 4, status = st, vv = c(NA, NA, 120, NA), layout = "short")
  a <- read_skyvue_csv(s$csv)
  expect_equal(attr(a, "layout"), "short")
  expect_equal(a$cloud_base1_m, c(2802.00000425904, NA, NA, 2802.00000425904))
  expect_equal(a$vertical_visibility_m, c(NA, NA, 120, NA))
  expect_true(all(is.na(a$mlh1_m)) && all(is.na(a$profile_csv)))
  r <- read_skyvue_day(s$csv, s$txt)
  m <- attr(r, "skyvue_meta")
  expect_equal(m$csv_layout, "short")
  expect_equal(m$profile_check, "not available (short CSV layout)")
  expect_true(is.na(m$n_profile_records_checked) && is.na(m$max_rel_diff_profile))
  expect_true(all(is.na(r$mlh1_m)) && all(is.na(r$mlh3_q)))
  expect_equal(r$detection_status, st)
  expect_equal(nrow(r), 4)
  # a vertical visibility that differs between CSV and TXT fails loudly
  s2 <- make_skyvue_day(withr::local_tempdir(), n = 4, status = st, vv = c(NA, NA, 120, NA), vv_csv = c(NA, NA, 150, NA),
                        layout = "short")
  expect_error(read_skyvue_day(s2$csv, s2$txt), "vertical visibility")
  # a shifted record fails (cloud base against a status-0 record)
  l <- readLines(s$csv)
  writeLines(c(l[1], l[3:5], l[5]), s$csv)
  expect_error(read_skyvue_day(s$csv, s$txt), "differ")
  # an unknown layout is refused
  writeLines(c("UTC Timestamp,Something( Meters )", "15.07.2026 00:00:09,1"), s$csv)
  expect_error(read_skyvue_csv(s$csv), "unknown SkyVUE CSV layout")
})

test_that("the real SkyVUE sample passes every pairing and decode check", {
  dd <- here::here("data", "raw", "measured", "DE-Tha", "samples", "07_2026")
  csv <- file.path(dd, "E1227_150726_Midnight.csv")
  skip_if_not(file.exists(csv), "SkyVUE sample not on this machine")
  r <- read_skyvue_day(csv, file.path(dd, "E1227_150726_Midnight.txt"))
  m <- attr(r, "skyvue_meta")
  expect_equal(m$n_records, 10063)
  expect_equal(m$n_crc_fail, 0L)
  expect_equal(m$n_profile_records_checked, 10063)
  expect_lt(m$max_rel_diff_profile, 1e-6)
  expect_equal(m$scale_values, 100)
  expect_equal(unique(r$n_gates), 2048L)
})

# --- BL-View --------------------------------------------------------------------------------------------------------

# Synthetic BL-View 3.2.2 file (as in test-measured_eprofile.R): hourly values repeated on time_aerosol, fill -999.
write_blview_nc <- function(path, tt, h, q) {
  dim_t <- ncdf4::ncdim_def("time_aerosol", "seconds since 1970-01-01 00:00:00.000", tt, unlim = FALSE)
  v <- list(h = ncdf4::ncvar_def("mean_mixing_layer_height", "m", list(dim_t), missval = -999, prec = "integer"),
            q = ncdf4::ncvar_def("mean_mixing_layer_quality", "", list(dim_t), missval = -999, prec = "integer"))
  nc <- ncdf4::nc_create(path, unname(v), force_v4 = TRUE)
  ncdf4::ncvar_put(nc, v$h, h)
  ncdf4::ncvar_put(nc, v$q, q)
  ncdf4::nc_close(nc)
  path
}

# Synthetic BL-View 3.1.x file (MLH as 1 x n, time as an int variable with attribute `Units`).
write_blview31_nc <- function(path, tt, h, q, v312 = FALSE, calc = NULL) {
  dt <- ncdf4::ncdim_def("timeDim", "", seq_along(tt), create_dimvar = FALSE)
  dh <- ncdf4::ncdim_def("Mean_Layer_HeightDim", "", 1L, create_dimvar = FALSE)
  dq <- ncdf4::ncdim_def("Mean_Layer_QualityIndexDim", "", 1L, create_dimvar = FALSE)
  v <- list(t = ncdf4::ncvar_def("time", "", list(dt), prec = "integer"),
            h = ncdf4::ncvar_def("Mean_Layer_Height", "", list(dh, dt), prec = "integer"),
            q = ncdf4::ncvar_def("Mean_Layer_QualityIndex", "", list(dq, dt), prec = "integer"))
  if (v312) {
    du <- ncdf4::ncdim_def("location_utc_offsetDim", "", 1L, create_dimvar = FALSE)
    dc <- ncdf4::ncdim_def("Mean_Layer_Calculation_TimeDim", "", 1L, create_dimvar = FALSE)
    v$u <- ncdf4::ncvar_def("location_utc_offset", "", list(du, dt), prec = "float")
    v$c <- ncdf4::ncvar_def("Mean_Layer_Calculation_Time", "", list(dc, dt), prec = "integer")
  }
  nc <- ncdf4::nc_create(path, unname(v))
  ncdf4::ncvar_put(nc, v$t, tt)
  ncdf4::ncatt_put(nc, "time", "Units", "seconds since 1970-01-01 00:00:00.000")
  ncdf4::ncvar_put(nc, v$h, h)
  ncdf4::ncvar_put(nc, v$q, q)
  if (v312) {
    ncdf4::ncvar_put(nc, v$u, rep(0, length(tt)))
    ncdf4::ncvar_put(nc, v$c, calc)
  }
  ncdf4::nc_close(nc)
  path
}

test_that("BL-View 3.1.1 / 3.1.2 files: version from the variables, hourly collapse, -999 to NA", {
  d0 <- as.numeric(as.POSIXct("2025-03-23", tz = "UTC"))
  tt <- d0 + 55 + 60 * (0:179)
  h <- c(rep(600L, 60), rep(-999L, 60), rep(1200L, 60))
  q <- c(rep(3L, 60), rep(-999L, 60), rep(1L, 60))
  calc <- rep(tt[c(1, 61, 121)], each = 60)
  f1 <- write_blview31_nc(tempfile(fileext = ".nc"), tt, h, q)
  f2 <- write_blview31_nc(tempfile(fileext = ".nc"), tt, h, q, v312 = TRUE, calc = calc)
  x1 <- read_blview_l3(f1, 0)
  x2 <- read_blview_l3(f2, 0)
  expect_equal(unique(x1$version), "3.1.1")
  expect_equal(unique(x2$version), "3.1.2")
  expect_equal(x1$mlh_m, c(600, NA, 1200))
  expect_equal(x1$quality, c(3L, NA, 1L))
  expect_equal(as.numeric(x1$hour_end_utc), d0 + 3600 * (0:2))
  expect_equal(as.numeric(x2$calc_time_utc) - as.numeric(x2$hour_end_utc), rep(55, 3))
  expect_equal(attr(x2, "blview_meta")$location_utc_offset, 0)
  # a file without the MLH variables
  dt <- ncdf4::ncdim_def("timeDim", "", 1:3, create_dimvar = FALSE)
  nc <- ncdf4::nc_create(f3 <- tempfile(fileext = ".nc"), list(ncdf4::ncvar_def("time", "", list(dt), prec = "integer")))
  ncdf4::ncvar_put(nc, "time", tt[1:3])
  ncdf4::nc_close(nc)
  x3 <- read_blview_l3(f3, 0)
  expect_equal(nrow(x3), 0)
  expect_equal(attr(x3, "blview_meta")$version, "no_mlh")
})

test_that("BL-View directory: version table, overlaps keep the newest version, standard table", {
  d <- withr::local_tempdir()
  d0 <- as.numeric(as.POSIXct("2025-10-13", tz = "UTC"))
  tt <- d0 + 6 + 60 * (0:119)
  write_blview31_nc(file.path(d, "a_311.nc"), tt, c(rep(500L, 60), rep(700L, 60)), c(rep(2L, 60), rep(3L, 60)))
  write_blview_nc(file.path(d, "b_322.nc"), tt[61:120] + 3600, c(rep(710L, 60)), rep(3L, 60))
  write_blview_nc(file.path(d, "c_322.nc"), tt[61:120], c(rep(720L, 60)), rep(3L, 60)) # same hour as a_311 hour 2
  r <- read_blview_dir(d, 0)
  expect_equal(sort(r$files$version), c("3.1.1", "3.2.2", "3.2.2"))
  vt <- blview_version_table(r)
  expect_equal(vt$version, c("3.1.1", "3.2.2"))
  expect_equal(vt$n_hours_valid, c(2L, 2L))
  expect_equal(vt$n_q3, c(1L, 2L))
  dd <- blview_dedup_hours(r$hours)
  expect_equal(nrow(dd$hours), 3)
  expect_equal(dd$hours$mlh_m[dd$hours$hour_end_utc == as.POSIXct(d0 + 3600, origin = "1970-01-01", tz = "UTC")], 720)
  expect_equal(dd$overlap$n_hours, 1L)
  expect_equal(dd$overlap$n_differ, 1L)
  s <- suppressMessages(standardize_blview_mlh(dd$hours, "XX-Tst", read_exclusions(fixture_path("exclusion_periods.csv")),
                                               1, time_offset_h = 0))
  expect_equal(s$data$pblh_m, c(500, 720, 710))
  expect_equal(as.numeric(s$data$timestamp_end_utc - s$data$timestamp_start_utc, units = "secs"), rep(3600, 3))
  expect_equal(s$data$x_version, c("3.1.1", "3.2.2", "3.2.2"))
  expect_true(all(is.na(s$data$qc_flag)))
  expect_no_error(validate_standard_pblh(s$data, s$meta))
  expect_error(standardize_blview_mlh(r$hours, "XX-Tst", read_exclusions(fixture_path("exclusion_periods.csv")), 1, 0),
               "duplicated")
})

test_that("the real BL-View samples: one file per version family is read", {
  dd <- here::here("data", "raw", "measured", "CZ-Lnz", "samples")
  skip_if_not(dir.exists(dd), "BL-View samples not on this machine")
  r <- read_blview_dir(dd, 0)
  expect_setequal(unique(r$files$version), c("3.1.1", "3.1.2", "3.2.2"))
  expect_true(all(is.na(r$files$error)))
  expect_equal(sum(r$files$n_inconsistent), 0L)
  expect_true(all(r$files$date_stamp_matches_utc[r$files$version != "3.2.2"]))
})

# --- Part C diagnostics: hour bins, morning bins, lag test ---------------------------------------------------------

test_that("hour_bins and instant_hour_bins put values on UTC hours", {
  t <- as.POSIXct("2025-07-01 10:00", tz = "UTC") + 1800 * 0:3
  b <- hour_bins(t, t + 1800, c(100, 200, NA, 400))
  expect_equal(format(b$bin_utc, "%H"), c("10", "11"))
  expect_equal(b$value, c(150, 400))
  expect_equal(b$n, c(2L, 1L))
  e <- instant_hour_bins(as.POSIXct("2025-07-01 10:00", tz = "UTC") + 3600 * 0:2, c(100, 300, NA))
  expect_equal(nrow(e), 1)
  expect_equal(e$value, 200)
})

test_that("lag_correlation finds the lag of a series stamped one hour late (CET)", {
  lat <- 50.96
  lon <- 13.57
  days <- as.Date("2025-06-01") + 0:19
  mb <- morning_bins(days, lat, lon)
  expect_true(all(mb$season == "JJA"))
  st <- sun_times(days[1], lat, lon)
  expect_true(all(as.numeric(mb$bin_utc[mb$day == days[1]]) + 1800 >= as.numeric(st$sunrise_utc)))
  # reference: a growth curve with day-to-day amplitude; measured = same curve stamped 1 h later
  hrs <- as.POSIXct(paste(rep(days, each = 24), sprintf("%02d:00", 0:23)), tz = "UTC")
  amp <- rep(seq(800, 1800, length.out = 20), each = 24)
  h <- as.numeric(format(hrs, "%H"))
  v <- pmax(100, amp * sin(pmin(pmax((h - 3) / 10, 0), 1) * pi / 2))
  ref <- dplyr::tibble(bin_utc = hrs, value = v, n = 1L)
  meas <- dplyr::tibble(bin_utc = hrs + 3600, value = v, n = 1L)
  lc <- lag_correlation(meas, ref, mb)
  b <- best_lag(lc)
  expect_equal(b$best_lag_h[b$season == "all"], 1)
  expect_equal(b$r_best[b$season == "all"], 1, tolerance = 1e-12)
  expect_true(all(is.na(lc$r[lc$season == "DJF"])))
  b0 <- best_lag(lag_correlation(ref, ref, mb))
  expect_equal(b0$best_lag_h[b0$season == "all"], 0)
})

test_that("the secondary sources have style rows", {
  expect_no_error(style_of(c("tud_skyvue_mlh", "vaisala_blview_mlh")))
  expect_no_error(style_of(c("skyvue_layer1", "skyvue_layer2", "skyvue_layer3", "blview_all", "blview_3.1.1",
                             "blview_3.1.2", "blview_3.2.2", "blview_q1", "blview_q2", "blview_q3"), group = "secondary"))
})

# --- standard files: explicit column types ---------------------------------------------------------------------------

test_that("read_standard_pblh: an x_ column empty in its first rows is read as numeric, not logical", {
  n <- 2500
  t <- as.POSIXct("2024-01-01", tz = "UTC") + 1800 * (seq_len(n) - 1)
  df <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t + 1800, pblh_m = 500, source_id = "src_x",
                      x_dlh_m = c(rep(NA, 2000), rep(1234.5, n - 2000)), x_regime = c(rep(NA, 2000), rep("stable", 500)),
                      x_consistent = c(rep(NA, 1500), rep(TRUE, 1000)), x_note = c(rep(NA, 2400), rep("a", 100)),
                      x_empty = NA_real_,
                      x_first_profile_utc = c(rep(NA, 2000), rep("2024-01-01T00:01:14Z", 500)))
  meta <- standard_meta("XX-Tst", "src_x", "aerosol", time_reference = "interval_mean")
  dir <- withr::local_tempdir()
  p <- write_standard_pblh(df, meta, dir)
  x <- read_standard_pblh(p[["data"]])
  expect_type(x$x_dlh_m, "double")
  expect_equal(sum(x$x_dlh_m == 1234.5, na.rm = TRUE), 500)
  expect_type(x$x_regime, "character")
  expect_type(x$x_consistent, "logical")
  expect_equal(sum(x$x_consistent, na.rm = TRUE), 1000)
  expect_type(x$x_note, "character") # unknown, not numeric: character
  expect_type(x$x_empty, "double")   # unknown, empty: double
  expect_s3_class(x$x_first_profile_utc, "POSIXct")
  expect_type(x$qc_flag, "integer")
})

# --- fix-14c ---------------------------------------------------------------------------------------------------------

# Synthetic morning growth on hourly bins: day-to-day amplitude, growth from 04 to 12 UTC.
growth_bins <- function(days, shift_h = 0) {
  hrs <- as.POSIXct(paste(rep(days, each = 24), sprintf("%02d:00", 0:23)), tz = "UTC")
  amp <- rep(seq(800, 1800, length.out = length(days)), each = 24)
  h <- as.numeric(format(hrs, "%H", tz = "UTC"))
  v <- pmax(100, amp * sin(pmin(pmax((h - 4) / 8, 0), 1) * pi / 2))
  dplyr::tibble(bin_utc = hrs + shift_h * 3600, value = v, n = 1L)
}

test_that("lag_table keeps the reference names and best_lag gives one row per series x reference x season", {
  days <- as.Date("2025-06-01") + 0:19
  mb <- morning_bins(days, 50.96, 13.57)
  ref <- growth_bins(days)
  lt <- lag_table(list(s1 = growth_bins(days, 1), s2 = growth_bins(days, 0)),
                  list(pblhmod_default = ref, era5_blh_center = ref), mb, "XX-Tst")
  expect_setequal(unique(lt$reference), c("pblhmod_default", "era5_blh_center")) # was: correlations (data masking)
  expect_equal(range(lt$lag_h), c(-3, 6))
  b <- best_lag(lt, by = c("site_id", "series", "reference"))
  expect_equal(nrow(b), nrow(unique(lt[!is.na(lt$r), c("site_id", "series", "reference", "season")])))
  expect_equal(anyDuplicated(b[c("series", "reference", "season")]), 0L)
  expect_equal(b$best_lag_h[b$series == "s1" & b$season == "all"], c(1L, 1L))
  expect_equal(b$best_lag_h[b$series == "s2" & b$season == "all"], c(0L, 0L))
  expect_equal(b$r_lag0[b$series == "s2" & b$season == "all"], c(1, 1), tolerance = 1e-12)
  expect_equal(reference_label(c("pblhmod_default", "era5_blh_center", "x")), c("PBLHmod", "ERA5", "x"))
  # a hand-made lag table: the maximum, r at lag 0, and a group whose lag 0 has no correlation
  lc <- dplyr::tibble(series = "s", reference = rep(c("a", "b"), each = 4), season = "all", lag_h = rep(c(-1, 0, 1, 2), 2),
                      n = 20L, r = c(0.1, 0.5, 0.9, 0.3, 0.2, NA, 0.4, 0.8))
  bb <- best_lag(lc, by = c("series", "reference"))
  expect_equal(bb$best_lag_h, c(1L, 2L))
  expect_equal(bb$r_best, c(0.9, 0.8))
  expect_equal(bb$r_lag0, c(0.5, NA))
  expect_equal(nrow(best_lag(lc[0, ], by = c("series", "reference"))), 0)
})

test_that("local-time stamps across a summer-time change: lag +1 / +2 before, 0 after the conversion", {
  # true hourly values in UTC; stored stamps = local clock time of Europe/Prague (CET +1, CEST +2), read as UTC
  days <- as.Date("2025-03-01") + 0:59 # spans the change of 2025-03-30
  truth <- growth_bins(days)
  stored <- truth
  stored$bin_utc <- as.POSIXct(format(truth$bin_utc, "%Y-%m-%d %H:%M:%S", tz = "Europe/Prague"), tz = "UTC")
  off <- (as.numeric(stored$bin_utc) - as.numeric(truth$bin_utc)) / 3600
  expect_equal(sort(unique(off)), c(1, 2))
  mb <- morning_bins(days, 50.96, 13.57)
  before <- best_lag(lag_correlation(stored, truth, mb[mb$day < as.Date("2025-03-30"), ]))
  expect_equal(before$best_lag_h[before$season == "all"], 1L) # winter part: CET
  after_st <- best_lag(lag_correlation(stored, truth, mb[mb$day > as.Date("2025-03-31"), ]))
  expect_equal(after_st$best_lag_h[after_st$season == "all"], 2L) # summer part: CEST
  # conversion of hour ends (stamps of the standard table) back to UTC
  back <- blview_local_to_utc(stored$bin_utc, "Europe/Prague")
  expect_equal(as.numeric(back), as.numeric(truth$bin_utc))
  conv <- truth
  conv$bin_utc <- back
  b <- best_lag(lag_correlation(conv, truth, mb))
  expect_equal(b$best_lag_h, rep(0L, nrow(b)))
  # the spring-forward clock hour does not exist: NA; standard rows of one version only
  expect_true(is.na(blview_local_to_utc(as.POSIXct("2025-03-30 02:30", tz = "UTC"), "Europe/Prague")))
  t <- as.POSIXct("2025-03-30 01:00", tz = "UTC") + 3600 * 0:3 # stored 01, 02 (nonexistent), 03, 04 local
  d <- dplyr::tibble(timestamp_start_utc = t - 3600, timestamp_end_utc = t, pblh_m = 1:4, x_version = rep("3.2.2", 4))
  d <- dplyr::bind_rows(d, dplyr::tibble(timestamp_start_utc = t[1] - 4 * 3600, timestamp_end_utc = t[1] - 3 * 3600, pblh_m = 0,
                                         x_version = "3.1.2"))
  cv <- blview_convert_version(d, "3.2.2", "Europe/Prague")
  expect_equal(cv$n_dropped, 1L)
  expect_equal(format(cv$data$timestamp_end_utc[cv$data$x_version == "3.2.2"], "%H:%M", tz = "UTC"),
               c("00:00", "01:00", "02:00"))
  expect_equal(cv$data$timestamp_end_utc[cv$data$x_version == "3.1.2"], t[1] - 3 * 3600) # other versions unchanged
  expect_equal(cv$n_duplicated, 0L)
})

test_that("the +1 h baseline: the BL-View chain matched to half-hourly references gives lag 0 and ratio 1", {
  # truth: half-hourly PBLH; BL-View: the mean of the hour ending at hh:00 carried by the profiles hh:00-hh:59
  days <- as.Date("2025-06-01") + 0:13
  hh <- as.POSIXct(paste(rep(days, each = 48), sprintf("%02d:%02d", rep(0:23, each = 2), c(0, 30))), tz = "UTC")
  hr <- as.numeric(format(hh, "%H", tz = "UTC")) + as.numeric(format(hh, "%M", tz = "UTC")) / 60 + 0.25
  amp <- rep(seq(800, 1800, length.out = length(days)), each = 48)
  truth <- dplyr::tibble(timestamp_start_utc = hh, timestamp_end_utc = hh + 1800,
                         pblh_m = pmax(100, amp * sin(pmin(pmax((hr - 4) / 8, 0), 1) * pi / 2)))
  hour_mean <- tapply(truth$pblh_m, floor(as.numeric(hh) / 3600) * 3600, mean) # [h, h + 1 h)
  h_start <- as.numeric(names(hour_mean))
  # profiles at hh:mm:50 of hour hh carry the mean of [hh - 1 h, hh)
  tt <- rep(h_start + 3600, each = 60) + 50 + 60 * (0:59)
  vals <- rep(as.integer(round(hour_mean)), each = 60)
  f <- write_blview_nc(tempfile(fileext = ".nc"), tt, vals, rep(3L, length(tt)))
  x <- read_blview_l3(f, 0)
  s <- suppressMessages(standardize_blview_mlh(x, "XX-Tst", read_exclusions(fixture_path("exclusion_periods.csv")),
                                               1, 0))$data
  expect_equal(as.numeric(s$timestamp_start_utc), h_start) # labelled by the hour it averages
  mb <- morning_bins(days, 50.96, 13.57)
  lc <- lag_correlation(hour_bins(s$timestamp_start_utc, s$timestamp_end_utc, s$pblh_m),
                        hour_bins(truth$timestamp_start_utc, truth$timestamp_end_utc, truth$pblh_m), mb)
  b <- best_lag(lc)
  expect_equal(b$best_lag_h[b$season == "all"], 0L)
  # window-metric matching (measured_time_base, hourly base): PBLHmod mean of the two half-hours = the hourly value
  mod <- dplyr::mutate(truth, x_regime = "convective")
  era_t <- as.POSIXct(h_start, origin = "1970-01-01", tz = "UTC")
  era <- list(era5_blh_center = dplyr::tibble(timestamp_start_utc = era_t, timestamp_end_utc = era_t, pblh_m = 500))
  tb <- measured_time_base(s, mod, era, step_s = 3600)
  pt <- pair_table(tb$mod, tb$meas[tb$meas$kind %in% "exact", c("timestamp_utc", "pblh_m")], "x")
  expect_gt(nrow(pt), 100)
  expect_lte(max(abs(pt$mod - pt$ref)), 0.5) # rounding of the integer MLH only
})

test_that("window metrics per SkyVUE layer although pblh_m is NA", {
  t <- as.POSIXct("2025-07-01 08:00", tz = "UTC") + 1800 * 0:11 # 08:00-14:00 UTC, in the window at DE-Tha
  meas <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t + 1800, pblh_m = NA_real_, qc_flag = NA_integer_,
                        x_layer1_m = 400, x_layer2_m = 1200)
  mod <- dplyr::tibble(timestamp_start_utc = t, timestamp_end_utc = t + 1800, pblh_m = 800, x_regime = "convective")
  era_t <- seq(min(t), max(t) + 1800, by = 3600)
  era <- list(era5_blh_center = dplyr::tibble(timestamp_start_utc = era_t, timestamp_end_utc = era_t, pblh_m = 1000))
  opts <- measured_compare_options()
  expect_null(secondary_window_metrics(meas, "x", mod, era, 13.57, opts, "XX-Tst")) # pblh_m NA: nothing
  w <- dplyr::bind_rows(lapply(1:2, function(k) {
    secondary_window_metrics(dplyr::mutate(meas, pblh_m = .data[[paste0("x_layer", k, "_m")]]),
                             paste0("skyvue_layer", k), mod, era, 13.57, opts, "XX-Tst")
  }))
  expect_equal(unique(w$quality), "all")
  p <- w[w$source == "pblhmod_default", ]
  expect_equal(p$median_log_ratio, c(log(800 / 400), log(800 / 1200)))
  expect_true(all(p$n > 0))
})

# --- fix-14a item 11: BL-View time correction ---------------------------------------------------------------------

test_that("blview_time_correction: 3.2.2 Prague local time to UTC, then -1 h for every version; DST hours counted", {
  u <- function(x) as.POSIXct(x, tz = "UTC")
  hours <- dplyr::tibble(
    hour_end_utc = u(c("2025-01-15 12:00", "2025-01-15 12:00", "2025-07-15 12:00", "2025-03-30 01:00",
                       "2025-03-30 02:00", "2025-03-30 03:00")),
    version = c("3.1.1", "3.2.2", "3.2.2", "3.2.2", "3.2.2", "3.2.2"), mlh_m = c(500, 510, 900, 300, 310, 320),
    file = "f.nc")
  tc <- blview_time_correction(hours, "Europe/Prague", -1)
  h <- tc$hours
  expect_equal(tc$n_dropped_dst, 1L) # 02:00 does not exist in Prague on 2025-03-30
  expect_equal(tc$dropped_dst, "2025-03-30 02:00")
  expect_equal(nrow(h), 5)
  key <- function(m) h[h$mlh_m == m, ]
  expect_equal(key(500)$hour_end_utc, u("2025-01-15 11:00")) # 3.1.1: -1 h only
  expect_equal(key(510)$hour_end_utc, u("2025-01-15 10:00")) # CET (UTC+1), then -1 h
  expect_equal(key(900)$hour_end_utc, u("2025-07-15 09:00")) # CEST (UTC+2), then -1 h
  expect_equal(key(300)$hour_end_utc, u("2025-03-29 23:00")) # 01:00 CET
  expect_equal(key(320)$hour_end_utc, u("2025-03-30 00:00")) # 03:00 CEST
  expect_equal(sort(unique(h$time_shift_h)), c(-3, -2, -1))
  expect_equal(key(510)$time_shift_h, -2)
  expect_equal(key(900)$time_shift_h, -3)
  # no correction: stored time kept
  t0 <- blview_time_correction(hours, NULL, 0)
  expect_equal(sort(t0$hours$hour_end_utc), sort(hours$hour_end_utc))
  expect_true(all(t0$hours$time_shift_h == 0))
  expect_error(blview_time_correction(hours, "Mars/Olympus", -1), "time zone")
})

test_that("blview_time_correction: autumn repeated hour does not create rows, overlapping files are not DST collisions", {
  u <- function(x) as.POSIXct(x, tz = "UTC")
  st <- u("2025-10-26 00:00") + 3600 * (0:5) # stored local clock 00:00 ... 05:00 across the repeated 02:00-03:00
  hours <- dplyr::tibble(hour_end_utc = c(st, st[1]), version = "3.2.2", mlh_m = c(1:6, 1), file = c(rep("a", 6), "b"))
  tc <- blview_time_correction(hours, "Europe/Prague", -1)
  expect_equal(tc$n_dropped_dst, 0L)
  expect_equal(nrow(tc$hours), 7)
  # the overlap a / b on the same stored hour is not counted as a DST duplicate (blview_dedup_hours() handles it)
  expect_lte(tc$n_duplicated_dst, 1L)
  expect_false(any(tc$duplicated_dst == format(u("2025-10-25 21:00"), "%Y-%m-%d %H:%M")))
  # every converted hour is in UTC and strictly earlier than the stored one
  expect_true(all(tc$hours$time_shift_h %in% c(-2, -3)))
})

test_that("standardize_blview_mlh records the time correction (x_time_shift_h, sidecar history and limitation)", {
  d <- withr::local_tempdir()
  tt <- as.numeric(as.POSIXct("2025-07-15 10:00:30", tz = "UTC")) + 60 * (0:119)
  write_blview_nc(file.path(d, "b_322.nc"), tt, c(rep(700L, 60), rep(800L, 60)), rep(3L, 120))
  r <- read_blview_dir(d, 0)
  tc <- blview_time_correction(r$hours, "Europe/Prague", -1)
  dd <- blview_dedup_hours(tc$hours)
  ex <- read_exclusions(fixture_path("exclusion_periods.csv"))
  s <- suppressMessages(standardize_blview_mlh(dd$hours, "XX-Tst", ex, 1, 0, time_correction = tc[names(tc) != "hours"]))
  expect_equal(s$data$timestamp_end_utc, as.POSIXct(c("2025-07-15 07:00", "2025-07-15 08:00"), tz = "UTC"))
  expect_equal(s$data$x_time_shift_h, c(-3, -3))
  expect_match(s$meta$processing_history, "Europe/Prague local time")
  expect_match(s$meta$processing_history, "-1 h, an empirical alignment")
  expect_match(s$meta$processing_history, "SON +3 h on 68 days", fixed = TRUE)
  expect_match(s$meta$known_limitations, "corrected empirically")
  s0 <- suppressMessages(standardize_blview_mlh(r$hours, "XX-Tst", ex, 1, 0))
  expect_false("x_time_shift_h" %in% names(s0$data))
  expect_match(s0$meta$processing_history, "NOT corrected")
  expect_error(standardize_blview_mlh(r$hours, "XX-Tst", ex, 1, 0, time_correction = list(shift_h = -1)), "not corrected")
})
