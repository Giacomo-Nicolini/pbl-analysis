api_sample <- function() jsonlite::fromJSON(fixture_path("cloudnet/api_files_sample.json"))

test_that("distance_class applies the distance rule at its boundaries", {
  d <- c(0, 4.99, 5, 25, 25.01, 50, 50.01, NA)
  expect_equal(distance_class(d), c("all_regimes", "all_regimes", "daytime_convective", "daytime_convective",
                                    "daytime_climatology", "daytime_climatology", "not_used", NA))
})

test_that("cloudnet_flatten_files returns one typed row per file", {
  f <- cloudnet_flatten_files(api_sample())
  expect_equal(nrow(f), 4)
  expect_equal(f$site, rep("testsite", 4))
  expect_equal(f$product, rep("epsilon-lidar", 4))
  expect_s3_class(f$measurement_date, "Date")
  expect_equal(f$size_bytes, c(1e6, 2e6, 5e5, 3e5))
  expect_equal(f$instrument_name, c("TEST HALO 1", "TEST HALO 1", "TEST HALO 1", "TEST HALO 2"))
  expect_equal(f$instrument_pid[4], "https://hdl.handle.net/x/halo2")
  expect_equal(f$volatile, c(FALSE, TRUE, FALSE, FALSE))
  expect_equal(f$pid[2], "") # volatile files have no PID
  expect_equal(nchar(f$checksum_sha256), rep(64L, 4))
})

test_that("cloudnet_flatten_files handles an empty response", {
  f <- cloudnet_flatten_files(list())
  expect_equal(nrow(f), 0)
  expect_true(all(c("uuid", "size_bytes", "instrument_pid") %in% names(f)))
  expect_equal(nrow(cloudnet_inventory_summary(f)), 0)
  expect_equal(nrow(cloudnet_coverage_periods(cloudnet_inventory_summary(f))), 0)
})

test_that("cloudnet_inventory_summary aggregates by instrument and month", {
  inv <- cloudnet_inventory_summary(cloudnet_flatten_files(api_sample()))
  expect_equal(nrow(inv), 3)
  h1 <- inv[inv$instrument_name == "TEST HALO 1" & inv$month == "2022-01", ]
  expect_equal(h1$n_files, 2)
  expect_equal(h1$size_mb, 3)
  expect_equal(h1$n_volatile, 1)
  expect_equal(h1$mean_coverage, 0.375)
  expect_equal(inv$n_legacy[inv$month == "2022-03"], 1)
})

test_that("cloudnet_coverage_periods splits runs at missing months", {
  inv <- cloudnet_inventory_summary(cloudnet_flatten_files(api_sample()))
  p <- cloudnet_coverage_periods(inv)
  expect_equal(nrow(p), 2) # 2022-01 and 2022-03 (2022-02 missing)
  expect_equal(p$first_month, c("2022-01", "2022-03"))
  expect_equal(p$n_files, c(3, 1)) # both instruments counted in 2022-01
  # consecutive months across a year boundary form one run
  inv2 <- dplyr::tibble(site = "s", product = "p", instrument_name = "i", instrument_pid = "x",
                        month = c("2021-11", "2021-12", "2022-01", "2022-05"), n_files = 1L, size_mb = 1,
                        n_legacy = 0L, n_volatile = 0L, mean_coverage = 1)
  p2 <- cloudnet_coverage_periods(inv2)
  expect_equal(p2$n_months, c(3, 1))
  expect_equal(p2$last_month, c("2022-01", "2022-05"))
})

test_that("read_measured_sites validates the pairing table", {
  ms <- read_measured_sites(here::here("config", "measured_sites.csv"))
  expect_true(all(c("FI-Hyy", "FI-Lom", "DE-RuS", "SE-Nor", "FR-Gri") %in% ms$icos_site))
  expect_true(all(ms$distance_km >= 0))
  # distance_km consistent with the coordinates in sites.csv
  sites <- read_sites(here::here("config", "sites.csv"))
  m <- merge(ms, sites, by.x = "icos_site", by.y = "site_id")
  expect_equal(m$distance_km, round(haversine_km(m$lat, m$lon, m$provider_lat, m$provider_lon), 2), tolerance = 1e-9)

  tmp <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(ms[c(1, 1), ], tmp)
  expect_error(read_measured_sites(tmp), "duplicated")
  readr::write_csv(ms[, -2], tmp)
  expect_error(read_measured_sites(tmp), "missing[[:space:]]+column")
})

test_that("download plan skips files present with the API checksum", {
  tmp <- withr::local_tempdir()
  writeBin(charToRaw("abc"), file.path(tmp, "a.nc"))
  writeBin(charToRaw("abd"), file.path(tmp, "b.nc"))
  sha_abc <- "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad" # SHA-256("abc")
  expect_equal(cloudnet_sha256(file.path(tmp, c("a.nc", "zz.nc"))), c(sha_abc, NA))
  files <- dplyr::tibble(filename = c("a.nc", "b.nc", "c.nc"), checksum_sha256 = sha_abc)
  p <- cloudnet_download_plan(files, tmp)
  expect_equal(p$present, c(TRUE, FALSE, FALSE)) # b.nc has another checksum, c.nc absent
  expect_equal(p$local_path[3], file.path(tmp, "c.nc"))
})

test_that("cloudnet_download verifies checksums and logs failures", {
  src <- withr::local_tempdir()
  dst <- withr::local_tempdir()
  writeBin(charToRaw("abc"), file.path(src, "good.nc"))
  writeBin(charToRaw("xyz"), file.path(src, "bad.nc"))
  sha_abc <- "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
  furl <- function(f) paste0("file:///", normalizePath(file.path(src, f), winslash = "/"))
  files <- dplyr::tibble(
    filename = c("good.nc", "bad.nc", "missing.nc"), checksum_sha256 = sha_abc,
    download_url = c(furl("good.nc"), furl("bad.nc"), paste0(furl("good.nc"), "_nope"))
  )
  p <- suppressMessages(cloudnet_download(cloudnet_download_plan(files, dst), retries = 1, pause_s = 0))
  expect_equal(p$status, c("downloaded", "failed", "failed"))
  expect_equal(p$message[2], "checksum mismatch")
  expect_true(file.exists(file.path(dst, "good.nc")))
  expect_false(any(file.exists(file.path(dst, c("bad.nc", "bad.nc.part", "missing.nc")))))
  # second run: nothing to download
  p2 <- suppressMessages(cloudnet_download(cloudnet_download_plan(files[1, ], dst), pause_s = 0))
  expect_equal(p2$status, "present")
})

# ---- MLH (spec 07 section 3) ----
f_eps <- function() fixture_path("cloudnet/20220101_testsite_epsilon-lidar_aaaa0001.nc")
f_cls <- function() fixture_path("cloudnet/20220101_testsite_classification.nc")

test_that("read_cloudnet_epsilon converts heights to agl with the file altitude", {
  e <- read_cloudnet_epsilon(f_eps())
  expect_equal(e$altitude_m, 150)
  expect_equal(e$height_asl, 165 + 30 * (0:5))
  expect_equal(e$height_agl, 15 + 30 * (0:5)) # height_asl - 150
  expect_equal(e$gate_m, 30)
  expect_equal(dim(e$eps), c(6L, 4L))
  expect_true(all(is.na(e$eps[1:3, ]))) # fill value -> NA
  expect_equal(as.numeric(e$time - as.POSIXct("2022-01-01", tz = "UTC"), units = "secs"),
               c(15, 30, 1200, 2400), tolerance = 1e-3)
  expect_equal(e$attrs$instrument_pid, "https://hdl.handle.net/x/halo1")
  expect_equal(e$attrs$rolling_window_s, 180)
  # other ground altitude (sensitivity): 181 m
  e2 <- read_cloudnet_epsilon(f_eps(), ground_altitude_m = 181)
  expect_equal(e2$height_agl, 165 + 30 * (0:5) - 181)
})

test_that("cloudnet_mlh_profiles: MLH at the top of the last turbulent gate, censoring and no signal", {
  e <- read_cloudnet_epsilon(f_eps())
  opts <- cloudnet_mlh_options(z_min_agl = 105)
  p <- cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, opts)
  expect_equal(p$status, c("measured", "no_signal", "lowest_gate_missing", "censored_weak"))
  expect_equal(p$mlh_m, c(120, NA, NA, NA)) # gate centre 105 m + 15 m
  expect_equal(p$censored, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(p$z_min_m, rep(105, 4))
  # default z_min: lowest gate with any valid epsilon in the file (105 m here)
  expect_equal(cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, cloudnet_mlh_options())$status, p$status)
  expect_error(cloudnet_mlh_profiles(e$eps, e$height_agl, e$gate_m, cloudnet_mlh_options(z_min_agl = 400)),
               "half a gate")
})

test_that("cloudnet_mlh_profiles: gap in the column and full-column cases", {
  z <- 15 + 30 * (0:5)
  eps <- cbind(
    c(NA, NA, NA, 1e-3, 1e-3, NA),     # stops at a missing gate: gap, MLH 150
    c(NA, NA, NA, 1e-3, 1e-3, 1e-3),   # turbulent up to the top gate: MLH 180
    c(NA, NA, NA, 1e-4, 9.9e-5, 1e-3), # threshold inclusive at the lowest gate: MLH 120
    c(NA, NA, NA, 1e-3, 1e-5, 1e-3)    # weak gate interrupts (surface-connected only): MLH 120
  )
  p <- cloudnet_mlh_profiles(eps, z, 30, cloudnet_mlh_options(z_min_agl = 105))
  expect_equal(p$status, rep("measured", 4))
  expect_equal(p$mlh_m, c(150, 180, 120, 120))
  expect_equal(p$gap_in_column, c(TRUE, FALSE, FALSE, FALSE))
  # gates below z_min are ignored
  eps2 <- cbind(c(1e-5, 1e-3, 1e-3, 1e-3, 1e-5, NA))
  expect_equal(cloudnet_mlh_profiles(eps2, z, 30, cloudnet_mlh_options(z_min_agl = 45))$mlh_m, 120)
  expect_equal(cloudnet_mlh_profiles(eps2, z, 30, cloudnet_mlh_options(z_min_agl = 15))$status, "censored_weak")
})

test_that("read_cloudnet_classification and the precipitation flag", {
  cls <- read_cloudnet_classification(f_cls())
  expect_equal(cls$height_agl, c(102, 150, 450))
  expect_equal(cls$cloud_base_agl_m, c(1000, 1000, NA, NA, NA))
  expect_equal(cls$rain_detected, c(FALSE, FALSE, FALSE, TRUE, FALSE))
  e <- read_cloudnet_epsilon(f_eps())
  pr <- cloudnet_precip_flag(e$time, cls, cloudnet_mlh_options())
  # profiles at 15 and 30 s: class 2 at 150 m agl at 60 s (within 90 s); 1200 and 2400 s: no classification
  expect_equal(pr$precip, c(TRUE, TRUE, NA, NA))
  expect_equal(pr$rain_detected, c(TRUE, TRUE, NA, NA))
  expect_equal(pr$cloud_base_agl_m, c(1000, 1000, NA, NA))
  # only gates up to 120 m agl checked: no precipitation (class 2 at 150 and 450 m agl ignored)
  expect_equal(cloudnet_precip_flag(e$time, cls, cloudnet_mlh_options(precip_top_agl = 120))$precip,
               c(FALSE, FALSE, NA, NA))
  # window of 10 s: class-2 profile at 60 s outside both windows; no classification profile within 10 s of 15 s
  expect_equal(cloudnet_precip_flag(e$time, cls, cloudnet_mlh_options(precip_window_s = 10))$precip,
               c(NA, FALSE, NA, NA))
  # no classification file: not checked
  expect_true(all(is.na(cloudnet_precip_flag(e$time, NULL, cloudnet_mlh_options())$precip)))
})

test_that("cloudnet_mlh_day: precipitation overrides the epsilon status", {
  d <- cloudnet_mlh_day(f_eps(), f_cls(), cloudnet_mlh_options(z_min_agl = 105))
  expect_equal(d$status, c("precipitation", "precipitation", "lowest_gate_missing", "censored_weak"))
  expect_true(all(is.na(d$mlh_m)))
  expect_equal(d$censored, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(d$n_expected, rep(2, 4)) # 3 profiles in the first half-hour, 1 in the second
  d0 <- cloudnet_mlh_day(f_eps(), NULL, cloudnet_mlh_options(z_min_agl = 105))
  expect_equal(d0$status[1:2], c("measured", "no_signal"))
  expect_true(all(is.na(d0$precip)))
})

test_that("cloudnet_expected_profiles = median count per half-hour", {
  t0 <- as.POSIXct("2022-01-01", tz = "UTC")
  t <- t0 + c(seq(0, 600, by = 15), 1800 + seq(0, 300, by = 15), 3600 + seq(0, 900, by = 15))
  expect_equal(cloudnet_expected_profiles(t), 41) # counts 41, 21, 61
  expect_true(is.na(cloudnet_expected_profiles(t0[0])))
})

# Synthetic profile table for the aggregation rules: n profiles in one half-hour.
mk_prof <- function(status, mlh = NA_real_, start = "2022-06-01 12:00:00", n_expected = 10) {
  n <- length(status)
  dplyr::tibble(
    time = as.POSIXct(start, tz = "UTC") + 14 * seq_len(n), status = status,
    mlh_m = ifelse(status == "measured", rep_len(mlh, n), NA_real_),
    censored = status %in% c("censored_weak", "lowest_gate_missing"), gap_in_column = FALSE, z_min_m = 105,
    precip = status == "precipitation", rain_detected = NA, cloud_base_agl_m = NA_real_,
    n_expected = n_expected, dt_s = 14, file_uuid = "u"
  )
}

test_that("cloudnet_mlh_halfhour: coverage, censored majority and partly censored", {
  o <- cloudnet_mlh_options(min_coverage = 0.5)
  # all measured: qc 0, median
  h <- cloudnet_mlh_halfhour(mk_prof(rep("measured", 6), mlh = c(300, 500, 400, 600, 200, 700)), o)
  expect_equal(nrow(h), 1)
  expect_equal(h$pblh_m, 450)
  expect_equal(h$qc_flag, 0L)
  expect_equal(h$n_samples, 6L)
  expect_equal(h$x_coverage, 0.6)
  expect_equal(h$x_stare_seconds, 84)
  expect_equal(h$timestamp_start_utc, as.POSIXct("2022-06-01 12:00:00", tz = "UTC"))
  expect_equal(as.numeric(h$timestamp_end_utc - h$timestamp_start_utc, units = "secs"), 1800)
  expect_equal(h$x_first_profile_utc, "2022-06-01T12:00:14Z")
  expect_equal(h$x_last_profile_utc, "2022-06-01T12:01:24Z")
  expect_true(h$x_precip_checked)
  # insufficient coverage: 4 valid of 10 expected (no_signal and precipitation not valid)
  h <- cloudnet_mlh_halfhour(mk_prof(c(rep("measured", 4), "no_signal", "precipitation"), mlh = 300), o)
  expect_true(is.na(h$pblh_m))
  expect_true(is.na(h$qc_flag))
  expect_equal(h$x_n_no_signal, 1)
  expect_equal(h$x_n_precipitation, 1)
  # censored majority (4 of 7 valid): qc 1
  h <- cloudnet_mlh_halfhour(mk_prof(c(rep("measured", 3), rep("censored_weak", 2), rep("lowest_gate_missing", 2)),
                                     mlh = 300), o)
  expect_true(is.na(h$pblh_m))
  expect_equal(h$qc_flag, 1L)
  expect_equal(h$x_n_censored_weak, 2)
  expect_equal(h$x_n_lowest_gate_missing, 2)
  expect_equal(h$x_censor_below_m, 105)
  # exactly half censored: not a majority, median of the measured, qc 2
  h <- cloudnet_mlh_halfhour(mk_prof(c(rep("measured", 3), rep("censored_weak", 3)), mlh = c(200, 300, 900)), o)
  expect_equal(h$pblh_m, 300)
  expect_equal(h$qc_flag, 2L)
  # coverage capped at 1
  h <- cloudnet_mlh_halfhour(mk_prof(rep("measured", 15), mlh = 300), o)
  expect_equal(h$x_coverage, 1)
  # no classification for the day: x_precip_checked FALSE
  p <- mk_prof(rep("measured", 6), mlh = 300)
  p$precip <- NA
  expect_false(cloudnet_mlh_halfhour(p, o)$x_precip_checked)
  # two half-hours
  p <- dplyr::bind_rows(mk_prof(rep("measured", 6), mlh = 300),
                        mk_prof(rep("measured", 6), mlh = 500, start = "2022-06-01 12:30:00"))
  expect_equal(cloudnet_mlh_halfhour(p, o)$pblh_m, c(300, 500))
})

test_that("standardize_cloudnet_mlh writes a valid standard source", {
  h <- cloudnet_mlh_halfhour(mk_prof(rep("measured", 6), mlh = 300), cloudnet_mlh_options(z_min_agl = 105))
  ex <- read_exclusions(fixture_path("exclusion_periods.csv"))
  info <- list(instrument = "TEST HALO 1 (https://hdl.handle.net/x/halo1)", citation = "test citation",
               processing_history = "test.", instrument_altitude_m = 150,
               original_resolution = "profiles every 14 s")
  s <- suppressMessages(standardize_cloudnet_mlh(h, "XX-Tst", ex, 1, info))
  expect_silent(validate_standard_pblh(s$data, s$meta))
  expect_equal(s$meta$source_id, "cloudnet_eps_mlh")
  expect_equal(s$meta$method_class, "turbulence")
  expect_equal(s$meta$height_reference, "agl")
  expect_equal(s$meta$time_reference, "interval_mean")
  tmp <- withr::local_tempdir()
  suppressMessages(write_standard_pblh(s$data, s$meta, tmp))
  r <- suppressMessages(read_standard_pblh(file.path(tmp, "XX-Tst__cloudnet_eps_mlh.csv.gz")))
  expect_equal(r$pblh_m, 300)
  expect_equal(attr(r, "meta")$source_id, "cloudnet_eps_mlh")
})
