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
  expect_error(read_measured_sites(tmp), "missing column")
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
