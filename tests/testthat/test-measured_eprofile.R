# Spec 13: MeteoSwiss E-PROFILE files (CABAM, STRATfinder new bits), pairing table, BL-View L3 prototype.
# Synthetic NetCDF files written to a temporary folder.

# CABAM client file: 15-min steps with `time` = centre and `time_bnds` = [start, end).
write_cabam_nc <- function(path, starts, mlh, qual, qcf, dlh = rep(NaN, length(starts)), nraw = rep(60, length(starts)),
                           centre_shift = 0, masks = TRUE, licence = "TO BE CONFIRMED. Test.") {
  dim_t <- ncdf4::ncdim_def("time", "seconds since 1970-01-01 00:00:00 UTC", starts + 450 + centre_shift, unlim = FALSE)
  dim_nv <- ncdf4::ncdim_def("nv", "", 1:2, create_dimvar = FALSE)
  fv <- function(name, u = "m") ncdf4::ncvar_def(name, u, list(dim_t), missval = NaN, prec = "float")
  iv <- function(name) ncdf4::ncvar_def(name, "1", list(dim_t), missval = NULL, prec = "short")
  sv <- function(name, u) ncdf4::ncvar_def(name, u, list(), prec = "float")
  v <- list(time_bnds = ncdf4::ncvar_def("time_bnds", "seconds since 1970-01-01 00:00:00 UTC", list(dim_nv, dim_t),
                                         missval = NULL, prec = "double"),
            mlh = fv("mlh"), dlh = fv("dlh"), quality_level_mlh = iv("quality_level_mlh"), qc_flag_mlh = iv("qc_flag_mlh"),
            homogeneous_segment = iv("homogeneous_segment"), n_samples = fv("n_samples", "1"),
            window_transmission_min = fv("window_transmission_min", "%"), tilt_angle = fv("tilt_angle", "degree"),
            station_altitude = sv("station_altitude", "m"), sensor_height = sv("sensor_height", "m"),
            lat = sv("lat", "degrees_north"), lon = sv("lon", "degrees_east"))
  nc <- ncdf4::nc_create(path, unname(v), force_v4 = TRUE)
  ncdf4::ncvar_put(nc, v$time_bnds, rbind(starts, starts + 900))
  put <- list(mlh = mlh, dlh = dlh, quality_level_mlh = qual, qc_flag_mlh = qcf,
              homogeneous_segment = rep(0L, length(starts)), n_samples = nraw,
              window_transmission_min = rep(95, length(starts)), tilt_angle = rep(12, length(starts)))
  for (k in names(put)) ncdf4::ncvar_put(nc, v[[k]], put[[k]])
  for (k in c("station_altitude", "sensor_height", "lat", "lon")) {
    ncdf4::ncvar_put(nc, v[[k]], c(station_altitude = 100, sensor_height = 0, lat = 48.6, lon = 2.7)[[k]])
  }
  ncdf4::ncatt_put(nc, "mlh", "height_reference", "above ground level (a.g.l.) at the ceilometer site")
  # masks: TRUE = the CABAM table, FALSE = no attributes, or a named integer vector (meaning = bit)
  b <- if (isTRUE(masks)) meteoswiss_qc_bits("CABAM") else if (is.numeric(masks)) masks else NULL
  if (!is.null(b)) {
    ncdf4::ncatt_put(nc, "qc_flag_mlh", "flag_masks", as.integer(unname(b)), prec = "int")
    ncdf4::ncatt_put(nc, "qc_flag_mlh", "flag_meanings", paste(names(b), collapse = " "))
  }
  ga <- list(title = "CABAM test", product_version = "1.0 (test)", date_created = "2026-11-05T00:00:00Z",
             license = licence, acknowledgement = "Test acknowledgement.", instrument = "Vaisala CL31 (test).",
             source = "CABAM 2021a (test)", institution = "MeteoSwiss (test)", data_owner = "Test operator (E-PROFILE L1)",
             station_name = "Test", station_code = "TST", wigos_station_id = "0-0-0-1",
             known_limitations = "Test limitation.", instrument_events = "none documented")
  for (k in names(ga)) ncdf4::ncatt_put(nc, 0, k, ga[[k]])
  ncdf4::nc_close(nc)
  path
}

# STRATfinder client file (as the Davos fixture), with the E-PROFILE bit 8192.
write_strat_nc <- function(path, tt, mlh, qual, qcf, events = "") {
  dim_t <- ncdf4::ncdim_def("time", "seconds since 1970-01-01 00:00:00 UTC", tt, unlim = FALSE)
  fv <- function(name) ncdf4::ncvar_def(name, "m", list(dim_t), missval = NaN, prec = "float")
  bv <- function(name) ncdf4::ncvar_def(name, "1", list(dim_t), missval = NULL, prec = "short")
  sv <- function(name, u) ncdf4::ncvar_def(name, u, list(), prec = "float")
  v <- list(mlh = fv("mlh"), ablh = fv("ablh"), quality_level_mlh = bv("quality_level_mlh"),
            quality_level_ablh = bv("quality_level_ablh"), qc_flag_mlh = bv("qc_flag_mlh"),
            homogeneous_segment = bv("homogeneous_segment"), station_altitude = sv("station_altitude", "m"),
            sensor_height = sv("sensor_height", "m"), search_upper_limit = sv("search_upper_limit", "m"),
            lat = sv("lat", "degrees_north"), lon = sv("lon", "degrees_east"))
  nc <- ncdf4::nc_create(path, unname(v), force_v4 = TRUE)
  n <- length(tt)
  put <- list(mlh = mlh, ablh = rep(900, n), quality_level_mlh = qual, quality_level_ablh = rep(0L, n),
              qc_flag_mlh = qcf, homogeneous_segment = rep(0L, n))
  for (k in names(put)) ncdf4::ncvar_put(nc, v[[k]], put[[k]])
  for (k in c("station_altitude", "sensor_height", "search_upper_limit", "lat", "lon")) {
    ncdf4::ncvar_put(nc, v[[k]], c(station_altitude = 111, sensor_height = 0, search_upper_limit = 4000, lat = 50.9,
                                   lon = 6.4)[[k]])
  }
  ncdf4::ncatt_put(nc, "mlh", "height_reference", "above ground level (a.g.l.) at the ceilometer site")
  b <- meteoswiss_qc_bits("STRATfinder")
  ncdf4::ncatt_put(nc, "qc_flag_mlh", "flag_masks", unname(b), prec = "short")
  ncdf4::ncatt_put(nc, "qc_flag_mlh", "flag_meanings", paste(names(b), collapse = " "))
  ga <- list(title = "STRATfinder test", product_version = "1.0 (test)", date_created = "2026-11-05T00:00:00Z",
             license = "TO BE CONFIRMED. Test.", acknowledgement = "Test.", source = "STRATfinder 2021a (test)",
             institution = "MeteoSwiss (test)", data_owner = "Test operator", station_name = "Test",
             station_code = "TST", instrument_events = events)
  for (k in names(ga)) ncdf4::ncatt_put(nc, 0, k, ga[[k]])
  ncdf4::nc_close(nc)
  path
}

t0 <- as.numeric(as.POSIXct("2023-06-01", tz = "UTC"))

test_that("read_meteoswiss_cabam: interval start from time_bnds, alignment checked, alh not needed", {
  starts <- t0 + 900 * (0:3)
  f <- write_cabam_nc(tempfile(fileext = ".nc"), starts, c(100, 200, 300, NaN), c(0L, 0L, 0L, 9L), c(0L, 0L, 8L, 1L))
  r <- read_meteoswiss_cabam(f)
  expect_equal(as.numeric(r$interval_start_utc), starts)
  expect_equal(as.numeric(r$interval_end_utc - r$interval_start_utc, units = "secs"), rep(900, 4))
  expect_true(is.na(r$mlh_m[4]))
  m <- attr(r, "meteoswiss_meta")
  expect_equal(m$method, "CABAM")
  expect_equal(m$data_owner, "Test operator (E-PROFILE L1)")
  expect_equal(m$wigos_station_id, "0-0-0-1")
  # time not at the interval centre (e.g. CABAM's end stamp) -> abort
  f2 <- write_cabam_nc(tempfile(fileext = ".nc"), starts, rep(100, 4), rep(0L, 4), rep(0L, 4), centre_shift = 450)
  expect_error(read_meteoswiss_cabam(f2), "quarter hour")
  # intervals not on the quarter hour -> abort
  f3 <- write_cabam_nc(tempfile(fileext = ".nc"), starts + 60, rep(100, 4), rep(0L, 4), rep(0L, 4))
  expect_error(read_meteoswiss_cabam(f3), "quarter hour")
  expect_identical(read_meteoswiss_mlh(file.path(tempdir(), "nothing.nc")) |> tryCatch(error = function(e) "err"), "err")
})

test_that("cabam_halfhour: two steps per half-hour, mean of the quality-0 steps, qc 0 both / 1 one / NA none", {
  # half-hour 1: 100, 200 good -> 150, qc 0; 2: 300 good + 900 incomplete (bit 16384, quality 1) -> 300, qc 1;
  # 3: both bad (dirty window, NaN) -> NA; 4: only the second step in the file, good -> 400, qc 1
  starts <- t0 + 900 * c(0, 1, 2, 3, 4, 5, 7)
  mlh <- c(100, 200, 300, 900, NaN, NaN, 400)
  qual <- c(0L, 0L, 0L, 1L, 2L, 2L, 0L)
  qcf <- c(8L, 0L, 0L, 16384L, 2L, 2L, 0L)
  dlh <- c(NaN, 500, NaN, NaN, NaN, NaN, 700)
  f <- write_cabam_nc(tempfile(fileext = ".nc"), starts, mlh, qual, qcf, dlh = dlh, nraw = c(60, 58, 60, 12, 60, 60, 60))
  h <- cabam_halfhour(read_meteoswiss_cabam(f))
  expect_equal(as.numeric(h$timestamp_start_utc), t0 + 1800 * (0:3))
  expect_equal(h$pblh_m, c(150, 300, NA, 400))
  expect_equal(h$qc_flag, c(0L, 1L, NA, 1L))
  expect_equal(h$n_samples, c(2L, 1L, 0L, 1L))
  expect_equal(h$x_coverage, c(1, 0.5, 0, 0.5))
  expect_equal(h$x_share_lowest, c(0.5, 0, NA, 0))
  expect_equal(h$x_dlh_m, c(500, NA, NA, 700))
  expect_equal(h$x_n_steps, c(2L, 2L, 2L, 1L))
  expect_equal(h$x_n_profiles, c(118, 60, 0, 60))
})

test_that("CABAM: unknown bits abort (declared or set)", {
  starts <- t0 + 900 * (0:1)
  expect_error(read_meteoswiss_cabam(write_cabam_nc(tempfile(fileext = ".nc"), starts, c(1, 2), c(0L, 0L), c(4L, 0L),
                                                    masks = FALSE)), "unknown bits")
  f <- write_cabam_nc(tempfile(fileext = ".nc"), starts, c(1, 2), c(0L, 0L), c(0L, 0L),
                      masks = c(no_retrieval = 1L, new_test = 32768L))
  expect_error(read_meteoswiss_cabam(f), "new_test")
})

test_that("standardize_meteoswiss_mlh (CABAM): source id, sidecar owner / processor / licence / citation", {
  starts <- t0 + 900 * (0:3)
  f <- write_cabam_nc(tempfile(fileext = ".nc"), starts, c(100, 200, 300, 400), rep(0L, 4), rep(0L, 4))
  ex <- read_exclusions(fixture_path("exclusion_periods.csv"))
  s <- standardize_meteoswiss_mlh(read_meteoswiss_cabam(f), "XX-Tst", ex, 1)
  expect_equal(unique(s$data$source_id), "meteoswiss_cabam_mlh")
  expect_equal(s$data$pblh_m, c(150, 350))
  expect_equal(s$meta$data_owner, "Test operator (E-PROFILE L1)")
  expect_match(s$meta$processor, "^MeteoSwiss \\(CABAM 2021a; client product version 1.0")
  expect_match(s$meta$licence, "^to be confirmed: internal use only, contact the provider before sharing results")
  expect_match(s$meta$citation, "Kotthaus2018 \\(Part 1\\) and Kotthaus2018b")
  expect_equal(s$meta$known_limitations, "Test limitation.")
  expect_equal(s$meta$method_class, "aerosol")
  expect_true(validate_standard_pblh(s$data, s$meta))
  p <- write_standard_pblh(s$data, s$meta, tempfile("std_"))
  y <- yaml::read_yaml(p[["meta"]])
  expect_equal(y$data_owner, "Test operator (E-PROFILE L1)")
  expect_equal(y$metadata_completeness, round(11 / 12, 4)) # contact unknown; the extra fields do not count
})

test_that("STRATfinder E-PROFILE: bit 8192 informational (kept, share reported), events parsed", {
  tt <- t0 + 30 + 60 * (0:29)
  qcf <- c(rep(8192L, 20), rep(8192L + 8L, 5), rep(8192L + 4L, 5)) # bit 4 -> quality 1
  qual <- c(rep(0L, 25), rep(1L, 5))
  ev <- paste("2021-03-18T00:00 to 2026-01-31T23:59 [no_temperature_overlap_model] no model (QC)",
              "2023-06-01T00:00 to 2023-06-01T23:59 [intervention] Window transmission below 60 % (L1)",
              "2025-10-21T00:00 [breakpoint] firmware 0.743 -> 1.13", sep = " | ")
  f <- write_strat_nc(tempfile(fileext = ".nc"), tt, c(rep(300, 25), rep(5000, 5)), qual, qcf, events = ev)
  r <- read_meteoswiss_stratfinder(f)
  h <- meteoswiss_halfhour(r)
  expect_equal(h$pblh_m, 300)
  expect_equal(h$n_samples, 25L)
  expect_equal(h$qc_flag, 0L)
  expect_equal(h$x_share_no_temp_overlap, 1)
  expect_equal(h$x_share_lowest, 0.2)
  e <- meteoswiss_events(attr(r, "meteoswiss_meta")$instrument_events)
  expect_equal(e$category, c("no_temperature_overlap_model", "intervention", "breakpoint"))
  expect_true(is.na(e$end_utc[3]))
  expect_equal(format(e$end_utc[1], "%Y-%m-%d %H:%M"), "2026-01-31 23:59")
  # CABAM-style events (no category, date only)
  e2 <- meteoswiss_events("2021-10-26 new correction period | 2025-11-05 firmware 205 -> 321")
  expect_equal(format(e2$start_utc, "%Y-%m-%d"), c("2021-10-26", "2025-11-05"))
  expect_true(all(is.na(e2$category)))
  expect_equal(nrow(meteoswiss_events("none documented")), 0)
  s <- standardize_meteoswiss_mlh(r, "XX-Tst", read_exclusions(fixture_path("exclusion_periods.csv")), 1)
  expect_equal(s$data$source_id, "meteoswiss_stratfinder_mlh")
  expect_match(s$meta$licence, "^to be confirmed")
  expect_match(s$meta$processor, "STRATfinder 2021a")
  # a bit outside the STRATfinder table aborts
  f2 <- write_strat_nc(tempfile(fileext = ".nc"), tt, rep(300, 30), rep(0L, 30), c(16384L, rep(0L, 29)))
  expect_error(read_meteoswiss_stratfinder(f2), "unknown bits")
})

test_that("distance classes and the pairing table (spec 13)", {
  expect_equal(distance_class(c(1.8, 4.99, 5, 5.3, 16.6, 24.87, 25, 25.01, 50, 50.1, NA)),
               c("all_regimes", "all_regimes", "daytime_convective", "daytime_convective", "daytime_convective",
                 "daytime_convective", "daytime_convective", "daytime_climatology", "daytime_climatology", "not_used",
                 NA))
  ms <- read_measured_sites(here::here("config", "measured_sites.csv"))
  me <- ms[ms$provider == "MeteoSwiss" & ms$icos_site != "CH-Dav", ]
  me <- me[match(c("DE-RuS", "DE-HoH", "FI-Lom", "FR-Fon", "FR-Lus"), me$icos_site), ]
  expect_equal(me$provider_site, c("JUE", "MAG", "KEN", "MEL", "POI"))
  expect_equal(me$distance_class, c("daytime_convective", "daytime_convective", "all_regimes", "daytime_convective",
                                    "daytime_convective"))
  # FI-Lom (KEN) inventory since G 2026-11-12: no overlap with the EC record (spec 14)
  expect_equal(me$primary_source_id, c(rep("meteoswiss_stratfinder_mlh", 2), NA, rep("meteoswiss_cabam_mlh", 2)))
  expect_equal(me$status, c("used", "used", "inventory", "used", "used"))
  # spec 14: secondary sources, never primary
  expect_equal(ms$status[ms$icos_site == "CZ-Lnz"], "secondary")
  expect_equal(ms$status[ms$icos_site == "DE-Tha"], "secondary")
  expect_true(all(is.na(ms$primary_source_id[ms$status == "secondary"])))
  expect_equal(ms$distance_class[ms$icos_site == "DE-Tha"], "all_regimes")
  expect_true(all(is.na(ms$primary_source_id[ms$provider == "cloudnet" & ms$icos_site %in% c("FI-Lom", "DE-RuS")])))
  # validation: wrong class, two primaries, primary on a pending row, unknown status
  bad <- function(edit) {
    x <- readr::read_csv(here::here("config", "measured_sites.csv"), col_types = readr::cols(.default = "c"),
                         na = c("NA", ""))
    x <- edit(x)
    p <- tempfile(fileext = ".csv")
    readr::write_csv(x, p, na = "NA")
    p
  }
  expect_error(read_measured_sites(bad(function(x) { x$distance_class[1] <- "not_used"; x })), "distance rule")
  expect_error(read_measured_sites(bad(function(x) {
    x$primary_source_id[x$icos_site == "DE-RuS" & x$provider == "cloudnet"] <- "cloudnet_eps_mlh"
    x$status[x$icos_site == "DE-RuS" & x$provider == "cloudnet"] <- "used"; x
  })), "primary source at")
  expect_error(read_measured_sites(bad(function(x) {
    x$primary_source_id[x$icos_site == "CZ-Lnz"] <- "blview_mlh"; x
  })), "not 'used'")
  expect_error(read_measured_sites(bad(function(x) { x$status[1] <- "maybe"; x })), "status")
  # borderline flag (G 2026-11-06): DE-HoH and FR-Lus only; a table without the column gives FALSE
  expect_equal(me$borderline, c(FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_equal(sum(ms$borderline), 2)
  nb <- suppressWarnings(read_measured_sites(bad(function(x) { x$borderline <- NULL; x })))
  expect_false(any(nb$borderline))
})

test_that("tag_measured_pairing adds the pairing distance and the borderline flag", {
  ms <- read_measured_sites(here::here("config", "measured_sites.csv"))
  row <- ms[ms$icos_site == "FR-Lus" & !is.na(ms$primary_source_id), ]
  x <- tag_measured_pairing(dplyr::tibble(pair = c("a", "b"), n = 1:2), "FR-Lus", row)
  expect_equal(names(x)[1], "site_id")
  expect_equal(x$borderline, c(TRUE, TRUE))
  expect_equal(x$measured_distance_km, c(24.87, 24.87))
  expect_equal(x$measured_provider_site, c("POI", "POI"))
  x0 <- tag_measured_pairing(dplyr::tibble(pair = character()), "FR-Lus", row)
  expect_equal(nrow(x0), 0)
  expect_error(tag_measured_pairing(x, "FR-Lus", ms[0, ]), "exactly one")
})

test_that("measured-source exclusions: table, validation, applied by the standard-file builder", {
  mex <- read_measured_exclusions(here::here("config", "measured_exclusions.csv"))
  expect_equal(sort(unique(mex$site_id)), c("FR-Fon", "FR-Lus"))
  expect_true(all(mex$source_id == "meteoswiss_cabam_mlh"))
  expect_equal(format(mex$date_start[mex$site_id == "FR-Fon"], "%Y-%m"),
               c("2023-05", "2023-06", "2025-09", "2026-08", "2026-09", "2026-10"))
  fl <- mex[mex$site_id == "FR-Lus", ]
  expect_equal(c(format(fl$date_start), format(fl$date_end)), c("2021-02-01", "2021-10-31"))
  wr <- function(x) { p <- tempfile(fileext = ".csv"); readr::write_csv(x, p); p }
  expect_error(read_measured_exclusions(wr(dplyr::mutate(mex, source_id = "nope"))), "unknown")
  expect_error(read_measured_exclusions(wr(dplyr::mutate(mex, date_end = date_start - 1))), "before")
  expect_error(read_measured_exclusions(wr(dplyr::mutate(mex, reason = NA))), "required")
  # interval [date_start 00:00, date_end + 1 d 00:00) UTC on timestamp_start_utc; other site / source untouched
  m1 <- dplyr::tibble(site_id = c("XX-Tst", "XX-Tst"), source_id = c("meteoswiss_cabam_mlh", "sirta_mlh_cl31"),
                      date_start = as.Date(c("2023-06-01", "2023-06-01")), date_end = as.Date(c("2023-06-01", "2023-06-02")),
                      reason = "test")
  d <- dplyr::tibble(timestamp_start_utc = as.POSIXct(c("2023-05-31 23:30", "2023-06-01 00:00", "2023-06-01 23:30",
                                                        "2023-06-02 00:00"), tz = "UTC"), pblh_m = 1:4)
  a <- suppressMessages(apply_measured_exclusions(d, "XX-Tst", "meteoswiss_cabam_mlh", m1))
  expect_equal(a$data$pblh_m, c(1, 4))
  expect_equal(a$n_excluded, 2)
  expect_match(a$history, "dropped \\(2; periods: 2023-06-01 to 2023-06-01 \\(test\\)\\)")
  expect_equal(suppressMessages(apply_measured_exclusions(d, "YY-Tst", "meteoswiss_cabam_mlh", m1))$n_excluded, 0)
  expect_equal(suppressMessages(apply_measured_exclusions(d, "XX-Tst", "meteoswiss_cabam_mlh", NULL))$periods, "none")
  # builder: the CABAM half-hours of 2023-06-01 dropped, counted, in the sidecar history
  starts <- t0 + 900 * (0:3)
  f <- write_cabam_nc(tempfile(fileext = ".nc"), starts, c(100, 200, 300, 400), rep(0L, 4), rep(0L, 4))
  s <- suppressMessages(standardize_meteoswiss_mlh(read_meteoswiss_cabam(f), "XX-Tst",
                                                   read_exclusions(fixture_path("exclusion_periods.csv")), 1, m1))
  expect_equal(nrow(s$data), 0)
  expect_equal(s$n_excluded_measured, 2)
  expect_match(s$meta$processing_history, "measured-source exclusion periods")
})

# BL-View L3 file: hourly values repeated on ~1-min profiles (time_aerosol), fill -999.
write_blview_nc <- function(path, tt, h, q) {
  dim_t <- ncdf4::ncdim_def("time_aerosol", "seconds since 1970-01-01 00:00:00.000", tt, unlim = FALSE)
  v <- list(h = ncdf4::ncvar_def("mean_mixing_layer_height", "m", list(dim_t), missval = -999, prec = "integer"),
            q = ncdf4::ncvar_def("mean_mixing_layer_quality", "", list(dim_t), missval = -999, prec = "integer"),
            o = ncdf4::ncvar_def("time_zone", "hours", list(dim_t), missval = -999, prec = "integer"))
  nc <- ncdf4::nc_create(path, unname(v), force_v4 = TRUE)
  ncdf4::ncvar_put(nc, v$h, h)
  ncdf4::ncvar_put(nc, v$q, q)
  ncdf4::ncvar_put(nc, v$o, rep(2L, length(tt)))
  ncdf4::ncatt_put(nc, 0, "sw_version", "2.3.2")
  ncdf4::ncatt_put(nc, 0, "instrument_serial_number", "T0")
  ncdf4::nc_close(nc)
  path
}

test_that("read_blview_l3: one row per hour, fill to NA, time offset applied, offset required", {
  d0 <- as.numeric(as.POSIXct("2026-10-07", tz = "UTC"))
  tt <- d0 + 50 + 60 * (0:179) # three hours of 1-min profiles at hh:mm:50
  h <- c(rep(990L, 60), rep(-999L, 60), rep(780L, 60))
  q <- c(rep(3L, 60), rep(-999L, 60), rep(2L, 60))
  f <- write_blview_nc(tempfile(fileext = ".nc"), tt, h, q)
  x <- read_blview_l3(f, time_offset_h = 0)
  # spec 14: profiles hh:00-hh:59 = the hour ending at hh:00
  expect_equal(as.numeric(x$hour_end_utc), d0 + 3600 * (0:2))
  expect_equal(x$mlh_m, c(990, NA, 780))
  expect_equal(x$quality, c(3L, NA, 2L))
  expect_equal(x$n_profiles, c(60L, 60L, 60L))
  expect_equal(x$version, rep("3.2.2", 3))
  expect_true(all(x$consistent))
  expect_equal(attr(x, "blview_meta")$sw_version, "2.3.2")
  x2 <- read_blview_l3(f, time_offset_h = 2) # stored time in CEST
  expect_equal(as.numeric(x2$hour_end_utc), d0 - 7200 + 3600 * (0:2))
  expect_error(read_blview_l3(f), "time_offset_h")
  # spec 14: a value that changes within an hour is reported (NA, consistent FALSE), not an error
  h[30] <- 500L
  y <- read_blview_l3(write_blview_nc(tempfile(fileext = ".nc"), tt, h, q), 0)
  expect_equal(y$consistent, c(FALSE, TRUE, TRUE))
  expect_equal(y$n_values, c(2L, 1L, 1L))
  expect_true(is.na(y$mlh_m[1]) && is.na(y$quality[1]))
  expect_equal(attr(y, "blview_meta")$n_inconsistent, 1L)
})
