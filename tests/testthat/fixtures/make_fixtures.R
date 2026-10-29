# Rebuild the synthetic ICOS L2, ERA5, Cloudnet and ONEFlux BIF test fixtures (tiny, no real data).
# Run from the project root: Rscript tests/testthat/fixtures/make_fixtures.R [icos] [era5] [cloudnet] [oneflux] [reobs]
# (no argument = all parts; rebuilding the ICOS zips changes their bytes even
# when the content is unchanged, so rebuild only the part you need).
# Requires a `zip` program for utils::zip() (on Windows: Rtools).
# Synthetic site XX-Tst, UTC+1.

parts <- commandArgs(trailingOnly = TRUE)
if (length(parts) == 0) parts <- c("icos", "era5", "cloudnet", "oneflux", "reobs")

fx <- here::here("tests", "testthat", "fixtures")
site <- "XX-Tst"

# Half-hour timestamp strings, local standard time.
ts <- function(from, n) {
  t0 <- as.POSIXct(from, format = "%Y%m%d%H%M", tz = "UTC")
  s <- t0 + (seq_len(n) - 1) * 1800
  data.frame(
    TIMESTAMP_START = format(s, "%Y%m%d%H%M"),
    TIMESTAMP_END = format(s + 1800, "%Y%m%d%H%M")
  )
}

write_zip <- function(df, dir, csv_name) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(tempdir(), csv_name)
  utils::write.csv(df, tmp, row.names = FALSE, quote = FALSE)
  zipfile <- file.path(dir, sub("\\.csv$", ".zip", csv_name))
  unlink(zipfile)
  utils::zip(zipfile, tmp, flags = "-j9Xq")
  unlink(tmp)
  zipfile
}

if ("icos" %in% parts) {
  l2_dir <- file.path(fx, "icos_l2", paste0(site, "_output_L2I_2026"))

  # FLUXES: 2020-01-01 00:00 to 03:00 local, the 01:00 record missing (6 rows).
  fl <- ts("202001010000", 7)[-3, ]
  n <- nrow(fl)
  fluxes_cols <- c(
    "H", "H_UNCLEANED", "H_SSITC_TEST", "USTAR", "MO_LENGTH", "ZL", "T_SONIC", "WS", "WD",
    "U_SIGMA", "V_SIGMA", "W_SIGMA", "TAU", "TAU_SSITC_TEST", "PBLH", "FETCH_MAX",
    "FETCH_50", "FETCH_70", "FETCH_80", "FETCH_90", "FOOTPRINT_TA_CONTR", "FOOTPRINT_FLAG"
  )
  fl$CO2 <- 400 # extra column, must not be read
  for (v in fluxes_cols) fl[[v]] <- seq_len(n)
  fl$H <- c(-9999, -50, 20, 100, 900, 10)        # -9999 -> NA; 900 out of range
  fl$T_SONIC <- c(283.15, 283.65, 284.15, 284.65, 285.15, 285.65)  # K (as in L2 FLUXES)
  fl$USTAR <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  fl$MO_LENGTH <- c(-100, 50, 0, -9999, 200, -20)
  fl$PBLH <- c(500, "-9999.0", 800, 1000, 5, 1200)
  invisible(write_zip(fl, l2_dir, paste0("ICOSETC_", site, "_FLUXES_INTERIM_L2.csv")))

  # FLUXNET: 00:00 to 01:30 local (4 rows, ends before FLUXES).
  fn <- ts("202001010000", 4)
  fn$TA_F_MDS <- 1 # extra column
  fn$TA_F <- c(5, 5.5, 6, 6.5)                   # degC
  fn$TA_F_QC <- 0
  fn$PA_F <- c(100.5, 100.4, 100.3, 100.2)       # kPa
  fn$PA_F_QC <- 0
  fn$H_F_MDS <- c(-30, -20, 10, 50)
  fn$H_F_MDS_QC <- c(0, 0, 1, 2)
  fn$H_CORR <- c(-35, -25, 12, -9999)
  fn$NIGHT <- c(1, 1, 1, 0)
  fn$SW_IN_POT <- c(0, 0, 0, 10)
  fn$WS_F <- c(2, 2.1, 2.2, 2.3)
  fn$P_F <- 0
  invisible(write_zip(fn, l2_dir, paste0("ICOSETC_", site, "_FLUXNET_FLUXMET_HH_2020-2020_v1.0_r1.csv")))

  # METEO: 00:00 to 03:30 local (8 rows), TA_10 absent.
  me <- ts("202001010000", 8)
  me$RH <- 80 # extra column
  me$TA <- c(2, 2.5, "-9999.0", 3.5, 4, 4.5, 5, 5.5) # degC
  me$PA <- 99.5                                    # kPa
  me$P <- 0
  me$SW_IN <- c(0, 0, 0, 0, 5, 10, 20, 30)
  me$WS <- 3
  for (i in 1:9) me[[paste0("TA_", i)]] <- 1 + i / 10
  invisible(write_zip(me, l2_dir, paste0("ICOSETC_", site, "_METEO_INTERIM_L2.csv")))
}

# ---- ERA5, format of the Oct 2026 delivery (CLAUDE.md; spec 04 section 10) ----
# 3x3 cells around the tower (XX-Tst at 50.0 N, 10.0 E = centre of the tower cell).
# Hourly single levels: 3 variables x 9 cells x 48 h (direction flags in the name);
# pressure levels: T and z x 9 cells x 8 times (6-hourly) x 3 levels; static z, lsm,
# sdor x 9 cells (signed offsets in the name). Known values used by the tests:
# - sshf at the tower cell, 2020-01-01 12:00 UTC = -360000 J m-2 (H = +100 W m-2);
# - T = 280 K at 850 hPa, tower cell, 2020-01-01 00:00 (theta = 293.30 K);
# - zg = 100 / 1500 / 3000 m at 1000 / 850 / 700 hPa (+ 2 m per cell row);
# - orography 100 + 10 row + 5 col m; lsm 0.6 in cell (row +1, col -1), else 1.
if ("era5" %in% parts) {
  g <- 9.80665
  lat0 <- 50.0
  lon0 <- 10.0
  era_root <- file.path(fx, "era5_2026")
  unlink(era_root, recursive = TRUE)
  sl_dir <- file.path(era_root, "sl", site)
  pl_dir <- file.path(era_root, "pl", site)
  st_dir <- file.path(era_root, "static")
  for (d in c(sl_dir, pl_dir, st_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  cells <- expand.grid(row = -1:1, col = -1:1)
  flags <- function(r, c) {
    sprintf("north_+%d__south_+%d__east_+%d__west_+%d", r == 1, r == -1, c == 1, c == -1)
  }
  signed <- function(r, c) sprintf("north_%+d__south_+0__east_%+d__west_+0", r, c)
  fmt_t <- function(t) format(t, "%Y-%m-%d %H:%M:%S")
  write_era5 <- function(path, fields, df, title = "# ==== METADATA RICHIESTA E DOWNLOAD ====") {
    writeLines(c(title, paste0("# ", names(fields), ": ", unlist(fields)),
                 "# ==== FINE METADATA (dati sotto) ===="), path)
    suppressWarnings(utils::write.table(df, path, sep = ",", row.names = FALSE, quote = FALSE,
                                        append = TRUE, na = ""))
  }
  arco_fields <- function(dataset, var, r, c, t) {
    list(
      script = "era5_download_arco.py", point_id = site, method = "arco", dataset = dataset,
      variable = var, latitude_requested = lat0, longitude_requested = lon0,
      longitude_conversion = "non necessaria (synthetic fixture)",
      latitude_grid = lat0 + 0.25 * r, longitude_grid = lon0 + 0.25 * c,
      grid_lat_index = 160 - r, grid_lon_index = 760 + c, grid_resolution_deg = 0.25,
      pixel_offset_north = as.integer(r == 1), pixel_offset_south = as.integer(r == -1),
      pixel_offset_east = as.integer(c == 1), pixel_offset_west = as.integer(c == -1),
      start_date_requested = format(min(t), "%Y-%m-%d"), end_date_requested = format(max(t), "%Y-%m-%d"),
      download_utc = "2026-10-14T09:00:00+00:00",
      n_rows = length(t), first_time = fmt_t(min(t)), last_time = fmt_t(max(t))
    )
  }

  # Hourly single levels, 2020-01-01 00:00 to 2020-01-02 23:00 UTC.
  t_h <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + 3600 * (0:47)
  solar_day <- function(t) { # 0..1, peak at local solar noon (UTC + lon / 15)
    hod <- as.numeric(format(t, "%H")) + as.numeric(format(t, "%M")) / 60
    pmax(0, cos(2 * pi * ((hod + lon0 / 15) - 12) / 24))
  }
  sl_name <- function(var, r, c) {
    file.path(sl_dir, sprintf("%s__reanalysis_era5_single_level__%s__20200101T0000__20200102T2300__%s.csv",
                              site, var, flags(r, c)))
  }
  for (i in seq_len(nrow(cells))) {
    r <- cells$row[i]
    c <- cells$col[i]
    blh <- round(200 + 800 * solar_day(t_h) + 40 * r + 10 * c, 3)
    cbh <- ifelse(seq_along(t_h) %% 4 == 0, NA, 1500 + 100 * r)
    h_up <- round(-20 + 220 * solar_day(t_h - 1800), 3) # W m-2, positive up, hour ending at t
    sshf <- -h_up * 3600                                 # J m-2, positive down
    if (r == 0 && c == 0) sshf[t_h == as.POSIXct("2020-01-01 12:00:00", tz = "UTC")] <- -360000
    vals <- list(boundary_layer_height = blh, cloud_base_height = cbh, surface_sensible_heat_flux = sshf)
    for (v in names(vals)) {
      df <- data.frame(time = fmt_t(t_h), x = vals[[v]])
      names(df)[2] <- v
      write_era5(sl_name(v, r, c), arco_fields("reanalysis_era5_single_level", v, r, c, t_h), df)
    }
  }
  # Ignored: a variable not in the registry, and a non-ERA5 file.
  df <- data.frame(time = fmt_t(t_h), "2m_temperature" = 280, check.names = FALSE)
  write_era5(sl_name("2m_temperature", 0, 0),
             arco_fields("reanalysis_era5_single_level", "2m_temperature", 0, 0, t_h), df)
  writeLines("not an ERA5 file", file.path(sl_dir, "README.txt"))

  # Pressure levels, 6-hourly, 2020-01-01 00:00 to 2020-01-02 18:00 UTC.
  t_p <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + 6 * 3600 * (0:7)
  lev <- c(700, 850, 1000)
  k <- seq_along(t_p) - 1
  for (i in seq_len(nrow(cells))) {
    r <- cells$row[i]
    c <- cells$col[i]
    tk <- c(`700` = 272, `850` = 280, `1000` = 285)
    zg <- c(`700` = 3000, `850` = 1500, `1000` = 100)
    for (v in c("temperature", "geopotential")) {
      df <- data.frame(time = fmt_t(t_p))
      for (p in lev) {
        pc <- as.character(p)
        df[[sprintf("%s_%dhPa", v, p)]] <- if (v == "temperature") {
          tk[[pc]] + 0.5 * r + 0.1 * k
        } else {
          (zg[[pc]] + 2 * r) * g
        }
      }
      fields <- c(arco_fields("reanalysis_era5_pressure_levels-pl", v, r, c, t_p),
                  list(levels = paste(lev, collapse = ","), level_dim = "pressureLevel"))
      write_era5(file.path(pl_dir, sprintf(
        "%s__reanalysis_era5_pressure_levels-pl__%s__20200101T0000__20200102T1800__%s.csv", site, v, flags(r, c)
      )), fields, df)
    }
  }

  # Static fields (valid time 1979-01-20 12:00), signed offsets in the name.
  st_units <- c(z = "m**2 s**-2", lsm = "(0 - 1)", sdor = "m")
  st_file <- function(s, v, r, c) {
    file.path(st_dir, sprintf("%s__single-levels__%s__197901201200__197901201200__%s.csv", s, v, signed(r, c)))
  }
  for (i in seq_len(nrow(cells))) {
    r <- cells$row[i]
    c <- cells$col[i]
    val <- c(z = (100 + 10 * r + 5 * c) * g, lsm = if (r == 1 && c == -1) 0.6 else 1, sdor = 5 + r + c)
    for (v in names(val)) {
      fields <- list(
        script = "extract_points.py", point_id = site, method = "eccodes_nearest", dataset = "single-levels",
        variable = v, units = st_units[[v]], latitude_requested = lat0, longitude_requested = lon0,
        longitude_conversion = "not needed (synthetic fixture)",
        latitude_grid = lat0 + 0.25 * r, longitude_grid = lon0 + 0.25 * c,
        distance_km = round(111.2 * sqrt((0.25 * r)^2 + (0.25 * c * cos(lat0 * pi / 180))^2), 3),
        grid_lat_index = 40 - r, grid_lon_index = 40 + c, grid_resolution_deg = 0.25,
        grid_index_origin = "north-west corner of the grid, 0-based",
        grid_scanning = "i: west->east, j: north->south",
        pixel_offset_north = 1, pixel_offset_south = 1, pixel_offset_east = 1, pixel_offset_west = 1,
        pixel_rows_from_center = sprintf("%+d (+ = north, - = south)", r),
        pixel_cols_from_center = sprintf("%+d (+ = east, - = west)", c),
        n_rows = 1, first_time = "1979-01-20 12:00:00", last_time = "1979-01-20 12:00:00",
        levels = "n/a", level_dim = "surface"
      )
      df <- data.frame(valid_time_utc = "1979-01-20 12:00", x = val[[v]])
      names(df)[2] <- v
      write_era5(st_file(site, v, r, c), fields, df, title = "# ==== METADATA RICHIESTA E ESTRAZIONE ====")
    }
  }
  # Another station in the shared static folder: silently ignored.
  writeLines("# other station", st_file("XX-Oth", "z", 0, 0))
}

# ---- Cloudnet API, /api/files response (spec 07): synthetic, structure as the real API ----
# 4 files: 2 epsilon-lidar days in 2022-01 (one volatile), 1 in 2022-03 (legacy), 1 of another instrument.
if ("cloudnet" %in% parts) {
  cn_dir <- file.path(fx, "cloudnet")
  dir.create(cn_dir, recursive = TRUE, showWarnings = FALSE)
  mk <- function(uuid, date, inst, pid, size, volatile, legacy, cov) {
    list(
      uuid = uuid, filename = paste0(gsub("-", "", date), "_testsite_epsilon-lidar_", substr(uuid, 1, 8), ".nc"),
      version = if (volatile) "" else "v1", pid = if (volatile) "" else paste0("https://hdl.handle.net/x/", uuid),
      volatile = volatile, legacy = legacy, measurementDate = date,
      checksum = strrep(substr(uuid, 1, 1), 64), size = as.character(size), coverage = cov,
      format = "HDF5 (NetCDF4)", errorLevel = "info",
      site = list(id = "testsite", latitude = 50.0, longitude = 10.0, altitude = 100),
      product = list(id = "epsilon-lidar", humanReadableName = "Turbulence (lidar)"),
      instrument = list(uuid = paste0(inst, "-uuid"), pid = pid, name = inst, type = "Doppler lidar"),
      downloadUrl = paste0("https://example.invalid/api/download/product/", uuid, ".nc")
    )
  }
  files <- list(
    mk("aaaa0001", "2022-01-01", "TEST HALO 1", "https://hdl.handle.net/x/halo1", 1000000, FALSE, FALSE, 0.5),
    mk("bbbb0002", "2022-01-02", "TEST HALO 1", "https://hdl.handle.net/x/halo1", 2000000, TRUE, FALSE, 0.25),
    mk("cccc0003", "2022-03-05", "TEST HALO 1", "https://hdl.handle.net/x/halo1", 500000, FALSE, TRUE, 1),
    mk("dddd0004", "2022-01-02", "TEST HALO 2", "https://hdl.handle.net/x/halo2", 300000, FALSE, FALSE, 0.1)
  )
  jsonlite::write_json(files, file.path(cn_dir, "api_files_sample.json"), auto_unbox = TRUE, pretty = TRUE)

  # epsilon-lidar file (structure of the Hyytiala files): 6 gates, centres 165..315 m a.s.l. (30 m spacing),
  # altitude 150 m (centres 15..165 m agl); 4 profiles at 00:00:15, 00:00:30, 00:20:00, 00:40:00 UTC.
  # Profile 1: eps >= 1e-4 at 15..105 m agl (gates 1-4), then 5e-5 -> MLH 120 m (gates below z_min used only
  # if z_min is lower); profile 2: all missing; profile 3: gate 4 (105 m) missing, gate 5 valid;
  # profile 4: gate 4 = 1e-5 (weak).
  fill <- 9.96920996838687e+36
  eps <- matrix(fill, nrow = 6, ncol = 4)
  eps[, 1] <- c(fill, fill, fill, 2e-4, 5e-5, 1e-3)
  eps[, 3] <- c(fill, fill, fill, fill, 3e-4, 3e-4)
  eps[, 4] <- c(fill, fill, fill, 1e-5, 2e-4, 2e-4)
  dim_t <- ncdf4::ncdim_def("time", "hours since 2022-01-01 00:00:00 +00:00", c(15, 30, 1200, 2400) / 3600)
  dim_z <- ncdf4::ncdim_def("height", "m", 165 + 30 * (0:5))
  v_eps <- ncdf4::ncvar_def("epsilon", "m2 s-3", list(dim_z, dim_t), missval = fill, prec = "float")
  v_alt <- ncdf4::ncvar_def("altitude", "m", list(), prec = "float")
  v_roll <- ncdf4::ncvar_def("rolling_window_period", "s", list(), prec = "float")
  f_eps <- file.path(cn_dir, "20220101_testsite_epsilon-lidar_aaaa0001.nc")
  unlink(f_eps)
  nc <- ncdf4::nc_create(f_eps, list(v_eps, v_alt, v_roll), force_v4 = TRUE)
  ncdf4::ncvar_put(nc, v_eps, eps)
  ncdf4::ncvar_put(nc, v_alt, 150)
  ncdf4::ncvar_put(nc, v_roll, 180)
  ncdf4::ncatt_put(nc, "height", "standard_name", "height_above_mean_sea_level")
  ncdf4::ncatt_put(nc, 0, "file_uuid", "aaaa0001-test")
  ncdf4::ncatt_put(nc, 0, "instrument_pid", "https://hdl.handle.net/x/halo1")
  ncdf4::ncatt_put(nc, 0, "references", "https://doi.org/10.1175/2010JTECHA1455.1")
  ncdf4::nc_close(nc)

  # classification file: 3 gates at 252, 300, 600 m a.s.l. (102, 150, 450 m agl), 5 profiles every 30 s from
  # 00:00:00; precipitation (class 2) at 150 m agl at 00:01:00 only; class 2 at 450 m agl at 00:00:00 (above
  # the checked layer); cloud base 1150 m a.s.l. at 00:00:00 and 00:00:30, none after; rain_detected at 00:01:30.
  tc <- matrix(0L, nrow = 3, ncol = 5)
  tc[3, 1] <- 2L
  tc[2, 3] <- 2L
  dim_tc <- ncdf4::ncdim_def("time", "hours since 2022-01-01 00:00:00 +00:00", c(0, 30, 60, 90, 120) / 3600)
  dim_zc <- ncdf4::ncdim_def("height", "m", c(252, 300, 600))
  v_tc <- ncdf4::ncvar_def("target_classification", "1", list(dim_zc, dim_tc), missval = -2147483647L, prec = "integer")
  v_cb <- ncdf4::ncvar_def("cloud_base_height_amsl", "m", list(dim_tc), missval = fill, prec = "float")
  v_rd <- ncdf4::ncvar_def("rain_detected", "1", list(dim_tc), missval = -2147483647L, prec = "integer")
  v_altc <- ncdf4::ncvar_def("altitude", "m", list(dim_tc), prec = "float")
  f_cls <- file.path(cn_dir, "20220101_testsite_classification.nc")
  unlink(f_cls)
  nc <- ncdf4::nc_create(f_cls, list(v_tc, v_cb, v_rd, v_altc), force_v4 = TRUE)
  ncdf4::ncvar_put(nc, v_tc, tc)
  ncdf4::ncvar_put(nc, v_cb, c(1150, 1150, NA, NA, NA))
  ncdf4::ncvar_put(nc, v_rd, c(0L, 0L, 0L, 1L, 0L))
  ncdf4::ncvar_put(nc, v_altc, rep(150, 5))
  ncdf4::ncatt_put(nc, "height", "standard_name", "height_above_mean_sea_level")
  ncdf4::ncatt_put(nc, 0, "file_uuid", "cls0001-test")
  ncdf4::nc_close(nc)
}

# ---- ONEFlux BIF (spec 08): synthetic GRP_ONEFLUX, GRP_FLUX_REF (CUT, VUT 2020/2021, a DD row), GRP_UST_THR ----
if ("oneflux" %in% parts) {
  of_dir <- file.path(fx, "oneflux")
  dir.create(of_dir, recursive = TRUE, showWarnings = FALSE)
  row <- function(g, grp, v, val) data.frame(SITE_ID = site, GROUP_ID = g, VARIABLE_GROUP = grp, VARIABLE = v, DATAVALUE = val)
  ref <- function(g, ver, thr, pct = NA, year = NA, res = "HH") {
    r <- rbind(row(g, "GRP_FLUX_REF", "REF_TIME_RESOLUTION", res),
               if (!is.na(pct)) row(g, "GRP_FLUX_REF", "REF_USTAR_PERCENTILE", pct),
               row(g, "GRP_FLUX_REF", "REF_USTAR_THRESHOLD", thr), row(g, "GRP_FLUX_REF", "REF_VARIABLE", "NEE"),
               row(g, "GRP_FLUX_REF", "REF_VERSION", ver))
    if (!is.na(year)) r <- rbind(r, row(g, "GRP_FLUX_REF", "REF_YEAR", year))
    r
  }
  b <- rbind(
    row(1, "GRP_ONEFLUX", "PRODUCT_FIRST_YEAR", "2020"), row(1, "GRP_ONEFLUX", "PRODUCT_LAST_YEAR", "2021"),
    row(1, "GRP_ONEFLUX", "PRODUCT_NAME", paste0("ICOSETC_", site, "_FLUXNET_2020-2021_v1.3_r1.zip")),
    ref(10, "CUT", "0.25", "48.75"), ref(11, "VUT", "0.2", "50", 2020), ref(12, "VUT", "0.3", "47.5", 2021),
    ref(13, "VUT", "0.9", NA, 2020, "DD"),
    row(20, "GRP_UST_THR", "UST_THR_YEAR", "2020"), row(20, "GRP_UST_THR", "UST_THR_VALUE", "0.21")
  )
  utils::write.csv(b, file.path(of_dir, paste0("ICOSETC_", site, "_FLUXNET_BIF_2020-2021_v1.3_r1.csv")),
                   row.names = FALSE, quote = FALSE)
}

if ("reobs" %in% parts) {
  # SIRTA 1D ReOBS structure (spec 10): hourly, time = minutes since 1970-01-01 (UTC, start of the hour),
  # time_bnds [t, t + 60 min), values float with NaN fill, counts int with -999 fill. 6 hours from 2022-06-01.
  # CHM15k (nominal 60): nb 60 (qc 0), 40 (qc 1), NaN value with nb 20, 29 (< 50 %, dropped), missing, 45 (qc 0).
  # CL31 (nominal 4): nb 4 (qc 0), 2 (qc 1), missing, 1 (dropped), 3 (qc 0), missing.
  ro_dir <- file.path(fx, "reobs")
  dir.create(ro_dir, recursive = TRUE, showWarnings = FALSE)
  tmin <- as.numeric(as.POSIXct("2022-06-01", tz = "UTC")) / 60 + 60 * (0:5)
  dim_t <- ncdf4::ncdim_def("time", "minutes since 1970-01-01", tmin, unlim = FALSE)
  dim_b <- ncdf4::ncdim_def("nbnds", "", 1:2, create_dimvar = FALSE)
  fv <- function(name, units) ncdf4::ncvar_def(name, units, list(dim_t), missval = NaN, prec = "float")
  iv <- function(name) ncdf4::ncvar_def(name, "1", list(dim_t), missval = -999L, prec = "integer")
  vals <- list(
    mlh_chm15k = c(500, 600, NaN, 700, NaN, 800), mlh_chm15k_std = c(20, 30, NaN, 40, NaN, 50),
    mlh_chm15k_nb = c(60L, 40L, 20L, 29L, -999L, 45L),
    ablh_chm15k = c(1500, 1400, NaN, NaN, NaN, 1600), ablh_chm15k_std = rep(NaN, 6), ablh_chm15k_nb = rep(-999L, 6),
    cbh_chm15k = c(NaN, 2000, NaN, NaN, NaN, 900), cbh_chm15k_std = rep(NaN, 6), cbh_chm15k_nb = rep(-999L, 6),
    mlh_cl31 = c(450, 550, NaN, 650, 300, NaN), mlh_cl31_std = c(10, 15, NaN, NaN, 5, NaN),
    mlh_cl31_nb = c(4L, 2L, -999L, 1L, 3L, -999L),
    cbh_cl31 = c(NaN, 1900, NaN, NaN, NaN, NaN), cbh_cl31_std = rep(NaN, 6), cbh_cl31_nb = rep(-999L, 6),
    tcf_cl31 = c(0, 4, NaN, 8, 1, NaN), tcf_cl31_std = rep(NaN, 6), tcf_cl31_nb = rep(-999L, 6)
  )
  vdef <- lapply(names(vals), function(k) if (is.integer(vals[[k]])) iv(k) else fv(k, if (grepl("^tcf", k)) "1" else "m"))
  names(vdef) <- names(vals)
  v_tb <- ncdf4::ncvar_def("time_bnds", "", list(dim_b, dim_t), prec = "integer")
  v_alt <- ncdf4::ncvar_def("station_alt", "m", list(), prec = "double")
  v_lat <- ncdf4::ncvar_def("station_lat", "degrees_north", list(), prec = "double")
  v_lon <- ncdf4::ncvar_def("station_lon", "degrees_east", list(), prec = "double")
  tmp_dir <- file.path(tempdir(), "reobs_fixture")
  dir.create(tmp_dir, showWarnings = FALSE)
  f_nc <- file.path(tmp_dir, "reobs_sirta_1D_request_test.nc")
  unlink(f_nc)
  nc <- ncdf4::nc_create(f_nc, c(unname(vdef), list(v_tb, v_alt, v_lat, v_lon)), force_v4 = TRUE)
  for (k in names(vals)) ncdf4::ncvar_put(nc, vdef[[k]], vals[[k]])
  ncdf4::ncvar_put(nc, v_tb, as.integer(rbind(tmin, tmin + 60)))
  ncdf4::ncvar_put(nc, v_alt, 156)
  ncdf4::ncvar_put(nc, v_lat, 48.718)
  ncdf4::ncvar_put(nc, v_lon, 2.207)
  ncdf4::ncatt_put(nc, "time", "bounds", "time_bnds")
  ga <- list(title = "SIRTA 1D ReOBS dataset (synthetic test)", license = "Inspired by CC-BY",
             acknowledgement = "Test acknowledgement.", references = "https://doi.org/10.5194/essd-10-919-2018",
             metadata_link = "https://reobs.aeris-data.fr/", reobs_version = "4.0.0", product_version = "2026.08.1",
             uuid = "00000000-test-reobs")
  for (k in names(ga)) ncdf4::ncatt_put(nc, 0, k, ga[[k]])
  ncdf4::nc_close(nc)
  zipfile <- file.path(ro_dir, "reobs_test.zip")
  unlink(zipfile)
  old <- setwd(tempdir())
  utils::zip(zipfile, file.path("reobs_fixture", basename(f_nc)), flags = "-9Xq")
  setwd(old)
}

message("Fixtures written to ", fx, " (parts: ", paste(parts, collapse = ", "), ")")
