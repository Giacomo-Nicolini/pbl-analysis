# Rebuild the synthetic ICOS L2 and ERA5 test fixtures (tiny, no real data).
# Run from the project root: Rscript tests/testthat/fixtures/make_fixtures.R [icos] [era5]
# (no argument = both parts; rebuilding the ICOS zips changes their bytes even
# when the content is unchanged, so rebuild only the part you need).
# Requires a `zip` program for utils::zip() (on Windows: Rtools).
# Synthetic site XX-Tst, UTC+1.

parts <- commandArgs(trailingOnly = TRUE)
if (length(parts) == 0) parts <- c("icos", "era5")

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

message("Fixtures written to ", fx, " (parts: ", paste(parts, collapse = ", "), ")")
