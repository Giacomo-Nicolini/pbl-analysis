utc <- function(x) as.POSIXct(x, tz = "UTC")

era_root <- normalizePath(fixture_path("era5_2026"))
era_cfg <- function(root = era_root) {
  list(
    era5_sl_dir = file.path(root, "sl", "{site}"),
    era5_pl_dir = file.path(root, "pl", "{site}"),
    era5_static_dir = file.path(root, "static")
  )
}
find_q <- function(family, cfg = era_cfg(), ...) {
  suppressMessages(find_era5_files("XX-Tst", cfg, family, ...))
}

# Copy the ERA5 fixtures to a temporary folder (for tests that alter files).
copy_era5 <- function(env = parent.frame()) {
  tmp <- withr::local_tempdir(.local_envir = env)
  file.copy(era_root, tmp, recursive = TRUE)
  file.path(tmp, "era5_2026")
}
# Apply `fun` to the lines of a file and write it back.
edit_lines <- function(path, fun) writeLines(fun(readLines(path)), path)
# Apply `fun` to the second column of the data block of an ERA5 file.
edit_values <- function(path, fun) {
  edit_lines(path, function(x) {
    start <- grep("FINE METADATA", x) + 2
    i <- start:length(x)
    parts <- strsplit(x[i], ",", fixed = TRUE)
    v <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
    x[i] <- paste(vapply(parts, `[`, "", 1), ifelse(is.na(v), "", format(fun(v), digits = 15)), sep = ",")
    x
  })
}
sl_file <- function(root, var, flags = "north_+0__south_+0__east_+0__west_+0") {
  file.path(root, "sl", "XX-Tst", paste0(
    "XX-Tst__reanalysis_era5_single_level__", var, "__20200101T0000__20200102T2300__", flags, ".csv"
  ))
}

test_that("era5_variables is a complete registry", {
  v <- era5_variables()
  expect_setequal(v$variable, c(
    "boundary_layer_height", "cloud_base_height", "surface_sensible_heat_flux",
    "temperature", "geopotential", "z", "lsm", "sdor"
  ))
  expect_false(anyDuplicated(v$column) > 0)
  expect_equal(v$column[v$variable == "surface_sensible_heat_flux"], "H_era5_Wm2")
})

test_that("parse_era5_header reads keys, converts numbers and levels", {
  h <- parse_era5_header(sl_file(era_root, "boundary_layer_height"))
  expect_equal(h$point_id, "XX-Tst")
  expect_identical(h$latitude_grid, 50)
  expect_identical(h$n_rows, 48)
  expect_type(h$longitude_conversion, "character") # longitude_* but text
  expect_equal(h$n_header_lines, 25L) # title + 23 keys + end marker
  expect_equal(h$first_time, "2020-01-01 00:00:00")

  pl <- list.files(file.path(era_root, "pl", "XX-Tst"), full.names = TRUE)[1]
  expect_identical(parse_era5_header(pl)$levels, c(700L, 850L, 1000L))

  st <- file.path(era_root, "static",
                  "XX-Tst__single-levels__z__197901201200__197901201200__north_-1__south_+0__east_+1__west_+0.csv")
  hs <- parse_era5_header(st)
  expect_identical(hs$pixel_rows_from_center, -1)
  expect_identical(hs$pixel_cols_from_center, 1)
  expect_type(hs$grid_scanning, "character")
  expect_identical(hs$levels, NA_integer_)
  expect_equal(hs$units, "m**2 s**-2")
})

test_that("parse_era5_header aborts without end marker", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("# ==== METADATA ====", "# point_id: XX-Tst", "time,x", "2020-01-01 00:00:00,1"), tmp)
  expect_error(parse_era5_header(tmp), "end marker")
})

test_that("find_era5_files takes cell positions from the header (direction flags)", {
  expect_message(f <- find_era5_files("XX-Tst", era_cfg(), "single_level"), "ignoring 2 files")
  expect_equal(nrow(f), 27)
  expect_false(any(grepl("2m_temperature|README", f$file)))
  b <- f[f$variable == "boundary_layer_height", ]
  expect_equal(nrow(b), 9)
  expect_setequal(paste(b$cell_row, b$cell_col), paste(rep(-1:1, 3), rep(-1:1, each = 3)))
  ne <- b[b$cell_row == 1 & b$cell_col == 1, ]
  expect_equal(c(ne$grid_lat, ne$grid_lon), c(50.25, 10.25))
  expect_match(ne$file, "north_+1__south_+0__east_+1__west_+0", fixed = TRUE)
  tower <- b[b$cell_row == 0 & b$cell_col == 0, ]
  expect_equal(tower$grid_dist_km, 0)
  expect_equal(tower$time_last, utc("2020-01-02 23:00:00"))
  expect_equal(tower$header_n_rows, 48L)
})

test_that("find_era5_files handles signed offsets of static files; other stations ignored silently", {
  expect_no_message(f <- find_era5_files("XX-Tst", era_cfg(), "static"))
  expect_equal(nrow(f), 27)
  se <- f[f$variable == "z" & f$cell_row == -1 & f$cell_col == 1, ]
  expect_equal(c(se$grid_lat, se$grid_lon), c(49.75, 10.25))
  expect_match(se$file, "north_-1__south_+0__east_+1__west_+0", fixed = TRUE)
})

test_that("a file whose name contradicts its header aborts (both conventions)", {
  root <- copy_era5()
  # NE cell file (north_+1, east_+1) with the header of the SE cell
  edit_lines(sl_file(root, "boundary_layer_height", "north_+1__south_+0__east_+1__west_+0"),
             function(x) sub("^# latitude_grid: 50.25$", "# latitude_grid: 49.75", x))
  expect_error(find_q("single_level", era_cfg(root)), "disagree")

  root2 <- copy_era5()
  st <- file.path(root2, "static",
                  "XX-Tst__single-levels__lsm__197901201200__197901201200__north_+1__south_+0__east_+0__west_+0.csv")
  edit_lines(st, function(x) sub("^# latitude_grid: 50.25$", "# latitude_grid: 49.75", x))
  expect_error(find_q("static", era_cfg(root2)), "Name: row 1, col 0; header .* row -1, col 0")
})

test_that("missing cell warns, duplicated cell and missing tower cell abort", {
  root <- copy_era5()
  file.remove(sl_file(root, "cloud_base_height", "north_+0__south_+1__east_+0__west_+0"))
  expect_warning(f <- find_q("single_level", era_cfg(root)), "1 of 9 cells missing")
  expect_equal(sum(f$variable == "cloud_base_height"), 8)

  root2 <- copy_era5()
  src <- sl_file(root2, "boundary_layer_height", "north_+1__south_+0__east_+0__west_+0")
  file.copy(src, sub("20200101T0000", "20190101T0000", src))
  expect_error(find_q("single_level", era_cfg(root2)), "duplicated cell")

  root3 <- copy_era5()
  file.remove(sl_file(root3, "boundary_layer_height"))
  expect_error(find_q("single_level", era_cfg(root3)), "containing the\\s+tower")
})

test_that("check_era5_site_coords warns on coordinate mismatch", {
  f <- find_q("single_level", variable = "boundary_layer_height")
  sites <- read_sites(fixture_path("sites.csv"))
  expect_no_warning(expect_true(check_era5_site_coords(f, sites)))
  sites$lat[sites$site_id == "XX-Tst"] <- 50.01
  expect_warning(expect_false(check_era5_site_coords(f, sites)), "differ from site metadata")
})

test_that("read_era5_single_level reads BLH of all cells or the tower cell", {
  f <- find_q("single_level")
  b <- read_era5_single_level(f, "boundary_layer_height")
  expect_equal(nrow(b), 9 * 48)
  expect_named(b, c("site_id", "cell_row", "cell_col", "grid_lat", "grid_lon", "grid_dist_km",
                    "timestamp_utc", "blh_m"))
  expect_equal(attr(b$timestamp_utc, "tzone"), "UTC")
  bc <- read_era5_single_level(f, "boundary_layer_height", cells = "center")
  expect_equal(nrow(bc), 48)
  # fixture: 200 + 800 * solar_day + 40 row + 10 col; at 00 UTC solar_day = 0
  expect_equal(b$blh_m[b$cell_row == 1 & b$cell_col == -1 & b$timestamp_utc == utc("2020-01-01 00:00:00")], 230)
})

test_that("cloud base NA (no cloud) is kept and reported", {
  f <- find_q("single_level")
  expect_message(cb <- read_era5_single_level(f, "cloud_base_height", cells = "center"), "12 of 48")
  expect_equal(sum(is.na(cb$cbh_m)), 12)
})

test_that("sshf: -360000 J m-2 -> +100 W m-2 over (t - 1 h, t]", {
  f <- find_q("single_level")
  h <- read_era5_single_level(f, "surface_sensible_heat_flux", cells = "center")
  r <- h[h$timestamp_utc == utc("2020-01-01 12:00:00"), ]
  expect_equal(r$H_era5_Wm2, 100)
  expect_equal(r$timestamp_start_utc, utc("2020-01-01 11:00:00"))
  expect_equal(r$timestamp_end_utc, utc("2020-01-01 12:00:00"))
  pk <- check_era5_h_sign(h)
  expect_equal(pk$peak_solar_hour, 11.5)
})

test_that("flipped sshf sign fails the diurnal sign check", {
  root <- copy_era5()
  p <- sl_file(root, "surface_sensible_heat_flux", "north_+0__south_+0__east_+1__west_+0")
  edit_values(p, function(v) -v)
  f <- find_q("single_level", era_cfg(root))
  expect_error(read_era5_single_level(f, "surface_sensible_heat_flux"), "solar\\s+noon")
})

test_that("implausible BLH unit aborts", {
  root <- copy_era5()
  for (fl in c("north_+0__south_+0__east_+0__west_+0", "north_+1__south_+0__east_+0__west_+0")) {
    edit_values(sl_file(root, "boundary_layer_height", fl), function(v) v * 1000)
  }
  f <- find_q("single_level", era_cfg(root))
  expect_error(read_era5_single_level(f, "boundary_layer_height", cells = "center"), "outside \\[20, 3000\\]")
})

test_that("n_rows and first/last time mismatches abort", {
  root <- copy_era5()
  p <- sl_file(root, "boundary_layer_height")
  edit_lines(p, function(x) sub("^# n_rows: 48$", "# n_rows: 47", x))
  f <- find_q("single_level", era_cfg(root), variable = "boundary_layer_height")
  expect_error(read_era5_single_level(f, "boundary_layer_height", cells = "center"), "n_rows")

  root2 <- copy_era5()
  p2 <- sl_file(root2, "boundary_layer_height")
  edit_lines(p2, function(x) sub("^# last_time: .*$", "# last_time: 2020-01-02 22:00:00", x))
  f2 <- find_q("single_level", era_cfg(root2), variable = "boundary_layer_height")
  expect_error(read_era5_single_level(f2, "boundary_layer_height", cells = "center"), "time range differs")
})

test_that("pressure levels: long format, theta known value, zg = z / g", {
  f <- find_q("pressure_level")
  expect_equal(nrow(f), 18)
  p <- read_era5_pressure_levels(f)
  expect_equal(nrow(p), 9 * 8 * 3)
  expect_named(p, c("site_id", "cell_row", "cell_col", "timestamp_utc", "p_hPa", "T_K", "z_m2s2",
                    "zg_m", "theta_K"))
  r <- p[p$cell_row == 0 & p$cell_col == 0 & p$timestamp_utc == utc("2020-01-01 00:00:00") & p$p_hPa == 850, ]
  expect_equal(r$T_K, 280)
  # theta = 280 * (1000 / 850)^(287.05 / 1005) = 280 * 1.176471^0.285622 = 280 * 1.047496 = 293.30 K
  expect_equal(r$theta_K, 293.30, tolerance = 0.01 / 293.30)
  expect_equal(r$zg_m, 1500, tolerance = 1e-9)
  expect_equal(r$zg_m, r$z_m2s2 / 9.80665)
  expect_equal(nrow(read_era5_pressure_levels(f, cells = "center")), 8 * 3)
})

test_that("pressure levels: header levels mismatch and non-monotonic zg abort", {
  root <- copy_era5()
  pf <- list.files(file.path(root, "pl", "XX-Tst"), pattern = "temperature.*north_\\+0__south_\\+0__east_\\+0__west_\\+0",
                   full.names = TRUE)
  edit_lines(pf, function(x) sub("^# levels: .*$", "# levels: 700,850", x))
  expect_error(read_era5_pressure_levels(find_q("pressure_level", era_cfg(root))), "levels")

  root2 <- copy_era5()
  zf <- list.files(file.path(root2, "pl", "XX-Tst"), pattern = "geopotential.*north_\\+0__south_\\+0__east_\\+0__west_\\+0",
                   full.names = TRUE)
  # swap the 700 and 1000 hPa columns' values -> zg decreasing upward
  edit_lines(zf, function(x) {
    i <- (grep("FINE METADATA", x) + 2):length(x)
    s <- strsplit(x[i], ",", fixed = TRUE)
    x[i] <- vapply(s, function(v) paste(v[1], v[4], v[3], v[2], sep = ","), "")
    x
  })
  expect_error(read_era5_pressure_levels(find_q("pressure_level", era_cfg(root2))), "not increasing")
})

test_that("read_era5_static: one row per cell, orography from z / g, units checked", {
  sites <- read_sites(fixture_path("sites.csv"))
  st <- read_era5_static("XX-Tst", era_cfg(), sites)
  expect_equal(nrow(st), 9)
  expect_named(st, c("site_id", "cell_row", "cell_col", "grid_lat", "grid_lon", "grid_dist_km",
                     "orog_m", "lsm", "sdor_m"))
  tw <- st[st$cell_row == 0 & st$cell_col == 0, ]
  expect_equal(tw$orog_m, 100, tolerance = 1e-9)
  expect_equal(st$orog_m[st$cell_row == -1 & st$cell_col == 1], 100 - 10 + 5, tolerance = 1e-9)
  expect_equal(st$lsm[st$cell_row == 1 & st$cell_col == -1], 0.6)
  expect_equal(tw$sdor_m, 5)

  root <- copy_era5()
  sf <- file.path(root, "static",
                  "XX-Tst__single-levels__sdor__197901201200__197901201200__north_+0__south_+0__east_+0__west_+0.csv")
  edit_lines(sf, function(x) sub("^# units: m$", "# units: km", x))
  expect_error(read_era5_static("XX-Tst", era_cfg(root)), "units")
})
