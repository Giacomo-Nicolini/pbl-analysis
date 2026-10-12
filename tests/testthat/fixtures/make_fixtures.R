# Rebuild the synthetic ICOS L2 and ERA5 test fixtures (tiny, no real data).
# Run from the project root: Rscript tests/testthat/fixtures/make_fixtures.R
# Requires a `zip` program for utils::zip() (on Windows: Rtools).
# Synthetic site XX-Tst, UTC+1.

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

# ERA5 BLH (plain CSVs).
era_dir <- file.path(fx, "era5", site)
dir.create(era_dir, recursive = TRUE, showWarnings = FALSE)
era_file <- function(product, year, lon = "+0", lat = "+0", var = "boundary_layer_height") {
  file.path(era_dir, sprintf(
    "%s__%s__reanalysis-era5-single-levels__%s__%d__lon%s__lat%s.csv",
    site, product, var, year, lon, lat
  ))
}
era_df <- function(year, values, lon_off = 0, lat_off = 0, gap = FALSE) {
  t <- as.POSIXct(sprintf("%d-01-01 00:00:00", year), tz = "UTC") + (seq_along(values) - 1) * 3600
  if (gap) t[length(t)] <- t[length(t)] + 3600 # one missing hour before the last value
  data.frame(
    value = values,
    time = format(t, "%Y-%m-%d %H:%M:%S"),
    latitude = 50.0, longitude = 10.0,
    location = sprintf("lon%+d__lat%+d", lon_off, lat_off),
    pos_my_lat = 1, pos_my_lon = 1,
    my_lat = 50.0 + 0.25 * lat_off, my_lon = 10.0 + 0.25 * lon_off,
    datetime = format(t, "%Y-%m-%d %H:%M:%S"),
    year = year
  )
}
w <- function(df, f) utils::write.csv(df, f, row.names = FALSE, quote = FALSE)
w(era_df(2020, c(500, 520, 540, 560)), era_file("ERA5", 2020))
w(era_df(2020, c(1, 1, 1, 1)), era_file("ERA5T", 2020))            # superseded by ERA5
w(era_df(2021, c(600, -5, 620, 640), gap = TRUE), era_file("ERA5T", 2021)) # negative value, one gap
w(era_df(2020, c(700, 710, 720, 730), lon_off = 1, lat_off = -1), era_file("ERA5", 2020, "+1", "-1"))
w(era_df(2020, c(280, 281, 282, 283)), era_file("ERA5", 2020, var = "2m_temperature")) # ignored
writeLines("not an ERA5 file", file.path(era_dir, "README.txt"))                       # ignored

message("Fixtures written to ", fx)
