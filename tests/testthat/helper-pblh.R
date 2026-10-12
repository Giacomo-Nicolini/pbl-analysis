# Synthetic inputs for the PBLH model tests (structure of read_icos_site() output).

test_site_meta <- function(lat = 51.31, lon = 4.52, altitude_m = 16, utc_offset_h = 1) {
  dplyr::tibble(site_id = "XX-Tst", lat = lat, lon = lon, altitude_m = altitude_m,
                utc_offset_h = utc_offset_h, zm_m = 30)
}

# Minimal read_icos_site()-like table: n half-hours from `start` with constant values.
make_icos_like <- function(n, start = "2024-06-01 00:00:00", ...) {
  t <- as.POSIXct(start, tz = "UTC") + (seq_len(n) - 1) * 1800
  d <- dplyr::tibble(
    site_id = "XX-Tst",
    timestamp_start_utc = t,
    timestamp_end_utc = t + 1800,
    H = 50, USTAR = 0.3, MO_LENGTH = -100, T_SONIC_K = 291,
    fn_H_F_MDS = 50, fn_H_F_MDS_QC = 0, fn_H_CORR = 60,
    fn_TA_F_K = 290, fn_PA_F_PA = 100000,
    met_TA_K = 289, met_PA_PA = 99000,
    zm_m = 30, excl_flag = FALSE
  )
  over <- list(...)
  for (nm in names(over)) d[[nm]] <- over[[nm]]
  d
}

# Synthetic multi-day series with a diurnal cycle: H follows the solar elevation
# (H = h_max * sin(elev) by day, -20 W m-2 at night), constant u*, L consistent
# with H: L = -u*^3 T rho cp / (k g H) (definition of the Obukhov length).
make_diurnal <- function(n_days = 3, start = "2024-06-01 00:00:00", site = test_site_meta(),
                         h_max = 300, ustar = 0.3, T_K = 293, p = 100000) {
  k <- pblh_constants()
  n <- n_days * 48
  t <- as.POSIXct(start, tz = "UTC") + (seq_len(n) - 1) * 1800
  elev <- solar_elevation(t + 900, site$lat, site$lon)
  H <- ifelse(elev > 0, h_max * sin(elev * pi / 180), -20)
  rho <- p / (k$r_d * T_K)
  L <- -ustar^3 * T_K * rho * k$cp / (k$k * k$g * H)
  make_icos_like(n, start = start, H = H, fn_H_F_MDS = H, USTAR = ustar, MO_LENGTH = L,
                 fn_TA_F_K = T_K, fn_PA_F_PA = p)
}
