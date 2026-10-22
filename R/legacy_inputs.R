# Inputs of the operational ETC call of the legacy calc_PBL_height(), rebuilt
# from the published ICOS L2 files only (spec 03 section 1 and Amendment 1).
# EddyPro ambient variables (not published in L2) follow the open-source
# EddyPro engine (github.com/LI-COR/eddypro-engine, commit 3cabe637ca,
# 2022-04-07): src/src_common/read_metadata_file.f90 (barometric pressure),
# src/src_rp/air_and_cell_parameters.f90 (AirAndCellParameters: ambient T and
# P, molar volume), src/src_rp/flux_params.f90 (FluxParams: water vapour,
# densities, humidity-corrected air temperature, heat capacities),
# src/src_rp/write_out_full.f90 (full-output air_temperature, air_pressure,
# air_density, air_heat_capacity), constants in
# src/src_common/m_common_global_var.f90.

#' EddyPro constants used for the ambient variables
#'
#' "EddyPro values" (src/src_common/m_common_global_var.f90): `Ru` universal
#' gas constant 8.314 J mol-1 K-1, `Rd` dry air 287.04 J kg-1 K-1, `Rw` water
#' vapour 461.5 J kg-1 K-1, `M_h2o` 18.02e-3 kg mol-1, `alpha` 0.51 (sonic
#' temperature humidity factor, flux_params.f90), `RH_max` 130 % (largest RH
#' kept). Used only for the rebuild of the EddyPro ambient variables.
#'
#' @return Named list.
eddypro_constants <- function() {
  list(Ru = 8.314, Rd = 287.04, Rw = 461.5, M_h2o = 18.02e-3, alpha = 0.51, RH_max = 130, error = -9999)
}

#' Rebuild the EddyPro ambient variables without biomet data
#'
#' Reproduces EddyPro's ambient air temperature, pressure, density and heat
#' capacity for a station without slow-response (biomet) temperature,
#' humidity and pressure (ETC setting, G 2026-10-17), from the sonic
#' temperature and the H2O mole fraction. Steps (engine routines in the file
#' header):
#' 1. `AP = 101300 exp(-altitude / 8200)` [Pa] (read_metadata_file.f90;
#'    Campbell & Norman 1998); `altitude` = site altitude only.
#' 2. `T = Ts`; molar volume `Va = Ru T / AP`; saturation vapour pressure
#'    `es = exp(77.345 + 0.0057 T - 7235 / T) / T^8.2` [Pa]; water vapour
#'    density `rho_w = (chi / Va) M_h2o 1e-3` (chi in mmol mol-1); `e = rho_w
#'    Rw T`; `RH = 100 e / es` (set to the error value outside 0..RH_max, to
#'    100 between 100 and RH_max); dry air `p_d = AP - e`, `rho_d = p_d / (Rd
#'    T)`; `rho_a = rho_d + rho_w`; `Q = rho_w / rho_a`.
#' 3. `Ta = Ts / (1 + alpha Q)` (humidity correction of the sonic temperature,
#'    van Dijk et al. 2004), then one iteration of step 2 with `T = Ta` (rho_w
#'    and Va unchanged) and `Ta = Ts / (1 + alpha Q)` again.
#' 4. `Cpd = 1005 + (Ta_pre - 273.15 + 23.12)^2 / 3364` (Ta_pre = temperature
#'    of the last density step), `Cpv = 1859 + 0.13 RH + (0.193 + 5.6e-3 RH)
#'    (Ta - 273.15) + (1e-3 + 5e-5 RH) (Ta - 273.15)^2`; `rho_cp = Cpv rho_w +
#'    Cpd rho_d`; heat capacity `cp = rho_cp / rho_a`.
#' Without H2O (or chi <= 0): `AT = Ts`, `rho = AP / (Rd Ts)`, `cp = Cpd(Ts)`
#' (the engine's fallback branches). Quirk reproduced: an RH outside 0..130 %
#' enters Cpv with the engine's error value (-9999).
#'
#' @param t_sonic_K Sonic temperature [K] (L2 `T_SONIC`).
#' @param h2o_mmol_mol Water vapour mole fraction [mmol mol-1] (L2 `H2O`).
#' @param altitude_m Site altitude [m a.s.l.] (scalar).
#' @param ep Output of [eddypro_constants()].
#' @return Tibble: `AT` air temperature [K], `AP` air pressure [Pa], `rho`
#'   moist air density [kg m-3], `cp` moist air heat capacity [J kg-1 K-1],
#'   `RH` relative humidity of the engine's iteration [%; -9999 (engine error
#'   value) outside 0..130 %, `NA` without H2O].
eddypro_ambient <- function(t_sonic_K, h2o_mmol_mol, altitude_m, ep = eddypro_constants()) {
  if (length(altitude_m) != 1 || is.na(altitude_m)) cli::cli_abort("{.arg altitude_m} must be one value [m].")
  n <- length(t_sonic_K)
  ts <- t_sonic_K
  ap <- rep(1e3 * 101.3 * exp(-altitude_m / 8200), n) # read_metadata_file.f90, l. 121
  es_fun <- function(t) exp(77.345 + 0.0057 * t - 7235 / t) / t^8.2
  cpd_fun <- function(t) 1005 + (t - 273.15 + 23.12)^2 / 3364
  rh_fun <- function(e, es) {
    rh <- e * 100 / es
    rh[rh < 0 | rh > ep$RH_max] <- ep$error
    rh[rh > 100 & rh < ep$RH_max] <- 100
    rh
  }
  va <- ep$Ru * ts / ap # air_and_cell_parameters.f90: Va = Ru T / Pr
  chi <- h2o_mmol_mol
  has_w <- !is.na(chi) & chi > 0 & !is.na(ts) & ts > 0
  rho_w <- ifelse(has_w, (chi / va) * ep$M_h2o * 1e-3, NA_real_)
  # first pass at T = Ts
  e <- rho_w * ep$Rw * ts
  p_d <- ifelse(!is.na(e) & e > 0, ap - e, ap)
  rho_d <- p_d / (ep$Rd * ts)
  rho_a <- rho_d + rho_w
  q <- rho_w / rho_a
  ta <- ts / (1 + ep$alpha * q)
  # one iteration at T = Ta (flux_params.f90, "Iterate the calculation of main quantities")
  t_pre <- ta
  e <- rho_w * ep$Rw * t_pre
  rh <- rh_fun(e, es_fun(t_pre))
  p_d <- ap - e
  rho_d <- p_d / (ep$Rd * t_pre)
  rho_a <- rho_d + rho_w
  cpd <- cpd_fun(t_pre)
  q <- rho_w / rho_a
  ta <- ts / (1 + ep$alpha * q)
  tc <- ta - 273.15
  cpv <- 1859 + 0.13 * rh + (0.193 + 5.6e-3 * rh) * tc + (1e-3 + 5e-5 * rh) * tc^2
  rho_cp <- cpv * rho_w + cpd * rho_d
  cp <- rho_cp / rho_a
  # no water vapour: engine fallback (Ta = Ts, dry air)
  nw <- !has_w & !is.na(ts) & ts > 0
  ta[nw] <- ts[nw]
  rho_a[nw] <- ap[nw] / (ep$Rd * ts[nw])
  cp[nw] <- cpd_fun(ts[nw])
  rh[nw] <- NA_real_
  dplyr::tibble(AT = ta, AP = ap, rho = rho_a, cp = cp, RH = rh)
}

#' Last observation carried forward
#'
#' `data.table::nafill(type = "locf")` semantics in base R: each `NA` takes the
#' last preceding non-missing value; leading `NA` stay `NA`.
#'
#' @param x Vector [any unit].
#' @return `x` with `NA` filled.
locf <- function(x) {
  ok <- !is.na(x)
  idx <- cumsum(ok)
  idx[idx == 0] <- NA # leading NA stay NA (a 0 index would drop elements)
  x[ok][idx]
}

#' Inputs of the operational legacy call, per processing run
#'
#' Builds the arguments of the ETC call of `calc_PBL_height()` (spec 03 section
#' 1, G 2026-10-17) from [read_icos_site()] output: `TIMESTAMP` =
#' `TIMESTAMP_END` (local standard time) as `"%Y-%m-%d %H:%M:%S"`; `MOL` = L2
#' `MO_LENGTH` with last-observation-carried-forward within each run
#' ([locf()]); `ustar`, `H` = L2 `USTAR`, `H`; `AT`, `AP`, `rho`, `cp` rebuilt
#' with [eddypro_ambient()] from L2 `T_SONIC` and `H2O`; `zm` = `zm_m`
#' (time-varying), `d` = 2/3 of `canopy_height_m` (they only enter `zdl` in
#' the legacy code, which is computed but never used). Runs (Amendment 1):
#' split at `TIMESTAMP_START` = `split_local` (local standard time); each run
#' starts at its first FLUXES row with data (any of `MO_LENGTH`, `USTAR`, `H`,
#' `T_SONIC` non-missing) and ends at the last row of the run.
#'
#' @param d Output of [read_icos_site()] (with `H2O`).
#' @param site_meta One row of [read_sites()] (`altitude_m`).
#' @param split_local Start of the second run, `YYYYMMDDHHMM` local standard
#'   time [character]; `NULL` = one run.
#' @param ambient `"rebuilt"` (default, [eddypro_ambient()]) or a sensitivity
#'   variant: `"AT_sonic"` (`AT = T_SONIC`), `"AP_altitude"` (`AP = NULL`: the
#'   legacy computes it from altitude + zm), `"rho_cp_legacy"` (`rho`, `cp` =
#'   `NULL`: computed inside the legacy function), `"amb_na"` (hypothesis test,
#'   spec 03: `AT`, `rho`, `cp` set to `NA` where L2 `H2O` is missing or the
#'   engine's RH is outside 0..130 %).
#' @return Tibble: `run_id` [integer], `timestamp_start_utc`,
#'   `timestamp_end_utc` [POSIXct UTC], `TIMESTAMP` [character], `MOL` [m],
#'   `ustar` [m s-1], `H` [W m-2], `AT` [K], `AP` [Pa], `rho` [kg m-3], `cp`
#'   [J kg-1 K-1] (the last three `NA` columns when the variant passes `NULL`,
#'   see attribute `pass_null`), `zm` [m], `d` [m].
legacy_inputs <- function(d, site_meta, split_local = "202601010000",
                          ambient = c("rebuilt", "AT_sonic", "AP_altitude", "rho_cp_legacy", "amb_na")) {
  ambient <- match.arg(ambient)
  for (k in c("TIMESTAMP_START", "TIMESTAMP_END", "MO_LENGTH", "USTAR", "H", "T_SONIC_K", "H2O")) {
    if (!k %in% names(d)) cli::cli_abort("Column {.field {k}} not found in {.arg d}.")
  }
  d <- d[order(d$timestamp_start_utc), ]
  run_id <- if (is.null(split_local)) rep(1L, nrow(d)) else ifelse(d$TIMESTAMP_START < split_local, 1L, 2L)
  has_data <- !is.na(d$MO_LENGTH) | !is.na(d$USTAR) | !is.na(d$H) | !is.na(d$T_SONIC_K)
  keep <- rep(FALSE, nrow(d))
  for (r in unique(run_id)) {
    i <- which(run_id == r)
    first <- i[which(has_data[i])[1]]
    if (!is.na(first)) keep[i[i >= first]] <- TRUE
  }
  d <- d[keep, ]
  run_id <- run_id[keep]
  amb <- eddypro_ambient(d$T_SONIC_K, d$H2O, site_meta$altitude_m)
  if (ambient == "AT_sonic") amb$AT <- d$T_SONIC_K
  if (ambient == "amb_na") {
    bad <- is.na(d$H2O) | d$H2O <= 0 | (!is.na(amb$RH) & amb$RH == -9999)
    amb$AT[bad] <- NA_real_
    amb$rho[bad] <- NA_real_
    amb$cp[bad] <- NA_real_
  }
  t_end <- as.POSIXct(d$TIMESTAMP_END, format = "%Y%m%d%H%M", tz = "UTC") # clock time, local standard
  mol <- d$MO_LENGTH
  for (r in unique(run_id)) mol[run_id == r] <- locf(mol[run_id == r])
  canopy <- if ("canopy_height_m" %in% names(d)) d$canopy_height_m else rep(NA_real_, nrow(d))
  out <- dplyr::tibble(
    run_id = run_id,
    timestamp_start_utc = d$timestamp_start_utc,
    timestamp_end_utc = d$timestamp_end_utc,
    TIMESTAMP = format(t_end, "%Y-%m-%d %H:%M:%S"),
    MOL = mol,
    ustar = d$USTAR,
    H = d$H,
    AT = amb$AT,
    AP = amb$AP,
    rho = amb$rho,
    cp = amb$cp,
    zm = if ("zm_m" %in% names(d)) d$zm_m else rep(NA_real_, nrow(d)),
    d = 2 / 3 * canopy
  )
  attr(out, "pass_null") <- switch(ambient,
    AP_altitude = "AP", rho_cp_legacy = c("rho", "cp"), character(0)
  )
  attr(out, "ambient") <- ambient
  out
}
