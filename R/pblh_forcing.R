# Forcing of PBLHmod from read_icos_site() output.
# Spec: docs/specs/02_pblh_model.md section 3 and Amendments 3-4.

# Column of `data`, or NA vector if absent (e.g. no FLUXNET file at a site).
col_or_na <- function(data, name) {
  if (name %in% names(data)) as.numeric(data[[name]]) else rep(NA_real_, nrow(data))
}

#' Linear interpolation over short gaps
#'
#' Fills interior runs of `NA` of at most `max_gap` values that have valid
#' values on both sides. Runs containing a `block` row (e.g. excluded records)
#' are never filled. No extrapolation at the series ends.
#'
#' @param x Numeric vector on a regular time grid [any unit].
#' @param max_gap Longest run filled [number of values].
#' @param block Logical vector; `TRUE` rows are never filled nor used to bridge.
#' @return List: `x` filled vector, `filled` logical (`TRUE` where filled).
interp_short_gaps <- function(x, max_gap, block = rep(FALSE, length(x))) {
  filled <- rep(FALSE, length(x))
  if (max_gap < 1 || length(x) < 3) return(list(x = x, filled = filled))
  r <- rle(is.na(x))
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  for (j in which(r$values & r$lengths <= max_gap)) {
    i0 <- starts[j] - 1
    i1 <- ends[j] + 1
    if (i0 < 1 || i1 > length(x)) next
    idx <- starts[j]:ends[j]
    if (any(block[c(i0, idx, i1)])) next
    w <- (idx - i0) / (i1 - i0)
    x[idx] <- x[i0] + w * (x[i1] - x[i0])
    filled[idx] <- TRUE
  }
  list(x = x, filled = filled)
}

#' Implausible friction velocity (input QC)
#'
#' Two tests (spec 05; spec 02 Amendment 12), applied to the measured u*
#' before any gap filling:
#' 1. ratio test: `ustar / ws > opts$ustar_ws_max` where `ws >=
#'    opts$ws_min_ratio_test` (not applied where `ws` is missing or low, where
#'    the ratio is unstable). Over a rough surface in neutral conditions the log
#'    law gives u*/U = k / ln((z - d) / z0), about 0.14-0.19 at a forest site
#'    such as BE-Bra (z - d ~ 24 m, z0 ~ 1-3 m); instability raises the ratio,
#'    but values above ~0.5 are not physical for a tower sonic (surface-layer
#'    similarity, e.g. Panofsky & Dutton 1984; Kaimal & Finnigan 1994);
#' 2. absolute cap: `ustar > opts$ustar_max` (u* of a few m s-1 over land
#'    requires storm-force winds).
#' Implausible u* also invalidates L (L = -u*^3 T_v / (k g w'theta_v')).
#'
#' @param ustar Friction velocity [m s-1].
#' @param ws Wind speed of the same sonic (FLUXES `WS`) [m s-1].
#' @param opts Output of [pblh_options()] (`ustar_ws_max` [-],
#'   `ws_min_ratio_test` [m s-1], `ustar_max` [m s-1]; `Inf` disables a test).
#' @return List of two logical vectors: `ratio` (rejected by the ratio test),
#'   `cap` (rejected by the cap); a row is implausible if either is `TRUE`.
ustar_implausible <- function(ustar, ws, opts = pblh_options()) {
  ok <- !is.na(ustar)
  ratio <- ok & !is.na(ws) & ws >= opts$ws_min_ratio_test & ustar / ws > opts$ustar_ws_max
  cap <- ok & ustar > opts$ustar_max
  list(ratio = ratio, cap = cap)
}

#' Obukhov length from the forcing
#'
#' `L = - ustar^3 T / (k g wtheta)` (Obukhov 1946; Stull 1988, Sect. 5.7,
#' Eq. 5.7c, with the air temperature in place of the virtual potential
#' temperature and the kinematic heat flux in place of the buoyancy flux: the
#' moisture contribution is neglected).
#'
#' @param ustar Friction velocity [m s-1].
#' @param wtheta Kinematic heat flux [K m s-1].
#' @param T_K Air temperature [K].
#' @param const Constants from [pblh_constants()] (`k`, `g`).
#' @return Obukhov length [m]; `NA` where not finite (e.g. `wtheta = 0`).
obukhov_length <- function(ustar, wtheta, T_K, const = pblh_constants()) {
  L <- -ustar^3 * T_K / (const$k * const$g * wtheta)
  L[!is.finite(L)] <- NA_real_
  L
}

#' Consistency of L2 Obukhov length with the forcing
#'
#' Spec 06 (spec 02 Amendment 13). A row is checked if `L`, `L_forcing` and
#' `wtheta0` are available and `|wtheta0| > opts$l_check_wtheta_min` (near
#' neutral, |L| is large and ill-determined in both estimates). A checked row
#' is inconsistent if the signs differ (`"sign"`) or, with the same sign,
#' `max(L / L_forcing, L_forcing / L) > opts$l_consistency_factor`
#' (`"magnitude"`).
#'
#' @param L L2 Obukhov length [m].
#' @param L_forcing Obukhov length from the forcing ([obukhov_length()]) [m].
#' @param wtheta0 Kinematic heat flux of the same row [K m s-1].
#' @param opts Output of [pblh_options()] (`l_consistency_factor` [-],
#'   `l_check_wtheta_min` [K m s-1]).
#' @return Character: `"not_checked"`, `"consistent"`, `"sign"`, `"magnitude"`.
l_consistency <- function(L, L_forcing, wtheta0, opts = pblh_options()) {
  checked <- !is.na(L) & !is.na(L_forcing) & !is.na(wtheta0) & abs(wtheta0) > opts$l_check_wtheta_min
  sign_bad <- checked & sign(L) != sign(L_forcing)
  ratio <- pmax(L / L_forcing, L_forcing / L)
  mag_bad <- checked & !sign_bad & ratio > opts$l_consistency_factor
  ifelse(!checked, "not_checked", ifelse(sign_bad, "sign", ifelse(mag_bad, "magnitude", "consistent")))
}

#' Prepare the PBLH model forcing
#'
#' Selects and combines the forcing variables (with fallbacks), derives air
#' density, kinematic heat flux and Coriolis parameter, sets excluded rows to
#' `NA` and fills short gaps by linear interpolation (L via 1/L).
#'
#' Rules (spec 02 section 3, Amendments 3-4, 12, 13). Order: excluded rows ->
#' u* filter -> L consistency check -> short-gap interpolation -> `rho`,
#' `wtheta`.
#' - `ustar` and `L` are set to `NA` where u* is implausible
#'   ([ustar_implausible()], FLUXES `WS`; `flag_ustar_implausible`), before
#'   the short-gap interpolation below; `H` and the other inputs are untouched.
#' - L consistency (Amendment 13): on the not-yet-interpolated values,
#'   `wtheta0 = H / (rho0 cp)`, `rho0 = p_Pa / (r_d T_K)`, `L_forcing` from
#'   [obukhov_length()], `L_check` from [l_consistency()] (computed for every
#'   `opts$l_source`). `"L2_checked"`: rows with `L_check` `"sign"` or
#'   `"magnitude"` get `L = L_forcing` (`flag_L_replaced`); `"L2"`: L2 L as is;
#'   `"forcing"`: `L = L_forcing` wherever available. Missing L2 L is never
#'   filled from the forcing (only by the short-gap interpolation).
#' - `H` by `opts$h_source`. `"H_F_MDS"`: `fn_H_F_MDS`, set to `NA` where
#'   `fn_H_F_MDS_QC > opts$h_qc_max` (`flag_H_qc_masked`); where `fn_H_F_MDS` is
#'   missing in the file, FLUXES `H` if `opts$h_fallback`. QC-masked rows stay `NA`.
#' - `T_K`: `fn_TA_F_K` -> `met_TA_K` -> `T_SONIC_K`.
#' - `p_Pa`: `fn_PA_F_PA` -> `met_PA_PA` -> `p0 exp(-altitude_m / scale_height)`
#'   (isothermal barometric formula, Wallace & Hobbs 2006).
#' - Gaps of at most `opts$max_gap_interp` half-hours in `ustar`, `1/L`, `H`,
#'   `T_K`, `p_Pa` are interpolated (never across excluded rows); `rho` and
#'   `wtheta` are derived after interpolation.
#'
#' @param data Output of [read_icos_site()] (30-min grid, sorted).
#' @param site_meta One row of [read_sites()] (`lat`, `lon`; `altitude_m` only
#'   needed for the barometric fallback).
#' @param opts Output of [pblh_options()].
#' @return Tibble, one row per input row: `timestamp_start_utc`,
#'   `timestamp_end_utc` [POSIXct UTC]; `ustar` [m s-1]; `L` Obukhov length [m];
#'   `H` sensible heat flux [W m-2]; `H_source`; `T_K` air temperature [K];
#'   `T_source`; `p_Pa` air pressure [Pa]; `p_source`; `rho` dry-air density
#'   [kg m-3]; `wtheta` kinematic heat flux [K m s-1]; `f_cor` Coriolis
#'   parameter [s-1, signed]; `is_day` [logical]; `excl_flag` [logical];
#'   `WS` wind speed (FLUXES) [m s-1]; `L_l2` L2 `MO_LENGTH` after the u*
#'   filter [m]; `L_forcing` [m]; `L_source` (`"L2"`, `"forcing"`,
#'   `"interpolated"`, `NA`); `L_check` (`"not_checked"`, `"consistent"`,
#'   `"sign"`, `"magnitude"`; `NA` on excluded rows); `flag_L_replaced`,
#'   `flag_ustar_implausible`, `flag_ustar_interp`, `flag_L_interp`,
#'   `flag_H_interp`, `flag_H_qc_masked` [logical].
prepare_pblh_forcing <- function(data, site_meta, opts = pblh_options()) {
  const <- pblh_constants()
  if (nrow(site_meta) != 1 || is.na(site_meta$lat) || is.na(site_meta$lon)) {
    cli::cli_abort("{.arg site_meta} must be one row with non-missing {.field lat} and {.field lon}.")
  }
  n <- nrow(data)
  excl <- if ("excl_flag" %in% names(data)) dplyr::coalesce(data$excl_flag, FALSE) else rep(FALSE, n)

  # --- sensible heat flux [W m-2] ---
  h_fx <- col_or_na(data, "H")
  qc_masked <- rep(FALSE, n)
  if (opts$h_source == "H_F_MDS") {
    h_mds <- col_or_na(data, "fn_H_F_MDS")
    qc <- col_or_na(data, "fn_H_F_MDS_QC")
    missing_in_file <- is.na(h_mds)
    qc_masked <- !missing_in_file & !is.na(qc) & qc > opts$h_qc_max
    H <- ifelse(qc_masked, NA_real_, h_mds)
    H_source <- ifelse(is.na(H), NA_character_, "H_F_MDS")
    if (opts$h_fallback) {
      use_fx <- missing_in_file & !is.na(h_fx)
      H[use_fx] <- h_fx[use_fx]
      H_source[use_fx] <- "H"
    }
  } else if (opts$h_source == "H") {
    H <- h_fx
    H_source <- ifelse(is.na(H), NA_character_, "H")
  } else {
    H <- col_or_na(data, "fn_H_CORR")
    H_source <- ifelse(is.na(H), NA_character_, "H_CORR")
  }

  # --- air temperature [K] ---
  t_cand <- list(TA_F = "fn_TA_F_K", TA_meteo = "met_TA_K", T_SONIC = "T_SONIC_K")
  T_K <- rep(NA_real_, n)
  T_source <- rep(NA_character_, n)
  for (s in names(t_cand)) {
    v <- col_or_na(data, t_cand[[s]])
    use <- is.na(T_K) & !is.na(v)
    T_K[use] <- v[use]
    T_source[use] <- s
  }

  # --- air pressure [Pa] ---
  p_cand <- list(PA_F = "fn_PA_F_PA", PA_meteo = "met_PA_PA")
  p_Pa <- rep(NA_real_, n)
  p_source <- rep(NA_character_, n)
  for (s in names(p_cand)) {
    v <- col_or_na(data, p_cand[[s]])
    use <- is.na(p_Pa) & !is.na(v)
    p_Pa[use] <- v[use]
    p_source[use] <- s
  }
  need_baro <- is.na(p_Pa) & !excl
  if (any(need_baro)) {
    if (is.null(site_meta$altitude_m) || is.na(site_meta$altitude_m)) {
      cli::cli_abort("{sum(need_baro)} rows need the barometric pressure fallback but {.field altitude_m} is missing for the site.")
    }
    # isothermal barometric formula p = p0 exp(-z / H_s); Wallace & Hobbs 2006, Sect. 3.2
    p_Pa[need_baro] <- const$p0 * exp(-site_meta$altitude_m / const$scale_height)
    p_source[need_baro] <- "barometric"
  }

  ustar <- col_or_na(data, "USTAR")
  L <- col_or_na(data, "MO_LENGTH")

  # --- excluded rows: no forcing ---
  ustar[excl] <- NA
  L[excl] <- NA
  H[excl] <- NA
  T_K[excl] <- NA
  p_Pa[excl] <- NA
  H_source[excl] <- NA
  T_source[excl] <- NA
  p_source[excl] <- NA
  ws <- col_or_na(data, "WS")
  ws[excl] <- NA

  # --- u* plausibility (input QC, not output clipping; spec 02 Amendment 12) ---
  rej <- ustar_implausible(ustar, ws, opts)
  bad_u <- rej$ratio | rej$cap
  ustar[bad_u] <- NA
  L[bad_u] <- NA # L ~ u*^3: unusable when u* is

  # --- Obukhov length: L2 vs forcing, before interpolation (spec 02 Amendment 13) ---
  rho0 <- p_Pa / (const$r_d * T_K) # ideal gas law, dry air [kg m-3]
  wtheta0 <- H / (rho0 * const$cp) # kinematic heat flux of the row [K m s-1]
  L_forcing <- obukhov_length(ustar, wtheta0, T_K, const)
  L_l2 <- L
  L_check <- l_consistency(L, L_forcing, wtheta0, opts)
  L_check[excl] <- NA
  repl <- switch(opts$l_source,
    L2_checked = !is.na(L_check) & L_check %in% c("sign", "magnitude"),
    L2 = rep(FALSE, n),
    forcing = !is.na(L_forcing)
  )
  if (opts$l_source == "forcing") L <- rep(NA_real_, n) # option D: L2 L not used at all
  L[repl] <- L_forcing[repl]
  L_source <- ifelse(repl, "forcing", ifelse(is.na(L), NA_character_, "L2"))

  # --- short gaps (never across excluded rows) ---
  mg <- opts$max_gap_interp
  iu <- interp_short_gaps(ustar, mg, excl)
  # L through 1/L: smooth through neutral (Amendment 3). L = 0 is invalid (not used).
  inv_l <- ifelse(is.finite(L) & L != 0, 1 / L, NA_real_)
  il <- interp_short_gaps(inv_l, mg, excl)
  fill_l <- il$filled & is.na(L)
  L[fill_l] <- 1 / il$x[fill_l] # 1/L = 0 gives L = Inf (neutral)
  L_source[fill_l] <- "interpolated"
  ih <- interp_short_gaps(H, mg, excl)
  it <- interp_short_gaps(T_K, mg, excl)
  ip <- interp_short_gaps(p_Pa, mg, excl)
  T_source[it$filled] <- "interpolated"
  p_source[ip$filled] <- "interpolated"
  H_source[ih$filled] <- "interpolated"

  T_K <- it$x
  p_Pa <- ip$x
  H <- ih$x
  rho <- p_Pa / (const$r_d * T_K)  # ideal gas law, dry air [kg m-3]
  wtheta <- H / (rho * const$cp)   # surface kinematic heat flux [K m s-1]
  f_cor <- 2 * const$omega * sin(site_meta$lat * pi / 180) # Coriolis parameter [s-1]

  dplyr::tibble(
    timestamp_start_utc = data$timestamp_start_utc,
    timestamp_end_utc = data$timestamp_end_utc,
    ustar = iu$x,
    WS = ws,
    L = L,
    L_l2 = L_l2,
    L_forcing = L_forcing,
    L_source = L_source,
    L_check = L_check,
    H = H,
    H_source = H_source,
    T_K = T_K,
    T_source = T_source,
    p_Pa = p_Pa,
    p_source = p_source,
    rho = rho,
    wtheta = wtheta,
    f_cor = f_cor,
    is_day = is_day(data$timestamp_start_utc, data$timestamp_end_utc, site_meta$lat, site_meta$lon),
    excl_flag = excl,
    flag_L_replaced = repl,
    flag_ustar_implausible = bad_u,
    flag_ustar_interp = iu$filled,
    flag_L_interp = fill_l,
    flag_H_interp = ih$filled,
    flag_H_qc_masked = qc_masked
  )
}
