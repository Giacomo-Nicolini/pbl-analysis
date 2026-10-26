# Inputs, validity and subsets for the footprint sensitivity to the PBLH source (spec 08, phase 5a).
# Footprint model: Kljun et al. (2015), `Kljun2015` (FFP). Only h differs between the PBLH cases.

#' Footprint heights from the EC height and the canopy height
#'
#' ETC / EddyPro convention (`legacy/FFP_analysis_dev.R`, lines 336-353): displacement
#' height d = 2/3 hc, roughness length z0 = 0.15 hc, measurement height above the
#' displacement height zm = EC height - d (Kljun et al. 2015: zm = z_receptor - z_d).
#'
#' @param ec_height_m EC system height above ground [m].
#' @param canopy_height_m Canopy height hc [m].
#' @return Tibble: `d_m`, `z0_m`, `zm_m` [m].
footprint_heights <- function(ec_height_m, canopy_height_m) {
  d <- 2 / 3 * canopy_height_m
  dplyr::tibble(d_m = d, z0_m = 0.15 * canopy_height_m, zm_m = ec_height_m - d)
}

#' ERA5 BLH at half-hour midpoints
#'
#' Linear interpolation in time between the two hourly instantaneous values that
#' bracket the half-hour midpoint (start + 15 min); `NA` unless both exist (spec 08).
#'
#' @param t_start Half-hour start times [POSIXct UTC].
#' @param era_time ERA5 times [POSIXct UTC, instantaneous].
#' @param era_h ERA5 BLH [m].
#' @return BLH at the half-hour midpoints [m].
era5_at_halfhour_mid <- function(t_start, era_time, era_h) {
  mid <- as.numeric(t_start) + 900
  ta <- floor(mid / 3600) * 3600
  te <- as.numeric(era_time)
  ha <- era_h[match(ta, te)]
  hb <- era_h[match(ta + 3600, te)]
  ha + (mid - ta) / 3600 * (hb - ha)
}

#' PBLH cases on the half-hour grid
#'
#' @param t_start Half-hour start times [POSIXct UTC].
#' @param pblhmod Standard table `pblhmod_default`.
#' @param era5_center,era5_bilinear Standard tables `era5_blh_center`, `era5_blh_bilinear`.
#' @param measured Standard table `cloudnet_eps_mlh`, or `NULL`.
#' @param oper Standard table `icos_oper`, or `NULL` (optional case, internal use only).
#' @return Tibble with one h column per case [m]: `pblhmod`, `era5`, `era5_bilinear`,
#'   `pblhmod_x0.8`, `pblhmod_x1.2` (Kljun et al. 2015 Table 5: +-20 %), and if
#'   given `measured` (exact values only, qc 0/2), `bound` (qc 3 lower bounds as h,
#'   labelled "bound"), `oper`.
footprint_h_cases <- function(t_start, pblhmod, era5_center, era5_bilinear, measured = NULL, oper = NULL) {
  at <- function(std) std$pblh_m[match(as.numeric(t_start), as.numeric(std$timestamp_start_utc))]
  p <- at(pblhmod)
  out <- dplyr::tibble(
    pblhmod = p,
    era5 = era5_at_halfhour_mid(t_start, era5_center$timestamp_start_utc, era5_center$pblh_m),
    era5_bilinear = era5_at_halfhour_mid(t_start, era5_bilinear$timestamp_start_utc, era5_bilinear$pblh_m),
    `pblhmod_x0.8` = 0.8 * p,
    `pblhmod_x1.2` = 1.2 * p
  )
  if (!is.null(measured)) {
    i <- match(as.numeric(t_start), as.numeric(measured$timestamp_start_utc))
    q <- measured$qc_flag[i]
    out$measured <- ifelse(q %in% c(0L, 2L), measured$pblh_m[i], NA_real_)
    out$bound <- ifelse(q %in% 3L, measured$pblh_m[i], NA_real_)
  }
  if (!is.null(oper)) out$oper <- at(oper)
  out
}

#' ETC roughness-sublayer multiplier
#'
#' `legacy/FFP_analysis_dev.R` (lines 1046-1066): starting from m = 5, a row is rejected
#' if zm < m z0; m is lowered by 0.1 while the rejected rows exceed half of the
#' retained ones, for at most 30 steps (minimum 2). Computed once per site, over
#' the full period, on all rows with zm and z0 (independent of h; spec 08).
#'
#' @param zm,z0 Measurement height above d and roughness length [m].
#' @return Multiplier m [-].
etc_z0_multiplier <- function(zm, z0) {
  ok <- !is.na(zm) & !is.na(z0)
  zm <- zm[ok]
  z0 <- z0[ok]
  m <- 5
  for (i in 1:30) {
    rej <- sum(zm < m * z0)
    ret <- length(zm) - rej
    if (rej == 0) break
    if (ret > 0 && rej <= 0.5 * ret) break
    m <- round(m - 0.1, 10)
  }
  m
}

#' Footprint validity flags per condition and set
#'
#' One logical column per condition (`TRUE` = condition met) and one per set.
#' Conditions:
#' - `inp`: inputs present (zm, z0, h, umean, L, sigma_v, u*, wind direction), not in
#'   an exclusion period, u* not rejected by the spec-05 plausibility rule;
#' - `code_*`: the FFP code's own checks (`Kljun2015` code): h >= 10 m, u* >= 0.1 m s-1,
#'   zm/L > -15.5 (climatology), u* >= 0, sigma_v >= 0;
#' - `domain`: zm < h (the FFP scaling needs 1 - zm/h > 0; the code does not enforce it:
#'   `docs/local/ffp_code_notes.md`);
#' - Eq. 27 (Kljun et al. 2015): `k_rsl` 20 z0 < zm, `k_el` zm < 0.8 h, `k_stab` zm/L >= -15.5;
#' - ETC (`legacy/FFP_analysis_dev.R` lines 1068-1095): `e_rsl` zm >= m z0, `e_el` zm <= 0.8 h,
#'   `e_tower` EC height <= h, `e_stab` -15.5 < zm/L <= 1, `e_ustar` u* >= 0.2 m s-1.
#'   Deviation from the ETC pipeline (G 2026-10-26): stability from zm/L with
#'   zm = EC height - 2/3 hc (as passed to the FFP; Kljun's definition), not from L2 `ZL`
#'   (whose implied z - d differs from the ICOS heights).
#' Sets (all include `inp`, the code checks and `domain`): `kljun` (Eq. 27), `kljun_norsl`
#' (Eq. 27 without the roughness condition), `etc`, `none`.
#'
#' @param zm,z0,ec_height Heights [m].
#' @param h Boundary layer height [m].
#' @param umean Mean wind speed [m s-1].
#' @param L Obukhov length [m].
#' @param sigmav Standard deviation of the lateral wind [m s-1].
#' @param ustar Friction velocity [m s-1].
#' @param wind_dir Wind direction [deg].
#' @param m_etc ETC multiplier ([etc_z0_multiplier()]).
#' @param ustar_implausible,excluded Logical flags (spec 05; exclusion periods).
#' @return Tibble of logical columns (conditions, then `valid_kljun`, `valid_kljun_norsl`,
#'   `valid_etc`, `valid_none`); `NA` inputs give `FALSE`.
footprint_validity <- function(zm, z0, ec_height, h, umean, L, sigmav, ustar, wind_dir, m_etc,
                               ustar_implausible = FALSE, excluded = FALSE) {
  t <- function(x) x %in% TRUE
  zl <- zm / L
  out <- dplyr::tibble(
    inp = t(!is.na(zm) & !is.na(z0) & !is.na(h) & !is.na(umean) & !is.na(L) & !is.na(sigmav) &
              !is.na(ustar) & !is.na(wind_dir) & !ustar_implausible & !excluded),
    code_h10 = t(h >= 10),
    code_ustar01 = t(ustar >= 0.1),
    code_stab = t(zl > -15.5),
    code_nonneg = t(sigmav >= 0 & ustar >= 0 & wind_dir >= 0 & wind_dir <= 360),
    domain = t(zm < h),
    k_rsl = t(zm > 20 * z0),
    k_el = t(zm < 0.8 * h),
    k_stab = t(zl >= -15.5),
    e_rsl = t(zm >= m_etc * z0),
    e_el = t(zm <= 0.8 * h),
    e_tower = t(ec_height <= h),
    e_stab = t(zl > -15.5 & zl <= 1),
    e_ustar = t(ustar >= 0.2)
  )
  base <- out$inp & out$code_h10 & out$code_ustar01 & out$code_stab & out$code_nonneg & out$domain
  out$valid_kljun <- base & out$k_rsl & out$k_el & out$k_stab
  out$valid_kljun_norsl <- base & out$k_el & out$k_stab
  out$valid_etc <- base & out$e_rsl & out$e_el & out$e_tower & out$e_stab & out$e_ustar
  out$valid_none <- base
  out
}

#' Footprint input table of a station
#'
#' Pure function: joins the ICOS L2 table, the PBLHmod output (regime, spec-05 u* flag)
#' and the PBLH cases.
#'
#' @param d Output of [read_icos_site()] (needs `WS`, `MO_LENGTH`, `V_SIGMA`, `USTAR`,
#'   `WD`, `zm_m` = EC height, `canopy_height_m`, `excl_flag`; optional `NEE`,
#'   `NEE_DATA_FLAG`, `fn_NEE_VUT_REF_QC`, `fn_SW_IN_POT`).
#' @param pblhmod_out PBLHmod output (`timestamp_start_utc`, `regime`, `flag_ustar_implausible`).
#' @param cases Output of [footprint_h_cases()] (same rows as `d`).
#' @return Tibble: `timestamp_start_utc`, `ec_height_m`, `canopy_height_m`, `d_m`, `z0_m`,
#'   `zm_m` [m], `umean` [m s-1], `L` [m], `sigmav`, `ustar` [m s-1], `wind_dir` [deg], `nee`,
#'   `nee_data_flag`, `nee_vut_ref_qc`, `day` [logical, SW_IN_POT > 0], `regime`,
#'   `ustar_implausible`, `excluded`, and `h_<case>` [m].
footprint_input_table <- function(d, pblhmod_out, cases) {
  col <- function(k) if (k %in% names(d)) d[[k]] else rep(NA_real_, nrow(d))
  i <- match(as.numeric(d$timestamp_start_utc), as.numeric(pblhmod_out$timestamp_start_utc))
  hz <- footprint_heights(d$zm_m, d$canopy_height_m)
  out <- dplyr::tibble(
    timestamp_start_utc = d$timestamp_start_utc, ec_height_m = d$zm_m, canopy_height_m = d$canopy_height_m,
    d_m = hz$d_m, z0_m = hz$z0_m, zm_m = hz$zm_m, umean = d$WS, L = d$MO_LENGTH, sigmav = d$V_SIGMA,
    ustar = d$USTAR, wind_dir = d$WD, nee = col("NEE"), nee_data_flag = col("NEE_DATA_FLAG"),
    nee_vut_ref_qc = col("fn_NEE_VUT_REF_QC"), day = col("fn_SW_IN_POT") > 0,
    regime = pblhmod_out$regime[i], ustar_implausible = pblhmod_out$flag_ustar_implausible[i] %in% TRUE,
    excluded = d$excl_flag %in% TRUE
  )
  for (k in names(cases)) out[[paste0("h_", k)]] <- cases[[k]]
  out
}

#' Validity and subset columns for every case and set
#'
#' @param tab Output of [footprint_input_table()].
#' @param m_etc ETC multiplier ([etc_z0_multiplier()]).
#' @param sets Validity sets.
#' @return List: `flags` (named list per case of [footprint_validity()] tables) and `subsets`
#'   (tibble with logical columns `<subset>__<set>__<case>`, subset in `valid`, `valid_qc`,
#'   `valid_qc_used`).
footprint_flags_all <- function(tab, m_etc, sets = c("kljun", "kljun_norsl", "etc", "none")) {
  cases <- sub("^h_", "", grep("^h_", names(tab), value = TRUE))
  flags <- list()
  sub <- list()
  for (k in cases) {
    v <- footprint_validity(tab$zm_m, tab$z0_m, tab$ec_height_m, tab[[paste0("h_", k)]], tab$umean, tab$L,
                            tab$sigmav, tab$ustar, tab$wind_dir, m_etc, tab$ustar_implausible, tab$excluded)
    flags[[k]] <- v
    for (s in sets) {
      ss <- footprint_subsets(v[[paste0("valid_", s)]], tab$nee, tab$nee_data_flag, tab$nee_vut_ref_qc)
      for (u in names(ss)) sub[[paste(u, s, k, sep = "__")]] <- ss[[u]]
    }
  }
  list(flags = flags, subsets = dplyr::as_tibble(sub))
}

#' Coverage table per case and validity set
#'
#' @param tab Output of [footprint_input_table()].
#' @param fl Output of [footprint_flags_all()].
#' @param sets Validity sets.
#' @return Long tibble per case, set and stratum (`all`, `regime=`, `season=`, `daynight=`):
#'   `n_rows`, `n_with_h`, `n_inputs`, `n_valid`, `n_valid_qc`, `n_valid_qc_used`, `n_h_lt_zm`,
#'   `n_h_le_10`, and `fail_<condition>` (rows with h and inputs failing the condition; rows
#'   can fail several).
footprint_coverage <- function(tab, fl, sets = c("kljun", "kljun_norsl", "etc", "none")) {
  season <- as.character(season_of(tab$timestamp_start_utc))
  dn <- ifelse(tab$day %in% TRUE, "day", ifelse(tab$day %in% FALSE, "night", NA))
  strata <- c(list(all = rep(TRUE, nrow(tab))),
              stats::setNames(lapply(sort(unique(stats::na.omit(tab$regime))), function(r) tab$regime %in% r),
                              paste0("regime=", sort(unique(stats::na.omit(tab$regime))))),
              stats::setNames(lapply(c("DJF", "MAM", "JJA", "SON"), function(s) season == s),
                              paste0("season=", c("DJF", "MAM", "JJA", "SON"))),
              list(`daynight=day` = dn %in% "day", `daynight=night` = dn %in% "night"))
  conds <- c("code_h10", "code_ustar01", "code_stab", "code_nonneg", "domain", "k_rsl", "k_el", "k_stab",
             "e_rsl", "e_el", "e_tower", "e_stab", "e_ustar")
  rows <- list()
  for (k in names(fl$flags)) {
    v <- fl$flags[[k]]
    h <- tab[[paste0("h_", k)]]
    for (s in sets) {
      for (st in names(strata)) {
        w <- strata[[st]]
        base <- w & !is.na(h) & v$inp
        r <- dplyr::tibble(case = k, set = s, stratum = st, n_rows = sum(w), n_with_h = sum(w & !is.na(h)),
                           n_inputs = sum(base),
                           n_valid = sum(w & fl$subsets[[paste("valid", s, k, sep = "__")]]),
                           n_valid_qc = sum(w & fl$subsets[[paste("valid_qc", s, k, sep = "__")]]),
                           n_valid_qc_used = sum(w & fl$subsets[[paste("valid_qc_used", s, k, sep = "__")]]),
                           n_h_lt_zm = sum(w & h < tab$zm_m, na.rm = TRUE), n_h_le_10 = sum(w & h <= 10, na.rm = TRUE))
        for (cn in conds) r[[paste0("fail_", cn)]] <- sum(base & !v[[cn]])
        rows[[length(rows) + 1]] <- r
      }
    }
  }
  dplyr::bind_rows(rows)
}

#' Data subsets of a valid footprint
#'
#' `valid`; `valid_qc`: and NEE QC good in L2 FLUXES (`NEE_DATA_FLAG == 0`, as ETC);
#' `valid_qc_used`: and NEE not gap-filled in FLUXNET (`NEE_VUT_REF_QC == 0`), i.e. the
#' data actually used for fluxes and budgets. The last step removes the u*-filtered rows
#' (ONEFlux VUT thresholds) and the other FLUXNET removals (G 2026-10-26).
#'
#' @param valid Valid footprint [logical].
#' @param nee FLUXES NEE [umol m-2 s-1].
#' @param nee_data_flag FLUXES `NEE_DATA_FLAG` [-].
#' @param nee_vut_ref_qc FLUXNET `NEE_VUT_REF_QC` [-].
#' @return Tibble: `valid`, `valid_qc`, `valid_qc_used` [logical, nested].
footprint_subsets <- function(valid, nee, nee_data_flag, nee_vut_ref_qc) {
  v <- valid %in% TRUE
  q <- v & !is.na(nee) & nee_data_flag %in% 0
  dplyr::tibble(valid = v, valid_qc = q, valid_qc_used = q & nee_vut_ref_qc %in% 0)
}
