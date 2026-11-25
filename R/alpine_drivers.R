# Spec 19: why the PBLH source changes the footprint at CH-Dav and IT-Ren (spec-15 exceptions) and not at other
# mountain sites. Site descriptors on the rows used for the footprints (terrain, tower, PBLH, turbulence, night-time
# flows, footprint), row-level R2 of the target-area difference at night, the night-time evidence of the existing
# measured comparisons, and the supplementary figures. Analysis only: no PBLHmod or FFP run.
#
# Conventions:
# - Rows: the spec-15 row set `common_pe` ([fp_common()]) of the reference year, main set `etc`, subsets
#   `valid_qc_used` (main) and `valid`. Where the chunks of a case are absent (local test without FFP results) the
#   row basis is `flags_only`: valid in the three cases by the step-08a flags, no FFP result required; the footprint
#   columns are then NA. The basis is a column of every output.
# - Strata: zm/L class of spec 08 (`zl_class()`), `neutral_stable` = neutral or stable, `all` = margin; day / night
#   from the input table (`day` = SW_IN_POT > 0).
# - Wind direction: meteorological (direction the wind blows from) [deg from N, 0-360).

# --- Circular statistics -------------------------------------------------------------------------------------------

#' Signed difference of two directions
#'
#' @param a,b Directions [deg].
#' @return a - b wrapped to (-180, 180] [deg].
circ_diff_deg <- function(a, b) {
  d <- (a - b) %% 360
  ifelse(d > 180, d - 360, d)
}

#' Circular mean direction (unit vectors, unweighted)
#'
#' @param wd Directions [deg from N]; non-finite values dropped.
#' @return Mean direction [deg, 0-360); `NA` without data or when the mean unit vector has zero length (< 1e-9).
circ_mean_deg <- function(wd) {
  wd <- wd[is.finite(wd)]
  if (length(wd) == 0) return(NA_real_)
  s <- mean(sin(wd * pi / 180))
  c0 <- mean(cos(wd * pi / 180))
  if (sqrt(s^2 + c0^2) < 1e-9) return(NA_real_)
  round(atan2(s, c0) * 180 / pi, 9) %% 360  # rounding: -1e-15 must give 0, not 360
}

#' Directional steadiness of the wind
#'
#' |vector mean wind| / scalar mean speed (1 = constant direction, 0 = no preferred direction).
#'
#' @param wd Wind direction [deg from N]; @param ws Wind speed [m s-1]; rows with both finite are used.
#' @return Steadiness [0-1]; `NA` without data or with zero mean speed.
wind_steadiness <- function(wd, ws) {
  ok <- is.finite(wd) & is.finite(ws)
  if (!any(ok) || mean(ws[ok]) <= 0) return(NA_real_)
  u <- -ws[ok] * sin(wd[ok] * pi / 180)
  v <- -ws[ok] * cos(wd[ok] * pi / 180)
  sqrt(mean(u)^2 + mean(v)^2) / mean(ws[ok])
}

#' Dominant 90-degree direction sector
#'
#' The sector centre (in steps of `step`) holding the largest share of the directions within +-`half`.
#'
#' @param wd Directions [deg from N]; @param half Half-width of the sector [deg]; @param step Centre step [deg].
#' @return List: `centre` [deg; first centre on ties], `share` [0-1]; `NA` without data.
dominant_sector <- function(wd, half = 45, step = 5) {
  wd <- wd[is.finite(wd)]
  if (length(wd) == 0) return(list(centre = NA_real_, share = NA_real_))
  centres <- seq(0, 360 - step, by = step)
  sh <- vapply(centres, function(c0) mean(abs(circ_diff_deg(wd, c0)) <= half), numeric(1))
  k <- which.max(sh)
  list(centre = centres[k], share = sh[k])
}

# --- Rows ----------------------------------------------------------------------------------------------------------

#' Columns of the step-08a input table copied to the spec-19 rows
#' @return Character vector.
alpine_tab_cols <- function() {
  c("ec_height_m", "canopy_height_m", "zm_m", "umean", "L", "ustar", "wind_dir", "h_pblhmod", "h_pblhmod_clim",
    "h_era5")
}

#' Row basis of the spec-19 rows
#'
#' @param r Rows of [fp_site_rows()] with the validity columns.
#' @param set,subset Validity set and subset [character].
#' @param require_ffp `TRUE`: `common_pe` (valid and an FFP result in pblhmod, pblhmod_clim, era5; [fp_common()]);
#'   `FALSE`: valid in the three cases by the flags only (`flags_only`, no chunks).
#' @return Logical, one per row.
alpine_basis <- function(r, set, subset, require_ffp = TRUE) {
  if (require_ffp) return(fp_common(r, set, subset))
  ok <- rep(TRUE, nrow(r))
  for (cs in c("pblhmod", "pblhmod_clim", "era5")) {
    v <- r[[paste("v", set, subset, cs, sep = "__")]]
    if (is.null(v)) return(rep(FALSE, nrow(r)))
    ok <- ok & v
  }
  ok
}

#' Derived row variables of spec 19
#'
#' @param r Rows of [fp_site_rows()] with [alpine_tab_cols()] (one site).
#' @param w_sigma Standard deviation of the vertical wind (L2 FLUXES `W_SIGMA`) on the rows of `r` [m s-1].
#' @param ustar_thr ONEFlux VUT u* threshold of each row's year [m s-1] (`NA` where not available).
#' @return `r` with `ln_h_pblhmod_era5` [-, ln(h PBLHmod / h ERA5)], `zm_h_era5`, `zm_h_pblhmod` [-, zm / h],
#'   `w_sigma` [m s-1], `sw_ustar` [-, sigma_w / u*], `below_ustar_thr` [logical], `ta_diff_pp` [pp, TA ERA5 - TA
#'   PBLHmod], `lnx80_era5_pblhmod` [-, ln(x_80 ERA5 / x_80 PBLHmod)].
alpine_row_vars <- function(r, w_sigma, ustar_thr) {
  pos <- function(x) ifelse(is.finite(x) & x > 0, x, NA_real_)
  r$ln_h_pblhmod_era5 <- log(pos(r$h_pblhmod) / pos(r$h_era5))
  r$zm_h_era5 <- r$zm_m / pos(r$h_era5)
  r$zm_h_pblhmod <- r$zm_m / pos(r$h_pblhmod)
  r$w_sigma <- w_sigma
  r$sw_ustar <- w_sigma / pos(r$ustar)
  r$below_ustar_thr <- ifelse(is.finite(r$ustar) & is.finite(ustar_thr), r$ustar < ustar_thr, NA)
  ta <- function(cs) if (is.null(r[[paste0("ta__", cs)]])) rep(NA_real_, nrow(r)) else r[[paste0("ta__", cs)]]
  x80 <- function(cs) if (is.null(r[[paste0("x80__", cs)]])) rep(NA_real_, nrow(r)) else r[[paste0("x80__", cs)]]
  r$ta_diff_pp <- 100 * (ta("era5") - ta("pblhmod"))
  r$lnx80_era5_pblhmod <- log(pos(x80("era5")) / pos(x80("pblhmod")))
  r
}

#' Strata of spec 19 (zm/L class x day/night, with margins and neutral + stable)
#'
#' @param zl zm/L class per row (`unstable`, `neutral`, `stable`, `NA`); @param dn `day` / `night` per row.
#' @return List of lists: `regime`, `daynight`, `i` (logical row index).
alpine_strata <- function(zl, dn) {
  regs <- c("all", "unstable", "neutral", "stable", "neutral_stable")
  out <- list()
  for (g in regs) for (t in c("all", "day", "night")) {
    ig <- switch(g, all = rep(TRUE, length(zl)), neutral_stable = zl %in% c("neutral", "stable"), zl %in% g)
    it <- if (t == "all") rep(TRUE, length(dn)) else dn %in% t
    out[[length(out) + 1]] <- list(regime = g, daynight = t, i = ig & it)
  }
  out
}

q_fin <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x)) unname(stats::quantile(x, p, names = FALSE)) else NA_real_
}
med_fin <- function(x) q_fin(x, 0.5)

#' Stratified statistics of one site's rows
#'
#' @param d Rows of one site and stratum ([alpine_row_vars()]).
#' @param min_n Smallest n of a flagged spec-15 exception cell [-]; @param threshold_pp Exception threshold [pp].
#' @return Tibble `block`, `variable`, `unit`, `value` (numeric), `value_chr` (character variables).
alpine_stratum_stats <- function(d, threshold_pp = 5, min_n = 50) {
  v <- function(block, variable, unit, value) dplyr::tibble(block = block, variable = variable, unit = unit,
                                                            value = as.numeric(value), value_chr = NA_character_)
  ta <- d$ta_diff_pp[is.finite(d$ta_diff_pp)]
  p10 <- q_fin(ta, 0.1)
  p90 <- q_fin(ta, 0.9)
  status <- if (length(ta) == 0) NA_character_ else if (abs(p10) > threshold_pp || abs(p90) > threshold_pp) {
    if (length(ta) >= min_n) "exception" else "below_min_n"
  } else "none"
  dplyr::bind_rows(
    v("rows", "n", "-", nrow(d)),
    v("pblh", "h_era5_median", "m", med_fin(d$h_era5)),
    v("pblh", "h_pblhmod_median", "m", med_fin(d$h_pblhmod)),
    v("pblh", "h_pblhmod_clim_median", "m", med_fin(d$h_pblhmod_clim)),
    v("pblh", "ln_h_pblhmod_era5_median", "-", med_fin(d$ln_h_pblhmod_era5)),
    v("pblh", "zm_h_era5_median", "-", med_fin(d$zm_h_era5)),
    v("pblh", "zm_h_era5_p90", "-", q_fin(d$zm_h_era5, 0.9)),
    v("pblh", "zm_h_pblhmod_median", "-", med_fin(d$zm_h_pblhmod)),
    v("pblh", "zm_h_pblhmod_p90", "-", q_fin(d$zm_h_pblhmod, 0.9)),
    v("turbulence", "ustar_median", "m s-1", med_fin(d$ustar)),
    v("turbulence", "ustar_p10", "m s-1", q_fin(d$ustar, 0.1)),
    v("turbulence", "ustar_p90", "m s-1", q_fin(d$ustar, 0.9)),
    v("turbulence", "sw_ustar_median", "-", med_fin(d$sw_ustar)),
    v("turbulence", "sw_ustar_p10", "-", q_fin(d$sw_ustar, 0.1)),
    v("turbulence", "sw_ustar_p90", "-", q_fin(d$sw_ustar, 0.9)),
    v("turbulence", "L_median", "m", med_fin(d$L)),
    v("turbulence", "L_p10", "m", q_fin(d$L, 0.1)),
    v("turbulence", "L_p90", "m", q_fin(d$L, 0.9)),
    v("turbulence", "share_below_ustar_thr", "0-1",
      if (any(!is.na(d$below_ustar_thr))) mean(d$below_ustar_thr, na.rm = TRUE) else NA_real_),
    v("flows", "wd_circ_mean", "deg", circ_mean_deg(d$wind_dir)),
    v("flows", "wd_steadiness", "0-1", wind_steadiness(d$wind_dir, d$umean)),
    v("flows", "wd_dominant_centre", "deg", dominant_sector(d$wind_dir)$centre),
    v("flows", "wd_share_dominant_90", "0-1", dominant_sector(d$wind_dir)$share),
    v("flows", "ws_median", "m s-1", med_fin(d$umean)),
    v("footprint", "ta_diff_n", "-", length(ta)),
    v("footprint", "ta_diff_p10", "pp", p10),
    v("footprint", "ta_diff_median", "pp", q_fin(ta, 0.5)),
    v("footprint", "ta_diff_p90", "pp", p90),
    dplyr::tibble(block = "footprint", variable = "spec15_status", unit = "-", value = NA_real_, value_chr = status)
  )
}

#' Site descriptors of one site (spec 19, section 1)
#'
#' @param r Rows of one site ([alpine_row_vars()]) with the validity columns.
#' @param bases Named logical list: row basis per subset (names = subsets, e.g. `valid_qc_used`, `valid`).
#' @param row_basis `common_pe` or `flags_only` [character].
#' @param site Site-level values: list with `cells` (one row of `09_network_era5_cells.csv`, may be `NULL`), `igbp`,
#'   `ustar_thr_source` (character, e.g. `VUT` or the reason it is missing).
#' @inheritParams alpine_stratum_stats
#' @return Long tibble: `site_id`, `block`, `variable`, `unit`, `subset`, `regime`, `daynight`, `value`, `value_chr`,
#'   `row_basis`. Site-level variables have `subset`, `regime`, `daynight` `NA` (tower: on the main subset rows).
alpine_site_descriptors <- function(r, bases, row_basis, site, threshold_pp = 5, min_n = 50) {
  sid <- r$site_id[1]
  out <- list()
  for (u in names(bases)) {
    b <- bases[[u]]
    for (s in alpine_strata(r$zl_class[b], r$daynight[b])) {
      st <- alpine_stratum_stats(r[b, ][s$i, ], threshold_pp, min_n)
      out[[length(out) + 1]] <- dplyr::mutate(st, subset = u, regime = s$regime, daynight = s$daynight)
    }
    # reversal: night minus day circular mean direction, all regimes
    rn <- r[b & r$daynight %in% "night", ]
    rd <- r[b & r$daynight %in% "day", ]
    out[[length(out) + 1]] <- dplyr::tibble(block = "flows", variable = "wd_reversal_night_minus_day", unit = "deg",
                                            value = circ_diff_deg(circ_mean_deg(rn$wind_dir), circ_mean_deg(rd$wind_dir)),
                                            value_chr = NA_character_, subset = u, regime = "all", daynight = "all")
  }
  b1 <- bases[[1]]
  cells <- site$cells
  cv <- function(k) if (is.null(cells) || !k %in% names(cells) || nrow(cells) == 0) NA_real_ else as.numeric(cells[[k]][1])
  sl <- function(block, variable, unit, value, chr = NA_character_) {
    dplyr::tibble(block = block, variable = variable, unit = unit, value = as.numeric(value), value_chr = chr,
                  subset = NA_character_, regime = NA_character_, daynight = NA_character_)
  }
  out[[length(out) + 1]] <- dplyr::bind_rows(
    sl("terrain", "site_altitude_m", "m a.s.l.", cv("site_altitude_m")),
    sl("terrain", "era5_orog_m", "m a.s.l.", cv("orog_m")),
    sl("terrain", "site_minus_orog_m", "m", cv("dalt_m")),
    sl("terrain", "sdor_max_m", "m", cv("sdor_max_m")),
    sl("terrain", "orog_range_9cells_m", "m", cv("orog_max_m") - cv("orog_min_m")),
    sl("tower", "ec_height_m", "m", med_fin(r$ec_height_m[b1])),
    sl("tower", "canopy_height_m", "m", med_fin(r$canopy_height_m[b1])),
    sl("tower", "zm_m", "m", med_fin(r$zm_m[b1])),
    sl("tower", "igbp", "-", NA, if (is.null(site$igbp)) NA_character_ else site$igbp),
    sl("turbulence", "ustar_thr_source", "-", NA, site$ustar_thr_source),
    sl("rows", "year", "-", as.integer(format(r$timestamp_start_utc[1], "%Y", tz = "UTC")))
  )
  dplyr::mutate(dplyr::bind_rows(out), site_id = sid, row_basis = row_basis, .before = 1)[
    c("site_id", "block", "variable", "unit", "subset", "regime", "daynight", "value", "value_chr", "row_basis")]
}

#' Header lines of the spec-19 tables
#' @param threshold_pp,min_n Spec-15 exception rule.
#' @return Named list of character vectors.
alpine_table_headers <- function(threshold_pp = 5, min_n = 50) {
  common <- c(
    "Spec 19 (scripts/19a_site_descriptors.R, R/alpine_drivers.R): reference year per site (config/reference_years.csv).",
    "Rows: row_basis common_pe = spec-15 rows (valid in set etc x subset, with an FFP result in pblhmod, pblhmod_clim and era5); flags_only = valid by the step-08a flags in the three cases, no FFP result (local test, footprint values NA).",
    "Strata: regime = zm/L class (unstable < -0.1 <= neutral <= 0.1 < stable; zm = EC height - 2/3 hc), neutral_stable = neutral or stable, all = margin; daynight from the input table (day = SW_IN_POT > 0).")
  list(
    descriptors = c(common,
      "Long table: one value per site x variable x subset x stratum; site-level variables (terrain, tower, u* threshold source) have subset/regime/daynight empty; character values in value_chr.",
      "terrain: from 09_network_era5_cells.csv (site altitude, ERA5 tower-cell orography = z / 9.80665, site minus orography, max sub-grid orography SD sdor of the 9 cells, orography max - min of the 9 cells) [m].",
      "tower: medians of the main-subset rows: EC height, canopy height hc, zm = EC height - 2/3 hc [m]; igbp from config/sites.csv.",
      "pblh: medians of h ERA5 (tower cell), PBLHmod, PBLHmod-clim [m]; ln_h_pblhmod_era5 = ln(h PBLHmod / h ERA5) [-]; zm / h median and p90 [-].",
      "turbulence: u* [m s-1], sigma_w / u* [-] (W_SIGMA of L2 FLUXES), L [m] (median, p10, p90); share_below_ustar_thr = share of rows with u* below the ONEFlux VUT reference threshold of the year (BIF file; ustar_thr_source says whether available); 0 by construction on valid_qc_used (FLUXNET NEE_VUT_REF_QC == 0 removes the u*-filtered rows), informative on valid.",
      "flows: wd_circ_mean = circular mean of the wind direction (unit vectors) [deg from N]; wd_steadiness = |vector mean wind| / scalar mean speed [0-1]; wd_dominant_centre = centre of the 90-degree sector (5-degree steps) holding most rows, wd_share_dominant_90 = its share [0-1]; ws_median [m s-1]; wd_reversal_night_minus_day = circular mean night minus day, wrapped to (-180, 180] [deg].",
      sprintf("footprint: TA difference ERA5 - PBLHmod [pp] (TA = footprint fraction in the ICOS target area): n, p10, median, p90; spec15_status = exception if |p10| or |p90| > %g pp and n >= %g, below_min_n if above the threshold with fewer rows, none otherwise (regime x day/night cells as in spec 15).", threshold_pp, min_n)),
    row_level = c(common,
      "Night rows of the selected sites, regime neutral, stable and neutral_stable; simple linear models (lm) of the response on the predictors; r2 = R squared, coef_1 / coef_2 = slopes in the order of the predictors.",
      "Responses: ta_diff_pp = TA ERA5 - TA PBLHmod [pp]; lnx80_era5_pblhmod = ln(x_80 ERA5 / x_80 PBLHmod) [-]. Predictors: ln_h_era5_pblhmod = ln(h ERA5 / h PBLHmod) [-]; zm_h_era5 = zm / h ERA5 [-]."),
    rows = c(common, "Night rows (neutral, stable) of the selected sites for figure 19s_row_level; columns as in 19_row_level.csv."),
    flows = c(common, "Wind-direction histogram per site x daynight (main subset, all regimes): 30-degree bins centred on 0, 30, ..., 330 deg from N; n and share [0-1] of the site x daynight rows."),
    closeness = c("Spec 19 section 3: night-time and daytime evidence from the existing measured comparisons (no rerun).",
      "FI-Hyy: 07_metrics_FI-Hyy.csv (Cloudnet epsilon MLH, valid at night), strata all and PBLHmod regime; CH-Dav and the other aerosol sites: 10_metrics_<site>.csv, stratum set window_convective (daytime window, PBLHmod convective), stratum all.",
      "median_log_ratio = median ln(source / measured) [-]; bias = mean(source - measured) [m]; mad_log_ratio = median |ln ratio| [-] where available. No measured night-time PBLH exists at the alpine sites (CH-Dav: STRATfinder night values on a floor, spec 11 C; IT-Ren: no measurement)."))
}

# --- Section 2: row level ------------------------------------------------------------------------------------------

#' Night rows of spec 19 section 2
#'
#' @param r Rows of one or more sites ([alpine_row_vars()]); @param basis Logical row basis (main subset).
#' @return Tibble: `site_id`, `zl_class` (neutral / stable), `ta_diff_pp` [pp], `lnx80_era5_pblhmod` [-],
#'   `ln_h_era5_pblhmod` [-, = -ln_h_pblhmod_era5], `zm_h_era5` [-].
alpine_night_rows <- function(r, basis) {
  i <- basis & r$daynight %in% "night" & r$zl_class %in% c("neutral", "stable")
  dplyr::tibble(site_id = r$site_id[i], zl_class = r$zl_class[i], ta_diff_pp = r$ta_diff_pp[i],
                lnx80_era5_pblhmod = r$lnx80_era5_pblhmod[i], ln_h_era5_pblhmod = -r$ln_h_pblhmod_era5[i],
                zm_h_era5 = r$zm_h_era5[i])
}

#' Row-level R2 of the night rows (spec 19 section 2)
#'
#' @param nr Output of [alpine_night_rows()].
#' @return Tibble: `site_id`, `stratum` (`neutral_night`, `stable_night`, `neutral_stable_night`), `response`,
#'   `predictors`, `n`, `r2`, `coef_1`, `coef_2`; `NA` with fewer than 3 complete rows.
alpine_row_level <- function(nr) {
  models <- list(ln_h_era5_pblhmod = "ln_h_era5_pblhmod", zm_h_era5 = "zm_h_era5",
                 both = c("ln_h_era5_pblhmod", "zm_h_era5"))
  out <- list()
  for (s in unique(nr$site_id)) for (st in c("neutral_night", "stable_night", "neutral_stable_night")) {
    cls <- switch(st, neutral_night = "neutral", stable_night = "stable", c("neutral", "stable"))
    d <- nr[nr$site_id == s & nr$zl_class %in% cls, ]
    for (y in c("ta_diff_pp", "lnx80_era5_pblhmod")) for (m in names(models)) {
      xs <- models[[m]]
      ok <- is.finite(d[[y]]) & Reduce(`&`, lapply(xs, function(k) is.finite(d[[k]])))
      n <- sum(ok)
      r2 <- NA_real_
      cf <- c(NA_real_, NA_real_)
      if (n >= 3 && stats::var(d[[y]][ok]) > 0) {
        fit <- stats::lm(stats::reformulate(xs, y), data = d[ok, ])
        yy <- d[[y]][ok]
        r2 <- 1 - sum(stats::residuals(fit)^2) / sum((yy - mean(yy))^2)
        b <- unname(stats::coef(fit)[-1])
        cf[seq_along(b)] <- b
      }
      out[[length(out) + 1]] <- dplyr::tibble(site_id = s, stratum = st, response = y, predictors = m, n = n,
                                              r2 = r2, coef_1 = cf[1], coef_2 = cf[2])
    }
  }
  dplyr::bind_rows(out)
}

#' Wind-direction histogram per site and day / night
#'
#' @param r Rows of one site; @param basis Logical row basis; @param width Bin width [deg].
#' @return Tibble: `site_id`, `daynight`, `dir_bin` [deg, bin centre], `n`, `share` [0-1].
alpine_wd_hist <- function(r, basis, width = 30) {
  centres <- seq(0, 360 - width, by = width)
  out <- list()
  for (t in c("day", "night")) {
    wd <- r$wind_dir[basis & r$daynight %in% t]
    wd <- wd[is.finite(wd)]
    k <- if (length(wd)) ((round(wd / width) %% length(centres)) + 1) else integer()
    n <- tabulate(k, length(centres))
    out[[t]] <- dplyr::tibble(site_id = r$site_id[1], daynight = t, dir_bin = centres, n = n,
                              share = if (sum(n) > 0) n / sum(n) else NA_real_)
  }
  dplyr::bind_rows(out)
}

# --- Section 3: existing measured comparisons ----------------------------------------------------------------------

#' Night-time and daytime evidence from the existing measured comparisons (spec 19 section 3)
#'
#' @param tables_dir Folder with `07_metrics_FI-Hyy.csv` and `10_metrics_<site>.csv` [character].
#' @return Tibble: `site_id`, `table`, `source` (`pblhmod_default`, `era5_center`), `stratum_set`, `stratum`, `n`,
#'   `median_log_ratio`, `mad_log_ratio`, `bias_m`; empty if no table is found.
alpine_night_closeness <- function(tables_dir) {
  rd <- function(p) {
    x <- readr::read_csv(p, comment = "#", show_col_types = FALSE, progress = FALSE)
    # tables written before the pairing columns (spec 13) have no site_id: take it from the file name
    if (!"site_id" %in% names(x)) x$site_id <- rep(sub("^.*_([A-Z]{2}-[A-Za-z0-9]{3})[.]csv$", "\\1", basename(p)), nrow(x))
    x
  }
  wide <- function(x, tab) {
    if (nrow(x) == 0) return(NULL)
    x$source <- ifelse(grepl("^1_pblhmod_default", x$pair), "pblhmod_default", "era5_center")
    k <- unique(x[c("site_id", "source", "stratum_set", "stratum")])
    g <- function(i, m) {
      v <- x$value[x$site_id == k$site_id[i] & x$source == k$source[i] & x$stratum_set == k$stratum_set[i] &
                     x$stratum == k$stratum[i] & x$metric == m]
      if (length(v)) v[1] else NA_real_
    }
    dplyr::tibble(site_id = k$site_id, table = tab, source = k$source, stratum_set = k$stratum_set, stratum = k$stratum,
                  n = vapply(seq_len(nrow(k)), g, numeric(1), m = "n"),
                  median_log_ratio = vapply(seq_len(nrow(k)), g, numeric(1), m = "median_log_ratio"),
                  mad_log_ratio = vapply(seq_len(nrow(k)), g, numeric(1), m = "mad_log_ratio"),
                  bias_m = vapply(seq_len(nrow(k)), g, numeric(1), m = "bias"))
  }
  out <- list()
  p7 <- file.path(tables_dir, "07_metrics_FI-Hyy.csv")
  if (file.exists(p7)) {
    x <- rd(p7)
    x <- x[x$pair %in% c("1_pblhmod_default_vs_measured", "2_era5_center_vs_measured") &
             (x$stratum == "all" | grepl("^regime=", x$stratum)), ]
    x$stratum_set <- "all_measured"
    out[[length(out) + 1]] <- wide(x, basename(p7))
  }
  for (p in list.files(tables_dir, "^10_metrics_.*[.]csv$", full.names = TRUE)) {
    x <- rd(p)
    if (!"stratum_set" %in% names(x)) next
    x <- x[grepl("^(1_pblhmod_default|2_era5_center)_vs_", x$pair) & x$stratum_set == "window_convective" &
             x$stratum == "all", ]
    out[[length(out) + 1]] <- wide(x, basename(p))
  }
  x <- dplyr::bind_rows(out)
  if (nrow(x) == 0) {
    return(dplyr::tibble(site_id = character(), table = character(), source = character(), stratum_set = character(),
                         stratum = character(), n = numeric(), median_log_ratio = numeric(), mad_log_ratio = numeric(),
                         bias_m = numeric()))
  }
  x
}

# --- Figures -------------------------------------------------------------------------------------------------------

#' Table of figure 19s_alpine_drivers
#'
#' @param desc Descriptor table ([alpine_site_descriptors()], several sites).
#' @param subset Subset [character]; @param label_sites Sites always labelled [character].
#' @param label_pp Label sites with |p10| above this [pp]; @param min_n Open symbol below this n [-].
#' @return Tibble per site x panel (`stable` night, `neutral` night): `site_id`, `panel`, `x` [-, night
#'   neutral_stable median ln(PBLHmod / ERA5)], `y` [pp, p10 of the TA difference], `n`, `ec_height_m`, `sdor_max_m`,
#'   `below_min_n`, `label` (site id or empty).
alpine_drivers_table <- function(desc, subset = "valid_qc_used", label_sites = c("CZ-BK1", "IT-MBo", "FR-Lqu"),
                                 label_pp = 2, min_n = 50) {
  g <- function(var, reg, dn = "night", sub = subset) {
    x <- desc[desc$variable == var & desc$subset %in% sub & desc$regime %in% reg & desc$daynight %in% dn, ]
    stats::setNames(x$value, x$site_id)
  }
  s1 <- function(var) {
    x <- desc[desc$variable == var & is.na(desc$subset), ]
    stats::setNames(x$value, x$site_id)
  }
  sites <- unique(desc$site_id)
  xv <- g("ln_h_pblhmod_era5_median", "neutral_stable")
  eh <- s1("ec_height_m")
  sd <- s1("sdor_max_m")
  out <- dplyr::bind_rows(lapply(c("stable", "neutral"), function(p) {
    y <- g("ta_diff_p10", p)
    n <- g("ta_diff_n", p)
    dplyr::tibble(site_id = sites, panel = paste(p, "night"), x = unname(xv[sites]), y = unname(y[sites]),
                  n = unname(n[sites]), ec_height_m = unname(eh[sites]), sdor_max_m = unname(sd[sites]))
  }))
  out$below_min_n <- !is.na(out$n) & out$n < min_n
  out$label <- ifelse(out$site_id %in% label_sites | (is.finite(out$y) & abs(out$y) > label_pp), out$site_id, "")
  out
}

#' Figure 19s_alpine_drivers
#'
#' @param ft Output of [alpine_drivers_table()]; @param threshold_pp Exception threshold [pp] (line at minus it).
#' @param caption Caption [character].
#' @return ggplot.
plot_alpine_drivers <- function(ft, threshold_pp = 5, caption = NULL) {
  ft <- ft[is.finite(ft$x) & is.finite(ft$y), ]
  ft$symbol <- ifelse(ft$below_min_n, "n < min", "n >= min")
  ggplot2::ggplot(ft, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = -threshold_pp, linetype = style_of("reference_line", "linetype", "element"),
                        colour = style_of("reference_line", "colour", "element")) +
    ggplot2::geom_point(ggplot2::aes(colour = ec_height_m, size = sdor_max_m, shape = symbol)) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.6, vjust = -1.1) +
    ggplot2::scale_shape_manual(values = c(`n >= min` = 16, `n < min` = 1), name = "rows in the stratum") +
    ggplot2::scale_colour_gradient(low = style_of("grad_low", "colour", "element"),
                                   high = style_of("grad_high", "colour", "element"), name = "EC height [m]") +
    ggplot2::scale_size_continuous(range = c(1.2, 5), name = "sub-grid orography SD [m]") +
    ggplot2::facet_wrap(~panel) +
    ggplot2::labs(x = "night (neutral + stable) median ln(PBLHmod / ERA5) [-]",
                  y = "p10 of TA(ERA5) - TA(PBLHmod) [pp]", caption = caption) +
    theme_pblh() + ggplot2::theme(legend.box = "vertical")
}

#' Figure 19s_row_level
#'
#' @param nr Night rows ([alpine_night_rows()]); @param caption Caption.
#' @return ggplot.
plot_alpine_row_level <- function(nr, caption = NULL) {
  nr <- nr[is.finite(nr$ta_diff_pp) & is.finite(nr$ln_h_era5_pblhmod), ]
  ggplot2::ggplot(nr, ggplot2::aes(x = ln_h_era5_pblhmod, y = ta_diff_pp, colour = zm_h_era5)) +
    ggplot2::geom_hline(yintercept = 0, colour = style_of("reference_line", "colour", "element")) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6) +
    ggplot2::scale_colour_gradient(low = style_of("grad_low", "colour", "element"),
                                   high = style_of("grad_high", "colour", "element"), name = "zm / h ERA5 [-]",
                                   trans = "log10") +
    ggplot2::facet_wrap(~site_id) +
    ggplot2::labs(x = "ln(h ERA5 / h PBLHmod) [-]", y = "TA(ERA5) - TA(PBLHmod) [pp]", caption = caption) +
    theme_pblh()
}

#' Figure 19s_night_flows: day and night wind-direction distributions (polar histograms)
#'
#' @param wh Output of [alpine_wd_hist()] (several sites); @param caption Caption.
#' @return ggplot.
plot_alpine_flows <- function(wh, caption = NULL) {
  wh$daynight <- factor(wh$daynight, levels = c("day", "night"))
  ggplot2::ggplot(wh, ggplot2::aes(x = dir_bin, y = share, fill = daynight)) +
    ggplot2::geom_col(width = 30, position = "identity", alpha = 0.55,
                      colour = style_of("reference_line", "colour", "element"), linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = style_of(c("day", "night"), "colour", "fp_share"),
                               labels = style_of(c("day", "night"), "label", "fp_share"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = c(0, 90, 180, 270), labels = c("N", "E", "S", "W"), limits = c(-15, 345)) +
    ggplot2::coord_polar(start = -15 * pi / 180) +
    ggplot2::facet_wrap(~site_id) +
    ggplot2::labs(x = NULL, y = "share of rows [-]", caption = caption) +
    theme_pblh()
}
