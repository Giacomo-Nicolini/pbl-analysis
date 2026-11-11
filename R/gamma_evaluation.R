# Gamma evaluation at the measured sites (spec 11 D and amendment 1, G 2026-11-11): tables for the decision on the
# EC-only constant gamma (spec 09 decision point 3b). Evidence only; no decision in the code.

#' Season rows of the gamma-case metrics
#'
#' @param wide Wide metrics of `scripts/11_gamma_evaluation.R` (one row per site x case x reference x stratum set;
#'   season sets named `<set>:<season>`; columns of [pblh_metrics()], `n_days`, `mlr_lo`, `mlr_hi`, `slope_sym`).
#' @param min_days Minimum number of days in the bootstrap; fewer days are flagged `low_n` [-].
#' @return Tibble, one row per site x case x reference x set x season: `site_id`, `case`, `reference`,
#'   `stratum_set` (without the season), `season` (DJF, MAM, JJA, SON), `n`, `n_days` [-], `median_log_ratio`,
#'   `mlr_lo`, `mlr_hi` [-, ln(source / reference), day-block bootstrap 95 %], `frac_within_factor2` [0-1], `fb`,
#'   `nrmsd`, `r_spearman`, `slope_sym` [-], `low_n` [logical].
gamma_season_table <- function(wide, min_days = 10L) {
  w <- wide[grepl(":", wide$stratum_set, fixed = TRUE), ]
  parts <- strsplit(w$stratum_set, ":", fixed = TRUE)
  out <- dplyr::tibble(site_id = w$site_id, case = w$case, reference = w$reference,
                       stratum_set = vapply(parts, `[`, "", 1), season = vapply(parts, `[`, "", 2),
                       n = w$n, n_days = w$n_days, median_log_ratio = w$median_log_ratio, mlr_lo = w$mlr_lo,
                       mlr_hi = w$mlr_hi, frac_within_factor2 = w$frac_within_factor2, fb = w$fb, nrmsd = w$nrmsd,
                       r_spearman = w$r_spearman, slope_sym = w$slope_sym)
  out$low_n <- is.na(out$n_days) | out$n_days < min_days
  out
}

#' Decision table of the gamma cases: one row per site
#'
#' Spec 11 amendment 1: per site, the median ln(PBLHmod / measured) of each case on the site's primary set (aerosol:
#' window 10-16 h apparent solar time, PBLHmod convective; FI-Hyy: convective, exact values) with its day-block
#' bootstrap 95 % interval, n, the season with the largest |median ln ratio| (seasons with at least `min_days` days),
#' and the median ln(PBLHmod / ERA5 tower cell) of the case (convective). Sites without a measured reference (BE-Bra)
#' get the ERA5 columns only.
#'
#' @param wide Wide metrics of `scripts/11_gamma_evaluation.R` (see [gamma_season_table()]); may include the rows of
#'   ERA5 against the measurement (case `era5_blh_center`).
#' @param site_ids Stations, in output order [character].
#' @param refs Named list (by station) of [measured_reference()] results (`NULL` = no measured reference).
#' @param cases Cases [character], in column order; `era5_blh_center` = ERA5 tower cell against the measurement.
#' @param min_days Minimum number of days in the bootstrap; fewer days are flagged `low_n_<case>` [-]. No minimum
#'   existed before (spec 11 amendment 1); default 10 days, flagged for G.
#' @return Tibble, one row per station: `site_id`, `reference`, `method_class`, `stratum_set`, `distance_km` [km],
#'   `distance_class`, `borderline`, and per case `<case>`: `n_<case>`, `n_days_<case>` [-], `mlr_<case>`,
#'   `mlr_lo_<case>`, `mlr_hi_<case>` [-, natural log], `low_n_<case>` [logical], `season_max_<case>`,
#'   `season_max_mlr_<case>` [-], `season_max_n_days_<case>` [-], `era5_mlr_<case>`, `era5_n_<case>` (PBLHmod case vs
#'   ERA5 tower cell, convective; not for `era5_blh_center`).
gamma_decision_table <- function(wide, site_ids, refs, cases, min_days = 10L) {
  seasons <- c("DJF", "MAM", "JJA", "SON")
  one <- function(s) {
    rf <- refs[[s]]
    row <- dplyr::tibble(site_id = s, reference = if (is.null(rf)) NA_character_ else rf$src,
                         method_class = if (is.null(rf)) NA_character_ else rf$type,
                         stratum_set = if (is.null(rf)) NA_character_ else rf$primary_set,
                         distance_km = if (is.null(rf)) NA_real_ else rf$distance_km,
                         distance_class = if (is.null(rf)) NA_character_ else rf$distance_class,
                         borderline = if (is.null(rf)) NA else rf$borderline)
    for (k in cases) {
      v <- list(n = NA_integer_, n_days = NA_integer_, mlr = NA_real_, mlr_lo = NA_real_, mlr_hi = NA_real_,
                low_n = NA, season_max = NA_character_, season_max_mlr = NA_real_, season_max_n_days = NA_integer_)
      if (!is.null(rf)) {
        w <- wide[wide$site_id == s & wide$case == k & wide$reference == rf$src, ]
        p <- w[w$stratum_set == rf$primary_set, ]
        if (nrow(p) > 1) cli::cli_abort("{s} {k}: {nrow(p)} rows for the primary set {rf$primary_set}.")
        if (nrow(p) == 1) {
          v[c("n", "n_days", "mlr", "mlr_lo", "mlr_hi")] <- list(as.integer(p$n), as.integer(p$n_days),
                                                                 p$median_log_ratio, p$mlr_lo, p$mlr_hi)
          v$low_n <- is.na(p$n_days) | p$n_days < min_days
        }
        ws <- w[w$stratum_set %in% paste0(rf$primary_set, ":", seasons) & !is.na(w$median_log_ratio) &
                  !is.na(w$n_days) & w$n_days >= min_days, ]
        if (nrow(ws) > 0) {
          j <- which.max(abs(ws$median_log_ratio))
          v[c("season_max", "season_max_mlr", "season_max_n_days")] <-
            list(sub("^.*:", "", ws$stratum_set[j]), ws$median_log_ratio[j], as.integer(ws$n_days[j]))
        }
      }
      names(v) <- paste0(names(v), "_", k)
      row <- dplyr::bind_cols(row, dplyr::as_tibble(v))
      if (k != "era5_blh_center") {
        e <- wide[wide$site_id == s & wide$case == k & wide$reference == "era5_blh_center" &
                    wide$stratum_set == "convective", ]
        row[[paste0("era5_mlr_", k)]] <- if (nrow(e) == 1) e$median_log_ratio else NA_real_
        row[[paste0("era5_n_", k)]] <- if (nrow(e) == 1) as.integer(e$n) else NA_integer_
      }
    }
    row
  }
  dplyr::bind_rows(lapply(site_ids, one))
}
