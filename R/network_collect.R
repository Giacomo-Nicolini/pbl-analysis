# Network collection (spec 09 A4): per-site outputs of steps 01, 11a, 02, 04, 08a bound into network tables, the
# reference-year proposal (spec 08 section 4c, restated in spec 09 A5) and the site flags. Pure functions; reading
# and writing in scripts/09_collect.R.

#' Usable half-hours per calendar year (reference-year rule)
#'
#' A half-hour is usable when both main cases (`pblhmod` and `era5`) have h and FFP inputs, are valid under the `etc`
#' set and in the `valid_qc_used` subset (spec 09 A5). The FLUXES reference period is the input table itself; ERA5
#' (from 2020) is implied by the `era5` h.
#'
#' @param tab Input table of `scripts/08a_footprint_flags.R` (`timestamp_start_utc`, `h_measured` if any).
#' @param subsets Subset columns (`valid_qc_used__etc__pblhmod`, `valid_qc_used__etc__era5`).
#' @return Tibble per year: `year`, `n_usable` [half-hours], `n_months` [months with at least one usable half-hour],
#'   `n_measured` [usable half-hours with a measured h].
footprint_year_counts <- function(tab, subsets) {
  get <- function(k) { v <- subsets[[k]]; if (is.null(v)) rep(FALSE, nrow(tab)) else v }
  ok <- get("valid_qc_used__etc__pblhmod") & get("valid_qc_used__etc__era5")
  meas <- if ("h_measured" %in% names(tab)) !is.na(tab$h_measured) else rep(FALSE, nrow(tab))
  yr <- as.integer(format(tab$timestamp_start_utc, "%Y"))
  mo <- format(tab$timestamp_start_utc, "%m")
  dplyr::bind_rows(lapply(sort(unique(yr)), function(y) {
    i <- yr == y
    dplyr::tibble(year = y, n_usable = sum(ok & i), n_months = length(unique(mo[ok & i])), n_measured = sum(ok & i & meas))
  }))
}

#' Reference-year proposal of one site
#'
#' Spec 09 A5 (spec 08 section 4c): among the years with at least `min_months` months of data, the one with the most
#' usable half-hours; years within `tol` of the maximum are equivalent, among them a year overlapping the site's
#' measured source is preferred (most measured half-hours), otherwise the most recent.
#'
#' @param counts Output of [footprint_year_counts()].
#' @param min_months Minimum months with usable data [-].
#' @param tol Relative tolerance on the number of usable half-hours [-].
#' @return One-row tibble: `year` (`NA` if no year qualifies), `n_usable`, `n_months`, `n_measured`, `n_candidates`
#'   (years within `tol`), `reason`.
propose_reference_year <- function(counts, min_months = 9, tol = 0.02) {
  ok <- counts[counts$n_months >= min_months & counts$n_usable > 0, ]
  if (nrow(ok) == 0) {
    return(dplyr::tibble(year = NA_integer_, n_usable = NA_integer_, n_months = NA_integer_, n_measured = NA_integer_,
                         n_candidates = 0L, reason = paste0("no year with >= ", min_months, " months of usable data")))
  }
  cand <- ok[ok$n_usable >= (1 - tol) * max(ok$n_usable), ]
  if (any(cand$n_measured > 0)) {
    pick <- cand[order(-cand$n_measured, -cand$year), ][1, ]
    why <- "most usable half-hours (within 2 %), overlapping the measured source"
  } else {
    pick <- cand[order(-cand$year), ][1, ]
    why <- if (nrow(cand) > 1) "most usable half-hours (within 2 %), most recent" else "most usable half-hours"
  }
  dplyr::tibble(year = as.integer(pick$year), n_usable = as.integer(pick$n_usable), n_months = as.integer(pick$n_months),
                n_measured = as.integer(pick$n_measured), n_candidates = nrow(cand), reason = why)
}

#' Read and bind per-site tables
#'
#' @param paths Named character vector (names = site ids) of csv paths; missing files are skipped.
#' @return Tibble with `site_id` first (from the file if present, else the name), or an empty tibble.
bind_site_tables <- function(paths) {
  out <- lapply(names(paths), function(s) {
    p <- paths[[s]]
    if (!file.exists(p)) return(NULL)
    x <- readr::read_csv(p, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
    if (!"site_id" %in% names(x)) x <- dplyr::mutate(x, site_id = s, .before = 1)
    x$site_id <- as.character(x$site_id)
    x
  })
  dplyr::bind_rows(out)
}

#' Network gamma summary per month
#'
#' @param monthly Site x month climatologies (`site_id`, `month`, `median` [K m-1], ...; script 11a).
#' @param summ Per-site summaries of script 11a (`site_id`, `layer`, `season`, `share_2_levels`).
#' @return Tibble per month 1-12 and `"all"`: `n_sites`, `median`, `p25`, `p75` [K m-1, over the site medians],
#'   `n_sites_2_levels` (sites whose fit uses only 2 levels on every day, default layer).
network_gamma_summary <- function(monthly, summ) {
  two <- summ$site_id[summ$layer == "default" & summ$season == "all" & summ$share_2_levels >= 1]
  one <- function(g, label) {
    g <- g[!is.na(g)]
    q <- if (length(g)) stats::quantile(g, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3)
    dplyr::tibble(month = label, n_sites = length(g), median = q[2], p25 = q[1], p75 = q[3],
                  n_sites_2_levels = length(unique(two)))
  }
  site_all <- tapply(monthly$median, monthly$site_id, stats::median, na.rm = TRUE)
  dplyr::bind_rows(c(lapply(1:12, function(m) one(monthly$median[monthly$month == m], as.character(m))),
                     list(one(as.numeric(site_all), "all"))))
}

#' Site flags of the network run
#'
#' @param sites Station codes.
#' @param status Run status table (`site_id`, `step`, `state`).
#' @param cases_tbl `02_cases` rows (`site_id`, `case`, `gt5km`).
#' @param coverage `08_coverage` rows (`site_id`, `case`, `set`, `stratum`, `n_rows`, `n_valid_qc_used`).
#' @param cells ERA5 cell metadata (`site_id`, `cell_dist_km`, `dalt_m`).
#' @param steps Steps expected to have run.
#' @param min_coverage Flag below this share of `valid_qc_used` half-hours (PBLHmod, set `etc`) [-].
#' @return Tibble: `site_id`, `flag`, `detail`.
network_flags <- function(sites, status, cases_tbl, coverage, cells, steps = c("01", "11a", "02", "04", "08a"),
                          min_coverage = 0.2) {
  f <- list()
  add <- function(s, flag, detail) f[[length(f) + 1]] <<- dplyr::tibble(site_id = s, flag = flag, detail = detail)
  for (s in sites) {
    st <- status[status$site_id == s, ]
    miss <- setdiff(steps, st$step[st$state %in% "ok"])
    if (length(miss)) {
      failed <- st$step[st$state %in% c("failed", "running")]
      add(s, if (length(failed)) "step_failed" else "step_missing",
          paste0("steps not ok: ", paste(miss, collapse = ", "),
                 if (length(failed)) paste0(" (failed or interrupted: ", paste(failed, collapse = ", "), "; see the log)") else ""))
    }
    big <- cases_tbl[cases_tbl$site_id == s & cases_tbl$case %in% c("default", "clim") & cases_tbl$gt5km > 0, ]
    for (j in seq_len(nrow(big))) add(s, "pblh_gt_5km", paste0(big$case[j], ": ", big$gt5km[j], " half-hours > 5 km"))
    cv <- coverage[coverage$site_id == s & coverage$case == "pblhmod" & coverage$set == "etc" & coverage$stratum == "all", ]
    if (nrow(cv) == 1 && cv$n_valid_qc_used / cv$n_rows < min_coverage) {
      add(s, "low_coverage", sprintf("PBLHmod etc valid_qc_used %.1f %% of half-hours", 100 * cv$n_valid_qc_used / cv$n_rows))
    }
    ce <- cells[cells$site_id == s, ]
    if (nrow(ce) == 1 && abs(ce$dalt_m) > 200) {
      add(s, "era5_orography", sprintf("ERA5 tower-cell orography %+.0f m from the site altitude", -ce$dalt_m))
    }
  }
  if (length(f) == 0) return(dplyr::tibble(site_id = character(), flag = character(), detail = character()))
  dplyr::bind_rows(f)
}
