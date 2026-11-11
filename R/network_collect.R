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
#' `n_valid_qc` counts the same half-hours in the `valid_qc` subset (NEE QC good in FLUXES, before the FLUXNET u*
#' filter): the data coverage of the `coverage` exception rule of [propose_reference_year()] (spec 09 amendment 2).
#'
#' @param tab Input table of `scripts/08a_footprint_flags.R` (`timestamp_start_utc`, `h_measured` if any).
#' @param subsets Subset columns (`valid_qc_used__etc__pblhmod`, `valid_qc_used__etc__era5`, `valid_qc__etc__pblhmod`,
#'   `valid_qc__etc__era5`).
#' @return Tibble per year: `year`, `n_usable` [half-hours], `n_months` [months with at least one usable half-hour],
#'   `n_measured` [usable half-hours with a measured h], `n_valid_qc` [half-hours], `n_months_valid_qc` [months].
footprint_year_counts <- function(tab, subsets) {
  get <- function(k) { v <- subsets[[k]]; if (is.null(v)) rep(FALSE, nrow(tab)) else v }
  ok <- get("valid_qc_used__etc__pblhmod") & get("valid_qc_used__etc__era5")
  qc <- get("valid_qc__etc__pblhmod") & get("valid_qc__etc__era5")
  meas <- if ("h_measured" %in% names(tab)) !is.na(tab$h_measured) else rep(FALSE, nrow(tab))
  yr <- as.integer(format(tab$timestamp_start_utc, "%Y"))
  mo <- format(tab$timestamp_start_utc, "%m")
  dplyr::bind_rows(lapply(sort(unique(yr)), function(y) {
    i <- yr == y
    dplyr::tibble(year = y, n_usable = sum(ok & i), n_months = length(unique(mo[ok & i])), n_measured = sum(ok & i & meas),
                  n_valid_qc = sum(qc & i), n_months_valid_qc = length(unique(mo[qc & i])))
  }))
}

#' Read the reference-year exceptions
#'
#' Spec 09 amendment 2 (G 2026-11-03): stations whose reference year may follow another rule when the standard rule
#' finds no year (e.g. GL-ZaF, remote site with large data loss, kept in the footprint run).
#'
#' @param path CSV with `site_id`, `rule` (`"coverage"`: most `valid_qc` half-hours, any number of months), `reason`,
#'   `decided`; a missing file means no exception.
#' @return Tibble `site_id`, `rule`, `reason`, `decided` (character).
read_reference_year_exceptions <- function(path) {
  cols <- c("site_id", "rule", "reason", "decided")
  if (is.null(path) || !file.exists(path)) return(dplyr::tibble(!!!stats::setNames(rep(list(character()), 4), cols)))
  x <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), comment = "#",
                       show_col_types = FALSE, progress = FALSE)
  miss <- setdiff(cols, names(x))
  if (length(miss)) cli::cli_abort("{.path {path}}: missing column{?s} {.field {miss}}.")
  if (anyNA(x$site_id) || anyDuplicated(x$site_id)) cli::cli_abort("{.path {path}}: empty or duplicated {.field site_id}.")
  bad <- setdiff(x$rule, reference_year_rules())
  if (length(bad)) cli::cli_abort("{.path {path}}: unknown rule{?s} {.val {bad}} (known: {.val {reference_year_rules()}}).")
  if (anyNA(x$reason)) cli::cli_abort("{.path {path}}: every exception needs a {.field reason}.")
  x[cols]
}

# Rules of propose_reference_year(): "standard" (spec 09 A5) and the exception "coverage" (spec 09 amendment 2).
reference_year_rules <- function() c("standard", "coverage")

#' Reference-year proposal of one site
#'
#' Spec 09 A5 (spec 08 section 4c), amended by spec 09 amendment 3 (G 2026-11-04): among the years with at least
#' `min_months` months of data, a year overlapping the site's measured source is preferred when its usable half-hours
#' are at least `measured_min_share` of the best year's (most measured half-hours; ties: most recent); otherwise the
#' year with the most usable half-hours, years within `tol` of the maximum being equivalent (most recent).
#'
#' Exception rule `coverage` (spec 09 amendment 2, G 2026-11-03; stations listed in
#' `config/reference_year_exceptions.csv`): only when the standard rule finds no year, the year with the most
#' `valid_qc` half-hours (`n_valid_qc`, any number of months; ties: most recent), the reason saying so.
#'
#' @param counts Output of [footprint_year_counts()].
#' @param min_months Minimum months with usable data [-].
#' @param tol Relative tolerance on the number of usable half-hours [-].
#' @param measured_min_share Minimum usable half-hours of a measured-overlap year, relative to the best year [-]
#'   (spec 09 amendment 3: 0.8).
#' @param rule `"standard"` or `"coverage"` (see [reference_year_rules()]).
#' @param rule_reason Why the site has the exception (written into `reason`) [character].
#' @return One-row tibble: `year` (`NA` if no year qualifies), `n_usable`, `n_months`, `n_measured`, `n_candidates`
#'   (measured-overlap years within `measured_min_share`, or years within `tol`), `reason`, `rule` (rule that chose
#'   the year).
propose_reference_year <- function(counts, min_months = 9, tol = 0.02, measured_min_share = 0.8, rule = "standard",
                                   rule_reason = NA_character_) {
  rule <- match.arg(rule, reference_year_rules())
  ok <- counts[counts$n_months >= min_months & counts$n_usable > 0, ]
  if (nrow(ok) == 0) {
    none <- paste0("no year with >= ", min_months, " months of usable data")
    out <- dplyr::tibble(year = NA_integer_, n_usable = NA_integer_, n_months = NA_integer_, n_measured = NA_integer_,
                         n_candidates = 0L, reason = none, rule = "standard")
    if (rule != "coverage") return(out)
    if (!"n_valid_qc" %in% names(counts)) {
      out$reason <- paste0(none, "; coverage exception: no n_valid_qc in the year counts (re-run step 08a)")
      return(out)
    }
    cv <- counts[counts$n_valid_qc > 0, ]
    if (nrow(cv) == 0) {
      out$reason <- paste0(none, "; coverage exception: no valid_qc half-hours")
      return(out)
    }
    pick <- cv[order(-cv$n_valid_qc, -cv$year), ][1, ]
    return(dplyr::tibble(year = as.integer(pick$year), n_usable = as.integer(pick$n_usable),
                         n_months = as.integer(pick$n_months), n_measured = as.integer(pick$n_measured),
                         n_candidates = 1L,
                         reason = paste0(none, "; exception (config/reference_year_exceptions.csv): most valid_qc ",
                                         "half-hours (", pick$n_valid_qc, ", ", pick$n_months_valid_qc,
                                         " months), any number of months: ", rule_reason),
                         rule = "coverage"))
  }
  best <- max(ok$n_usable)
  meas <- ok[ok$n_measured > 0 & ok$n_usable >= measured_min_share * best, ]
  if (nrow(meas)) {
    cand <- meas
    pick <- cand[order(-cand$n_measured, -cand$year), ][1, ]
    why <- sprintf("overlapping the measured source: most measured half-hours (%d), usable %.0f %% of the best year (>= %.0f %%)",
                   as.integer(pick$n_measured), 100 * pick$n_usable / best, 100 * measured_min_share)
  } else {
    cand <- ok[ok$n_usable >= (1 - tol) * best, ]
    pick <- cand[order(-cand$year), ][1, ]
    why <- if (nrow(cand) > 1) sprintf("most usable half-hours (within %g %%), most recent", 100 * tol) else "most usable half-hours"
  }
  dplyr::tibble(year = as.integer(pick$year), n_usable = as.integer(pick$n_usable), n_months = as.integer(pick$n_months),
                n_measured = as.integer(pick$n_measured), n_candidates = nrow(cand), reason = why, rule = "standard")
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

#' Read the list of stations that are not run
#'
#' Spec 09 amendment 1: stations kept in `config/sites.csv` (metadata) but left out of the network run, with the
#' reason (e.g. IT-SR2: no data in the L2I 2026 release, G 2026-11-03).
#'
#' @param path CSV with `site_id`, `reason` (short label, e.g. `"no L2 data"`), `detail`, `decided`; a missing file
#'   means no station is excluded.
#' @return Tibble `site_id`, `reason`, `detail`, `decided` (character).
read_sites_not_run <- function(path) {
  cols <- c("site_id", "reason", "detail", "decided")
  if (is.null(path) || !file.exists(path)) return(dplyr::tibble(!!!stats::setNames(rep(list(character()), 4), cols)))
  x <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()), comment = "#",
                       show_col_types = FALSE, progress = FALSE)
  miss <- setdiff(cols, names(x))
  if (length(miss)) cli::cli_abort("{.path {path}}: missing column{?s} {.field {miss}}.")
  if (anyNA(x$site_id) || anyDuplicated(x$site_id)) cli::cli_abort("{.path {path}}: empty or duplicated {.field site_id}.")
  if (anyNA(x$reason)) cli::cli_abort("{.path {path}}: every station needs a {.field reason}.")
  x[cols]
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
#' @param not_run Output of [read_sites_not_run()]: listed as `"not_run"` (detail `"not run (<reason>): ..."`),
#'   never as a missing step, even if also in `sites`.
#' @param nee_qc `08_nee_qc` rows (`site_id`, `year`, `nee_qc_source`, `status`, ...; [nee_qc_status()]; the
#'   whole-period row is used): flag `"nee_qc_cut"` when the source is `"CUT"` (`NEE_VUT_REF_QC` entirely missing) and
#'   `"nee_qc"` when the status is not `"ok"` (`valid_qc_used` empty; spec 09 amendment 2).
#' @return Tibble: `site_id`, `flag`, `detail`.
network_flags <- function(sites, status, cases_tbl, coverage, cells, steps = c("01", "11a", "02", "04", "08a"),
                          min_coverage = 0.2, not_run = NULL, nee_qc = NULL) {
  f <- list()
  add <- function(s, flag, detail) f[[length(f) + 1]] <<- dplyr::tibble(site_id = s, flag = flag, detail = detail)
  if (!is.null(not_run)) {
    for (j in seq_len(nrow(not_run))) {
      add(not_run$site_id[j], "not_run", paste0("not run (", not_run$reason[j], ")",
                                                if (!is.na(not_run$detail[j])) paste0(": ", not_run$detail[j]) else "",
                                                if (!is.na(not_run$decided[j])) paste0(" [", not_run$decided[j], "]") else ""))
    }
    sites <- setdiff(sites, not_run$site_id)
  }
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
    nq <- if (is.null(nee_qc)) NULL else nee_qc[nee_qc$site_id == s, ]
    if (!is.null(nq) && "year" %in% names(nq)) nq <- nq[is.na(nq$year), ] # whole-period row
    src <- if (!is.null(nq) && "nee_qc_source" %in% names(nq)) nq$nee_qc_source else rep("VUT", NROW(nq))
    if (!is.null(nq) && nrow(nq) == 1 && src == "CUT") {
      add(s, "nee_qc_cut", sprintf("NEE_VUT_REF_QC entirely missing: valid_qc_used from NEE_CUT_REF_QC (%d FLUXNET rows, %d QC 0, %d used)",
                                   as.integer(nq$n_fluxnet), as.integer(nq$n_qc0), as.integer(nq$n_used)))
    }
    if (!is.null(nq) && nrow(nq) == 1 && nq$status != "ok") {
      add(s, "nee_qc", sprintf("FLUXNET NEE_%s_REF_QC %s: %d FLUXNET rows, %d QC values, %d QC 0 (valid_qc_used empty)",
                               if (src == "none") "VUT/CUT" else src, nq$status, as.integer(nq$n_fluxnet),
                               as.integer(nq$n_qc), as.integer(nq$n_qc0)))
    }
    ce <- cells[cells$site_id == s, ]
    if (nrow(ce) == 1 && abs(ce$dalt_m) > 200) {
      add(s, "era5_orography", sprintf("ERA5 tower-cell orography %+.0f m from the site altitude", -ce$dalt_m))
    }
  }
  if (length(f) == 0) return(dplyr::tibble(site_id = character(), flag = character(), detail = character()))
  dplyr::bind_rows(f)
}

#' Step-04 pairs of the network log-ratio figures
#'
#' @param supp `FALSE`: main figure (PBLHmod, PBLHmod-clim); `TRUE`: supplementary figure, also the side-result
#'   constants `pblhmod_gamma_med` and `pblhmod_gamma_site` (spec 09 amendment 3).
#' @return Named character: source id -> pair.
network_logratio_pairs <- function(supp = FALSE) {
  p <- c(pblhmod_default = "1_pblhmod_default_vs_era5_center", pblhmod_clim = "6_pblhmod_clim_vs_era5_center")
  if (supp) p <- c(p, pblhmod_gamma_med = "7_pblhmod_gamma_med_vs_era5_center",
                   pblhmod_gamma_site = "8_pblhmod_gamma_site_vs_era5_center")
  p
}

#' Network log-ratio table for the figure
#'
#' Median and interquartile range of ln(PBLHmod / ERA5 tower cell) per site x regime, PBLHmod (pair 1) and
#' PBLHmod-clim (pair 6) (spec 09 amendment 2, G 2026-11-03); for the supplementary figure also the side-result
#' constants gamma_med (pair 7) and gamma_site (pair 8) (spec 09 amendment 3, G 2026-11-04).
#'
#' @param met Step-04 metrics, long (`site_id`, `pair`, `stratum`, `metric`, `value`); `p25_log_ratio` and
#'   `p75_log_ratio` are `NA` when the metrics predate spec 09 amendment 2.
#' @param regimes Regimes to keep (`"all"` = whole table).
#' @param pairs Named character: source id -> step-04 pair ([network_logratio_pairs()]).
#' @return Tibble: `site_id`, `source` (names of `pairs`), `regime`, `median`, `p25`, `p75` [-, natural log],
#'   `n_log` [-].
network_logratio_table <- function(met, regimes = c("all", "stable", "unstable_diag", "convective"),
                                   pairs = network_logratio_pairs()) {
  d <- met[met$pair %in% pairs & (met$stratum == "all" | startsWith(met$stratum, "regime=")), ]
  d$regime <- ifelse(d$stratum == "all", "all", sub("^regime=", "", d$stratum))
  d <- d[d$regime %in% regimes, ]
  key <- unique(d[c("site_id", "pair", "regime")])
  get <- function(k) {
    x <- d[d$metric == k, ]
    x$value[match(paste(key$site_id, key$pair, key$regime), paste(x$site_id, x$pair, x$regime))]
  }
  dplyr::tibble(site_id = key$site_id, source = names(pairs)[match(key$pair, pairs)], regime = key$regime,
                median = get("median_log_ratio"), p25 = get("p25_log_ratio"), p75 = get("p75_log_ratio"),
                n_log = get("n_log"))
}

#' Network table of the ETC z0 multiplier
#'
#' Spec 08 amendment 2: per site the loop value, the value used and its source (`loop` or `override`,
#' `config/etc_z0_multiplier.csv`), from the 08a tables `08_etc_multiplier_<site>.csv`. Tables written before the
#' amendment have only `m_etc` (always the loop value then): `m_loop = m_used = m_etc`, `m_source = "loop"`.
#'
#' @param x Bound 08a tables ([bind_site_tables()]).
#' @param ref_years Optional approved reference years ([read_reference_years()]): adds `reference_year`.
#' @return Tibble `site_id`, `m_loop`, `m_used` [-], `m_source`, `n_rows`, `n_rejected` (with `m_used`),
#'   `n_rejected_loop` [-], `zm_over_z0_median`, `zm_over_z0_min` [-], `override_reason`, `reference_year`.
network_etc_multiplier <- function(x, ref_years = NULL) {
  if (!nrow(x)) return(dplyr::tibble())
  fill <- function(k, v) if (k %in% names(x)) x[[k]] else rep(v, nrow(x))
  m_loop <- fill("m_loop", NA_real_)
  old <- is.na(m_loop)
  m_loop[old] <- x$m_etc[old]
  m_used <- fill("m_used", NA_real_)
  m_used[is.na(m_used)] <- x$m_etc[is.na(m_used)]
  src <- fill("m_source", NA_character_)
  src[is.na(src)] <- "loop"
  nrl <- fill("n_rejected_loop", NA_real_)
  nrl[old] <- x$n_rejected[old]
  out <- dplyr::tibble(site_id = x$site_id, m_loop = m_loop, m_used = m_used, m_source = src, n_rows = x$n_rows,
                       n_rejected = x$n_rejected, n_rejected_loop = nrl, zm_over_z0_median = x$zm_over_z0_median,
                       zm_over_z0_min = fill("zm_over_z0_min", NA_real_),
                       override_reason = fill("override_reason", NA_character_))
  out$reference_year <- if (is.null(ref_years)) NA_integer_ else ref_years$year[match(out$site_id, ref_years$site_id)]
  out
}
