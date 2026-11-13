# Network footprint summary for the paper (spec 15, G 2026-11-12): target-area fraction difference, fetch ratio,
# most affected rows, ICOS site-design criteria (Rebmann et al. 2018) and the ICOS labelling rule (ETC), exceptions;
# gamma supplementary tables. Built by scripts/09_collect.R from the per-half-hour chunk metrics of step 08b (the
# per-site tables of step 08 do not hold the crossed strata nor the clim - ERA5 pair).
#
# Conventions (all tables):
# - Case pairs (difference = first minus second): era5 - pblhmod, pblhmod_clim - pblhmod, era5 - pblhmod_clim.
# - Row set `common_pe`: half-hours valid (validity set x subset) and with an FFP result in every case compared
#   (pblhmod, pblhmod_clim, era5), so the three pairs use the same rows.
# - Target-area fraction (TA) = footprint integral over the ICOS target area relative to the full footprint [%];
#   difference in percentage points [pp]. Fetch: ln(x_80 first / x_80 second) [-], x_80 = distance of the 80 %
#   crosswind-integrated footprint [m].
# - Stability regime = zm/L class of spec 08 (`zl_class()`: unstable < -0.1 <= neutral <= 0.1 < stable); day / night
#   from the input table (`day`); season of the half-hour midpoint.

#' Case pairs of the network footprint summary
#'
#' @return Tibble: `pair` (id), `case_a`, `case_b` (difference = a minus b), `label`.
fp_pairs <- function() {
  dplyr::tibble(pair = c("era5_vs_pblhmod", "pblhmod_clim_vs_pblhmod", "era5_vs_pblhmod_clim"),
                case_a = c("era5", "pblhmod_clim", "era5"), case_b = c("pblhmod", "pblhmod", "pblhmod_clim"),
                label = c("ERA5 - PBLHmod", "PBLHmod-clim - PBLHmod", "ERA5 - PBLHmod-clim"))
}

#' Per-half-hour rows of one site for the network footprint summary
#'
#' Reads the site's footprint input table (`<site>_footprint_inputs.rds`, step 08a) and the per-half-hour chunk
#' metrics (`chunks/<site>_<case>_<month>[a|b]_hh.csv.gz`, step 08b) of the reference year.
#'
#' @param site_id Station code [character].
#' @param year Reference year [integer or character].
#' @param fp_dir `<derived_dir>/footprint` [character].
#' @param cases Footprint cases [character].
#' @param sets,subsets Validity sets and subsets kept as flags [character].
#' @return Tibble, one row per half-hour of the year: `site_id`, `timestamp_start_utc` [POSIXct UTC], `zl_class`
#'   (unstable / neutral / stable / `NA`), `daynight` (day / night), `season`, `regime` (PBLHmod regime), and per case
#'   `<cs>`: `h__<cs>` [m], `ta__<cs>` [0-1, target-area fraction], `x80__<cs>` [m], `ok__<cs>` [logical, FFP result
#'   with a finite target-area fraction]; per set x subset x case `v__<set>__<subset>__<cs>` [logical]. Attribute
#'   `missing`: case-months without chunks.
fp_site_rows <- function(site_id, year, fp_dir, cases = c("pblhmod", "pblhmod_clim", "era5"),
                         sets = c("etc", "kljun_norsl", "kljun"), subsets = c("valid", "valid_qc", "valid_qc_used")) {
  sd <- file.path(fp_dir, site_id)
  ip <- file.path(sd, paste0(site_id, "_footprint_inputs.rds"))
  if (!file.exists(ip)) cli::cli_abort("{site_id}: footprint input table {.path {ip}} not found (step 08a).")
  inp <- readRDS(ip)
  sel <- format(inp$tab$timestamp_start_utc, "%Y", tz = "UTC") == as.character(year)
  if (!any(sel)) cli::cli_abort("{site_id}: no half-hour of {year} in the footprint input table.")
  tab <- inp$tab[sel, ]
  sub <- inp$subsets[sel, ]
  t <- tab$timestamp_start_utc
  out <- dplyr::tibble(site_id = site_id, timestamp_start_utc = t, zl_class = as.character(zl_class(tab$zm_m, tab$L)),
                       daynight = ifelse(tab$day %in% TRUE, "day", ifelse(tab$day %in% FALSE, "night", NA_character_)),
                       season = as.character(season_of(t + 900)), regime = as.character(tab$regime))
  months <- sort(unique(format(t, "%Y-%m", tz = "UTC")))
  cdir <- file.path(sd, "chunks")
  missing <- character()
  for (cs in cases) {
    parts <- unlist(lapply(months, function(mo) {
      p <- chunk_parts_on_disk(cdir, site_id, cs, mo)
      if (length(p) == 0) missing <<- c(missing, paste(cs, mo))
      p
    }))
    files <- file.path(cdir, paste0(site_id, "_", cs, "_", parts, "_hh.csv.gz"))
    files <- files[file.exists(files)]
    x <- if (length(files)) {
      dplyr::bind_rows(lapply(files, function(f) readr::read_csv(
        f, col_select = c("timestamp_start_utc", "h", "x_80", "ffp_ok", "f_target_area"), show_col_types = FALSE,
        progress = FALSE)))
    } else {
      dplyr::tibble(timestamp_start_utc = as.POSIXct(character(), tz = "UTC"), h = numeric(), x_80 = numeric(),
                    ffp_ok = logical(), f_target_area = numeric())
    }
    m <- match(as.numeric(t), as.numeric(x$timestamp_start_utc))
    out[[paste0("h__", cs)]] <- x$h[m]
    out[[paste0("ta__", cs)]] <- x$f_target_area[m]
    out[[paste0("x80__", cs)]] <- x$x_80[m]
    out[[paste0("ok__", cs)]] <- x$ffp_ok[m] %in% TRUE & is.finite(x$f_target_area[m])
    for (s in sets) for (u in subsets) {
      v <- sub[[paste(u, s, cs, sep = "__")]]
      out[[paste("v", s, u, cs, sep = "__")]] <- if (is.null(v)) rep(FALSE, nrow(out)) else v %in% TRUE
    }
  }
  attr(out, "missing") <- missing
  out
}

#' Rows of the `common_pe` row set
#'
#' @param r Rows ([fp_site_rows()], one or more sites).
#' @param set,subset Validity set and subset [character].
#' @param cases Cases that must all be valid with an FFP result [character].
#' @return Logical, one per row.
fp_common <- function(r, set, subset, cases = c("pblhmod", "pblhmod_clim", "era5")) {
  ok <- rep(TRUE, nrow(r))
  for (cs in cases) {
    v <- r[[paste("v", set, subset, cs, sep = "__")]]
    o <- r[[paste0("ok__", cs)]]
    if (is.null(v) || is.null(o)) return(rep(FALSE, nrow(r)))
    ok <- ok & v & o
  }
  ok
}

# Long table of one value per row with the regime x day/night strata, including the margins ("all").
fp_strata_long <- function(site_id, zl, dn, ...) {
  base <- dplyr::tibble(site_id = site_id, zl = zl, dn = dn, ...)
  dplyr::bind_rows(
    dplyr::mutate(base, regime = zl, daynight = dn),
    dplyr::mutate(base, regime = "all", daynight = dn),
    dplyr::mutate(base, regime = zl, daynight = "all"),
    dplyr::mutate(base, regime = "all", daynight = "all")
  ) |> dplyr::filter(!is.na(regime), !is.na(daynight)) |> dplyr::select(-zl, -dn)
}

q_at <- function(x, p) if (length(x)) unname(stats::quantile(x, p, names = FALSE)) else NA_real_

# Per-site, pooled and median-of-site-medians summaries of `value` by regime x day/night.
fp_scope_summary <- function(long, extra = function(d) list()) {
  stat <- function(d, k) dplyr::as_tibble(c(list(n = nrow(d),
                                                 n_sites = if ("site_id" %in% names(d)) length(unique(d$site_id)) else 1L,
                                                 p10 = q_at(d$value, 0.1), median = stats::median(d$value),
                                                 p90 = q_at(d$value, 0.9)), extra(d)))
  site <- dplyr::ungroup(dplyr::group_modify(dplyr::group_by(long, site_id, regime, daynight), stat))
  site$n_sites <- NULL
  pooled <- dplyr::ungroup(dplyr::group_modify(dplyr::group_by(long, regime, daynight), stat))
  sites <- dplyr::summarise(dplyr::group_by(site, regime, daynight), n = sum(n), n_sites = dplyr::n(),
                            p10 = stats::median(p10), median = stats::median(median), p90 = stats::median(p90),
                            .groups = "drop")
  dplyr::bind_rows(dplyr::mutate(site, scope = "site", n_sites = 1L),
                   dplyr::mutate(pooled, scope = "network_pooled", site_id = "network"),
                   dplyr::mutate(sites, scope = "network_site_median", site_id = "network"))
}

fp_regime_order <- function(x) {
  x$regime <- factor(x$regime, levels = c("all", "unstable", "neutral", "stable"))
  x$daynight <- factor(x$daynight, levels = c("all", "day", "night"))
  x <- x[order(x$scope != "site", x$site_id, x$set, x$subset, x$pair, x$regime, x$daynight), ]
  x$regime <- as.character(x$regime)
  x$daynight <- as.character(x$daynight)
  x
}

#' Network table 1: target-area fraction difference between cases
#'
#' @param r Pooled rows of all sites ([fp_site_rows()]).
#' @param sets,subsets Validity sets and subsets [character].
#' @param pairs Output of [fp_pairs()] (pairs whose cases are in `r`).
#' @return Tibble: `scope` (`site`; `network_pooled` = quantiles of all rows of all sites; `network_site_median` =
#'   median over the sites of the site p10, median and p90), `site_id`, `set`, `subset`, `pair`, `case_a`, `case_b`,
#'   `regime`, `daynight` (`all` = margin), `n` [-, half-hours], `n_sites`, `p10_pp`, `median_pp`, `p90_pp` [pp, TA of
#'   a minus TA of b], `median_ta_a_pct`, `median_ta_b_pct` [% , not for `network_site_median`].
fp_target_area_table <- function(r, sets, subsets, pairs = fp_pairs()) {
  out <- list()
  for (s in sets) for (u in subsets) {
    cm <- fp_common(r, s, u)
    if (!any(cm)) next
    for (k in seq_len(nrow(pairs))) {
      a <- pairs$case_a[k]
      b <- pairs$case_b[k]
      ta_a <- 100 * r[[paste0("ta__", a)]][cm]
      ta_b <- 100 * r[[paste0("ta__", b)]][cm]
      long <- fp_strata_long(r$site_id[cm], r$zl_class[cm], r$daynight[cm], value = ta_a - ta_b, ta_a = ta_a, ta_b = ta_b)
      sm <- fp_scope_summary(long, function(d) list(median_ta_a_pct = stats::median(d$ta_a),
                                                   median_ta_b_pct = stats::median(d$ta_b)))
      out[[length(out) + 1]] <- dplyr::mutate(sm, set = s, subset = u, pair = pairs$pair[k], case_a = a, case_b = b)
    }
  }
  x <- dplyr::bind_rows(out)
  if (nrow(x) == 0) return(x)
  x <- dplyr::rename(x, p10_pp = p10, median_pp = median, p90_pp = p90)
  fp_regime_order(x[c("scope", "site_id", "set", "subset", "pair", "case_a", "case_b", "regime", "daynight", "n",
                      "n_sites", "p10_pp", "median_pp", "p90_pp", "median_ta_a_pct", "median_ta_b_pct")])
}

#' Network table 2: fetch, ln(x_80 ratio) between cases and x_80 per case
#'
#' @inheritParams fp_target_area_table
#' @return Tibble as [fp_target_area_table()] with `p10_lnratio`, `median_lnratio`, `p90_lnratio` [-, ln(x_80 a /
#'   x_80 b)] and `median_x80_a_m`, `p90_x80_a_m`, `median_x80_b_m`, `p90_x80_b_m` [m].
fp_fetch_table <- function(r, sets, subsets, pairs = fp_pairs()) {
  out <- list()
  for (s in sets) for (u in subsets) {
    cm <- fp_common(r, s, u)
    if (!any(cm)) next
    for (k in seq_len(nrow(pairs))) {
      a <- pairs$case_a[k]
      b <- pairs$case_b[k]
      xa <- r[[paste0("x80__", a)]][cm]
      xb <- r[[paste0("x80__", b)]][cm]
      lr <- fp_log_ratio(xa, xb)
      ok <- !is.na(lr)
      long <- fp_strata_long(r$site_id[cm][ok], r$zl_class[cm][ok], r$daynight[cm][ok], value = lr[ok], xa = xa[ok],
                             xb = xb[ok])
      sm <- fp_scope_summary(long, function(d) list(median_x80_a_m = stats::median(d$xa), p90_x80_a_m = q_at(d$xa, 0.9),
                                                   median_x80_b_m = stats::median(d$xb), p90_x80_b_m = q_at(d$xb, 0.9)))
      out[[length(out) + 1]] <- dplyr::mutate(sm, set = s, subset = u, pair = pairs$pair[k], case_a = a, case_b = b)
    }
  }
  x <- dplyr::bind_rows(out)
  if (nrow(x) == 0) return(x)
  x <- dplyr::rename(x, p10_lnratio = p10, median_lnratio = median, p90_lnratio = p90)
  fp_regime_order(x[c("scope", "site_id", "set", "subset", "pair", "case_a", "case_b", "regime", "daynight", "n",
                      "n_sites", "p10_lnratio", "median_lnratio", "p90_lnratio", "median_x80_a_m", "p90_x80_a_m",
                      "median_x80_b_m", "p90_x80_b_m")])
}

#' Top `frac` of rows by |target-area difference| (logical selection within each group)
#'
#' @param d Absolute differences [pp].
#' @param group Grouping (e.g. site) [character]; one selection per group.
#' @param frac Share of rows selected [-].
#' @return List: `top` [logical], `threshold` (named by group) [pp]: rows with |d| >= the (1 - frac) quantile.
fp_top_rows <- function(d, group, frac = 0.1) {
  g <- unique(group)
  thr <- stats::setNames(vapply(g, function(k) q_at(d[group == k], 1 - frac), numeric(1)), g)
  list(top = unname(d >= thr[group]), threshold = thr)
}

#' Network table 3: the most affected rows and what removes them
#'
#' Spec 15 A3: per site (and network, pooled), the 10 % of the `common_pe` rows of subset `valid` with the largest
#' |TA difference| of each pair; their shares by stability class, PBLHmod regime, day / night and season, and the
#' shares that fail the `common_pe` rows of `valid_qc` (NEE QC) and of `valid_qc_used` (QC and the FLUXNET NEE used,
#' including the u* filter). `network_site_median` = median of the site shares.
#'
#' @inheritParams fp_target_area_table
#' @param frac Share of rows selected [-].
#' @return Tibble: `scope`, `site_id`, `set`, `pair`, `n_rows`, `n_top`, `threshold_pp`, `share_unstable`,
#'   `share_neutral`, `share_stable` (zm/L class), `share_regime_stable`, `share_regime_unstable_diag`,
#'   `share_regime_convective` (PBLHmod regime), `share_night`, `share_DJF` .. `share_SON`, `share_removed_by_qc`,
#'   `share_removed_by_qc_or_used` [0-1].
fp_topaffected_table <- function(r, sets, pairs = fp_pairs(), frac = 0.1) {
  out <- list()
  shares <- function(i, rows, qc_ok, used_ok) {
    z <- rows$zl_class[i]
    dplyr::tibble(n_top = sum(i), share_unstable = mean(z %in% "unstable"), share_neutral = mean(z %in% "neutral"),
                  share_stable = mean(z %in% "stable"), share_regime_stable = mean(rows$regime[i] %in% "stable"),
                  share_regime_unstable_diag = mean(rows$regime[i] %in% "unstable_diag"),
                  share_regime_convective = mean(rows$regime[i] %in% "convective"),
                  share_night = mean(rows$daynight[i] %in% "night"), share_DJF = mean(rows$season[i] %in% "DJF"),
                  share_MAM = mean(rows$season[i] %in% "MAM"), share_JJA = mean(rows$season[i] %in% "JJA"),
                  share_SON = mean(rows$season[i] %in% "SON"), share_removed_by_qc = mean(!qc_ok[i]),
                  share_removed_by_qc_or_used = mean(!used_ok[i]))
  }
  for (s in sets) {
    cm <- fp_common(r, s, "valid")
    if (!any(cm)) next
    rows <- r[cm, ]
    qc_ok <- fp_common(rows, s, "valid_qc")
    used_ok <- fp_common(rows, s, "valid_qc_used")
    for (k in seq_len(nrow(pairs))) {
      d <- abs(100 * (rows[[paste0("ta__", pairs$case_a[k])]] - rows[[paste0("ta__", pairs$case_b[k])]]))
      st <- fp_top_rows(d, rows$site_id, frac)
      site <- dplyr::bind_rows(lapply(names(st$threshold), function(si) {
        i <- st$top & rows$site_id == si
        dplyr::bind_cols(dplyr::tibble(scope = "site", site_id = si, n_rows = sum(rows$site_id == si),
                                       threshold_pp = unname(st$threshold[si])), shares(i, rows, qc_ok, used_ok))
      }))
      nt <- fp_top_rows(d, rep("network", length(d)), frac)
      pooled <- dplyr::bind_cols(dplyr::tibble(scope = "network_pooled", site_id = "network", n_rows = length(d),
                                               threshold_pp = unname(nt$threshold)), shares(nt$top, rows, qc_ok, used_ok))
      sh <- grep("^share_", names(site), value = TRUE)
      med <- dplyr::bind_cols(dplyr::tibble(scope = "network_site_median", site_id = "network", n_rows = sum(site$n_rows),
                                            threshold_pp = stats::median(site$threshold_pp), n_top = sum(site$n_top)),
                              dplyr::as_tibble(lapply(site[sh], stats::median)))
      out[[length(out) + 1]] <- dplyr::mutate(dplyr::bind_rows(site, pooled, med), set = s, pair = pairs$pair[k],
                                              .after = site_id)
    }
  }
  dplyr::bind_rows(out)
}

#' ICOS site-design criteria (Rebmann et al. 2018) and the ICOS labelling rule (ETC) per half-hour set
#'
#' Rebmann et al. (2018): >= 90 % of the flux from the area of interest, also under stable conditions, for >= 50 % of
#' the time (p. 475); cumulative footprint >= 80 % for an acceptable measurement (p. 487). ICOS labelling rule (ETC;
#' G 2026-11-12): among the half-hours of the "good data" subset, the share with a cumulative target-area
#' contribution >= 70 %; the target is met if that share is >= 70 %.
#'
#' @param f_ta Target-area fractions of the half-hours [0-1].
#' @return One-row tibble: `n`, `share_ta90`, `share_ta80`, `share_ta70` [0-1], `pass_ta90_50pct`, `pass_lab_70_70`
#'   [logical; `NA` without rows].
icos_ta_criteria <- function(f_ta) {
  f <- f_ta[!is.na(f_ta)]
  n <- length(f)
  s <- function(x) if (n) mean(f >= x) else NA_real_
  dplyr::tibble(n = n, share_ta90 = s(0.9), share_ta80 = s(0.8), share_ta70 = s(0.7),
                pass_ta90_50pct = if (n) s(0.9) >= 0.5 else NA, pass_lab_70_70 = if (n) s(0.7) >= 0.7 else NA)
}

#' Network table 4: ICOS criteria per site and case
#'
#' Rebmann criteria on the `common_pe` rows of (`set`, `subset`): all rows and stable rows (zm/L class); ICOS
#' labelling rule on the `common_pe` rows of (`set`, `valid_qc`) ("good data", G 2026-11-12).
#'
#' @param r Pooled rows.
#' @param set,subset Validity set and subset of the Rebmann criteria [character].
#' @param cases Cases [character].
#' @return Tibble: `site_id`, `set`, `subset`, `case`, `n`, `share_ta90`, `pass_ta90`, `n_stable`, `share_ta90_stable`,
#'   `pass_ta90_stable`, `share_ta80` [0-1], `n_lab`, `icos_lab_share`, `icos_lab_pass`, and per site
#'   `ta90_differs`, `ta90_stable_differs`, `icos_lab_differs` [logical: pass / fail differs between cases].
fp_icos_criteria_table <- function(r, set, subset, cases = c("pblhmod", "pblhmod_clim", "era5")) {
  cm <- fp_common(r, set, subset, cases)
  cq <- fp_common(r, set, "valid_qc", cases)
  out <- dplyr::bind_rows(lapply(unique(r$site_id), function(si) dplyr::bind_rows(lapply(cases, function(cs) {
    ta <- r[[paste0("ta__", cs)]]
    i <- cm & r$site_id == si
    a <- icos_ta_criteria(ta[i])
    b <- icos_ta_criteria(ta[i & r$zl_class %in% "stable"])
    l <- icos_ta_criteria(ta[cq & r$site_id == si])
    dplyr::tibble(site_id = si, set = set, subset = subset, case = cs, n = a$n, share_ta90 = a$share_ta90,
                  pass_ta90 = a$pass_ta90_50pct, n_stable = b$n, share_ta90_stable = b$share_ta90,
                  pass_ta90_stable = b$pass_ta90_50pct, share_ta80 = a$share_ta80, n_lab = l$n,
                  icos_lab_share = l$share_ta70, icos_lab_pass = l$pass_lab_70_70)
  }))))
  differs <- function(x) length(unique(x[!is.na(x)])) > 1
  flags <- dplyr::summarise(dplyr::group_by(out, site_id), ta90_differs = differs(pass_ta90),
                            ta90_stable_differs = differs(pass_ta90_stable), icos_lab_differs = differs(icos_lab_pass),
                            .groups = "drop")
  dplyr::left_join(out, flags, by = "site_id")
}

#' Network table 5: exceptions, |p10| or |p90| of the TA difference above a threshold in a regime x day/night stratum
#'
#' The threshold (5 pp) is a reporting choice (spec 15). Only the crossed strata (no `all` margin), per site.
#'
#' @param ta Output of [fp_target_area_table()].
#' @param r Pooled rows (for the median h per case in the stratum).
#' @param threshold_pp Threshold [pp].
#' @return Tibble: `site_id`, `set`, `subset`, `pair`, `regime`, `daynight`, `n`, `p10_pp`, `median_pp`, `p90_pp`,
#'   `smaller_case` (case with the smaller TA: sign of the median difference; `equal` at 0), `median_h_a_m`,
#'   `median_h_b_m` [m], `threshold_pp`.
fp_exceptions_table <- function(ta, r, threshold_pp = 5) {
  x <- if (nrow(ta) == 0) ta else ta[ta$scope == "site" & ta$regime != "all" & ta$daynight != "all" &
                                       (abs(ta$p10_pp) > threshold_pp | abs(ta$p90_pp) > threshold_pp), ]
  if (nrow(x) == 0) {
    return(dplyr::tibble(site_id = character(), set = character(), subset = character(), pair = character(),
                         regime = character(), daynight = character(), n = integer(), p10_pp = numeric(),
                         median_pp = numeric(), p90_pp = numeric(), smaller_case = character(), median_h_a_m = numeric(),
                         median_h_b_m = numeric(), threshold_pp = numeric()))
  }
  h <- t(vapply(seq_len(nrow(x)), function(k) {
    i <- fp_common(r, x$set[k], x$subset[k]) & r$site_id == x$site_id[k] & r$zl_class %in% x$regime[k] &
      r$daynight %in% x$daynight[k]
    c(stats::median(r[[paste0("h__", x$case_a[k])]][i]), stats::median(r[[paste0("h__", x$case_b[k])]][i]))
  }, numeric(2)))
  dplyr::tibble(site_id = x$site_id, set = x$set, subset = x$subset, pair = x$pair, regime = x$regime,
                daynight = x$daynight, n = x$n, p10_pp = x$p10_pp, median_pp = x$median_pp, p90_pp = x$p90_pp,
                smaller_case = ifelse(x$median_pp < 0, x$case_a, ifelse(x$median_pp > 0, x$case_b, "equal")),
                median_h_a_m = h[, 1], median_h_b_m = h[, 2], threshold_pp = threshold_pp)
}

#' Write a table with a `#` comment header
#'
#' @param x Table.
#' @param path Output path [character].
#' @param header Lines of the header (without `#`) [character].
#' @return `path`, invisibly. Read back with `readr::read_csv(path, comment = "#")`.
write_csv_header <- function(x, path, header) {
  writeLines(paste("#", header), path)
  readr::write_csv(x, path, append = TRUE, col_names = TRUE)
  invisible(path)
}

#' Supplementary table: gamma per site
#'
#' @param daily_by_site Named list of daily gamma tables (`date`, `gamma`; [gamma_daily()]).
#' @param igbp Named character (names = site ids): land cover [IGBP].
#' @return Tibble: `site_id`, `igbp`, `gamma_site` [K m-1, median of the valid daily gamma, = [gamma_site_constant()]],
#'   `p25`, `p75`, `p10`, `p90` [K m-1], `n_days`, `n_months` [-, year-months with a valid day], `months_used`
#'   (calendar months present, e.g. `1-12` or `3,4,5`).
gamma_site_table <- function(daily_by_site, igbp = NULL) {
  dplyr::bind_rows(lapply(names(daily_by_site), function(s) {
    d <- daily_by_site[[s]]
    ok <- is.finite(d$gamma)
    g <- d$gamma[ok]
    mo <- as.integer(format(as.Date(d$date[ok]), "%m"))
    cm <- sort(unique(mo))
    mu <- if (length(cm) == 12) "1-12" else paste(cm, collapse = ",")
    dplyr::tibble(site_id = s, igbp = if (is.null(igbp)) NA_character_ else unname(igbp[s]),
                  gamma_site = if (length(g)) stats::median(g) else NA_real_, p25 = q_at(g, 0.25), p75 = q_at(g, 0.75),
                  p10 = q_at(g, 0.1), p90 = q_at(g, 0.9), n_days = length(g),
                  n_months = length(unique(format(as.Date(d$date[ok]), "%Y-%m"))), months_used = mu)
  }))
}

#' Supplementary table: network gamma with a site bootstrap
#'
#' `gamma_med` = median of the site medians ([gamma_network_median()]); 95 % interval from `n_boot` resamples of the
#' sites with replacement (percentile interval), reproducible with `seed`.
#'
#' @param daily_by_site Named list of daily gamma tables.
#' @param n_boot Bootstrap resamples [-].
#' @param seed Random seed [integer].
#' @return One-row tibble: `gamma_med`, `gamma_med_lo95`, `gamma_med_hi95`, `gamma_med_pooled`, `site_median_min`,
#'   `site_median_max` [K m-1], `n_sites`, `n_days`, `n_boot`, `seed`.
gamma_network_table <- function(daily_by_site, n_boot = 2000L, seed = 20261112L) {
  net <- gamma_network_median(daily_by_site)
  med <- vapply(daily_by_site, function(x) {
    g <- x$gamma[is.finite(x$gamma)]
    if (length(g)) stats::median(g) else NA_real_
  }, numeric(1))
  med <- med[!is.na(med)]
  boot <- withr::with_seed(seed, vapply(seq_len(n_boot), function(i) stats::median(sample(med, replace = TRUE)), numeric(1)))
  ci <- unname(stats::quantile(boot, c(0.025, 0.975)))
  dplyr::tibble(gamma_med = net$gamma_med, gamma_med_lo95 = ci[1], gamma_med_hi95 = ci[2],
                gamma_med_pooled = net$gamma_med_pooled, site_median_min = min(med), site_median_max = max(med),
                n_sites = net$n_sites, n_days = net$n_days, n_boot = n_boot, seed = seed)
}

#' Step-08 top-affected shares aggregated over the network (reproduction of the per-site tables)
#'
#' The per-site tables `08_topaffected_<site>_<year>.csv` select the top 10 % by |ln(x_80 ERA5 / x_80 PBLHmod)| on the
#' rows valid in both cases (subset `valid`), shares by PBLHmod regime (not the zm/L class). Aggregated here as the
#' median of the site shares and pooled (weighted by n), for comparison with the spec-15 selection.
#'
#' @param tabs Bound per-site tables with `site_id` ([bind_site_tables()]).
#' @return Tibble in the layout of [fp_topaffected_table()] (rows `group == "top10"`), `selection` =
#'   `abs_lnx80_era5_pblhmod (step 08)`, zm/L shares and seasons `NA`.
fp_topaffected_step08 <- function(tabs) {
  t <- tabs[tabs$group %in% "top10", ]
  if (nrow(t) == 0) return(dplyr::tibble())
  map <- function(x, scope, site) dplyr::tibble(
    scope = scope, site_id = site, set = x$set, pair = "era5_vs_pblhmod", n_rows = NA_integer_, threshold_pp = NA_real_,
    n_top = x$n, share_regime_stable = x$share_stable, share_regime_unstable_diag = x$share_unstable_diag,
    share_regime_convective = x$share_convective, share_night = x$share_night,
    share_removed_by_qc = x$share_removed_by_qc, share_removed_by_qc_or_used = x$share_removed_by_qc_or_used)
  site <- map(t, "site", t$site_id)
  sh <- c("share_regime_stable", "share_regime_unstable_diag", "share_regime_convective", "share_night",
          "share_removed_by_qc", "share_removed_by_qc_or_used")
  agg <- dplyr::bind_rows(lapply(split(site, site$set), function(x) dplyr::bind_rows(
    dplyr::bind_cols(dplyr::tibble(scope = "network_site_median", site_id = "network", set = x$set[1],
                                   pair = "era5_vs_pblhmod", n_top = sum(x$n_top)), dplyr::as_tibble(lapply(x[sh], stats::median))),
    dplyr::bind_cols(dplyr::tibble(scope = "network_pooled", site_id = "network", set = x$set[1], pair = "era5_vs_pblhmod",
                                   n_top = sum(x$n_top)),
                     dplyr::as_tibble(lapply(x[sh], function(v) sum(v * x$n_top) / sum(x$n_top)))))))
  dplyr::mutate(dplyr::bind_rows(site, agg), selection = "abs_lnx80_era5_pblhmod (step 08)")
}

#' The five network footprint tables of spec 15
#'
#' @param rows Pooled rows ([fp_site_rows()] of every site).
#' @param sets,subsets Validity sets and subsets of tables 1, 2 [character] (main: `etc`, `valid_qc_used`).
#' @param top_sets Validity sets of table 3 [character].
#' @param main_set,main_subset Set and subset of tables 4, 5 [character].
#' @param threshold_pp Exception threshold [pp] (reporting choice).
#' @param step08_top Output of [fp_topaffected_step08()] (added to table 3) or `NULL`.
#' @return Named list: `target_area`, `fetch`, `topaffected`, `icos_criteria`, `exceptions`; the ICOS table also for
#'   the sensitivity sets x subsets (column `set`, `subset`).
fp_network_summary <- function(rows, sets = c("etc", "kljun_norsl"), subsets = c("valid", "valid_qc_used"),
                               top_sets = c("etc", "kljun_norsl", "kljun"), main_set = "etc",
                               main_subset = "valid_qc_used", threshold_pp = 5, step08_top = NULL) {
  ta <- fp_target_area_table(rows, sets, subsets)
  fe <- fp_fetch_table(rows, sets, subsets)
  top <- dplyr::mutate(fp_topaffected_table(rows, top_sets), selection = "abs_ta_diff (spec 15)")
  if (!is.null(step08_top) && nrow(step08_top)) top <- dplyr::bind_rows(top, step08_top)
  ic <- dplyr::bind_rows(lapply(sets, function(s) dplyr::bind_rows(lapply(subsets, function(u) {
    fp_icos_criteria_table(rows, s, u)
  }))))
  ex <- fp_exceptions_table(ta[ta$set == main_set & ta$subset == main_subset, ], rows, threshold_pp)
  list(target_area = ta, fetch = fe, topaffected = top, icos_criteria = ic, exceptions = ex)
}

#' Header lines of the network footprint tables (definitions, sign convention, units)
#'
#' @param threshold_pp Exception threshold [pp].
#' @return Named list of character vectors (one per table of [fp_network_summary()]).
fp_table_headers <- function(threshold_pp = 5) {
  common <- c(
    "Spec 15 network footprint summary (scripts/09_collect.R, R/footprint_network.R); reference year per site.",
    "Pairs: difference = first case minus second (era5 - pblhmod, pblhmod_clim - pblhmod, era5 - pblhmod_clim).",
    "Rows: common_pe = half-hours valid (set x subset) and with an FFP result in pblhmod, pblhmod_clim and era5.",
    "Regime = zm/L class (unstable < -0.1 <= neutral <= 0.1 < stable, zm = EC height - 2/3 hc); daynight from the input table; 'all' = margin.",
    "Scopes: site; network_pooled = quantiles of all rows of all sites; network_site_median = median over the sites of the site p10, median, p90.")
  list(
    target_area = c(common, "TA = footprint integral over the ICOS target area / full footprint [%]; p10_pp, median_pp, p90_pp = TA(a) - TA(b) [percentage points]."),
    fetch = c(common, "ln ratio = ln(x_80 a / x_80 b) [-], x_80 = distance of the 80 % crosswind-integrated footprint [m]; x_80 per case: median and p90 [m]."),
    topaffected = c(common,
      "selection 'abs_ta_diff (spec 15)': per site (and network, pooled), the 10 % of the common_pe rows of subset valid with the largest |TA(a) - TA(b)|; threshold_pp = the 90 % quantile.",
      "Shares [0-1] by zm/L class, PBLHmod regime, night, season; share_removed_by_qc = not in common_pe of valid_qc (NEE QC); share_removed_by_qc_or_used = not in common_pe of valid_qc_used (QC and FLUXNET NEE used, u* filter).",
      "selection 'abs_lnx80_era5_pblhmod (step 08)': aggregation of the per-site 08_topaffected tables (top 10 % by |ln x_80 ratio| ERA5 vs PBLHmod, rows valid in both, regime = PBLHmod regime), for comparison."),
    icos_criteria = c(common,
      "Rebmann et al. (2018), Int. Agrophys. 32: >= 90 % of the flux from the area of interest, also under stable conditions, for >= 50 % of the time (p. 475); cumulative footprint >= 80 % for an acceptable measurement (p. 487).",
      "share_ta90 = share of the common_pe rows with TA >= 90 % (all rows; stable = zm/L class stable); pass = share >= 0.5; share_ta80 = share with TA >= 80 %.",
      "ICOS labelling rule (ETC): per site and case, among the half-hours of the 'good data' subset, the share with a cumulative target-area contribution >= 70 %; target met if that share is >= 70 % (G 2026-11-12). Computed on the common_pe rows of valid_qc: icos_lab_share, icos_lab_pass.",
      "*_differs = TRUE where pass / fail differs between the cases at the site."),
    exceptions = c(common,
      sprintf("Exception = |p10| or |p90| of TA(a) - TA(b) above %g pp in a regime x day/night stratum of a site (threshold: a reporting choice, spec 15); main set and subset only.", threshold_pp),
      "smaller_case = case with the smaller TA (sign of the median difference); median_h_a_m, median_h_b_m = median PBLH of the two cases on the stratum rows [m]."))
}

#' Classify the case-months without chunks (fix-15)
#'
#' A month absent from `chunk_list.csv` has no valid row for the case (08c lists only months with rows to run): its
#' missing chunks are expected (e.g. SE-Myc 2025-07). A month listed there but without a chunk on disk is a real gap
#' (task not run, failed or not copied).
#'
#' @param missing Attribute `missing` of [fp_site_rows()] ("<case> <YYYY-MM>") [character].
#' @param site_id Station code [character].
#' @param chunk_list `<derived_dir>/footprint/chunk_list.csv` as read (columns `site_id`, `case`, `start`), or `NULL`
#'   when the file is not available.
#' @return Tibble: `case`, `month` [YYYY-MM], `kind` (`not_listed` = no valid rows, expected; `gap` = listed but
#'   missing; `unknown` = no chunk list).
fp_classify_missing <- function(missing, site_id, chunk_list = NULL) {
  if (length(missing) == 0) return(dplyr::tibble(case = character(), month = character(), kind = character()))
  cs <- sub(" .*$", "", missing)
  mo <- sub("^.* ", "", missing)
  kind <- if (is.null(chunk_list)) {
    rep("unknown", length(missing))
  } else {
    cl <- chunk_list[chunk_list$site_id == site_id, ]
    listed <- paste(cl$case, substr(as.character(cl$start), 1, 7))
    ifelse(paste(cs, mo) %in% listed, "gap", "not_listed")
  }
  dplyr::tibble(case = cs, month = mo, kind = kind)
}

#' Messages for [fp_classify_missing()]
#'
#' @param cm Output of [fp_classify_missing()].
#' @return Named character vector for `cli::cli_inform()`-style output: names `i` (expected), `x` (gaps), `!`
#'   (unknown); empty when nothing is missing.
fp_missing_messages <- function(cm) {
  fmt <- function(k) {
    x <- cm[cm$kind == k, ]
    if (nrow(x) == 0) return(NULL)
    g <- tapply(x$case, x$month, paste, collapse = ", ")
    paste(paste0(names(g), " (", g, ")"), collapse = "; ")
  }
  out <- c(i = if (!is.null(fmt("not_listed"))) paste0("no chunks, as expected (no valid rows: not in chunk_list.csv): ",
                                                       fmt("not_listed")),
           x = if (!is.null(fmt("gap"))) paste0("chunks MISSING although listed in chunk_list.csv (real gap): ",
                                                fmt("gap")),
           `!` = if (!is.null(fmt("unknown"))) paste0("no chunks (chunk_list.csv not found, cannot tell why): ",
                                                      fmt("unknown")))
  out
}
