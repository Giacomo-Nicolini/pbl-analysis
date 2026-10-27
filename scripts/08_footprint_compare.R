# Spec 08 (phase 5a), section 4: footprint sensitivity to the PBLH source, analysis of the 08b chunks for one station
# and period. Paired comparisons on common rows (valid in both cases, FFP computed for both); coverage separately
# (08a; here only the monthly figure from the full-period input table).
# Output (<period> = e.g. _2022-09; test outputs labelled "test, one month"):
#   tables  08_halfhour_metrics_<site><period>.csv   ln(x_r case / x_r ref), elasticity, target-area difference, by stratum
#           08_topaffected_<site><period>.csv        top 10 % / 25 % |ln x_80 ratio| (ERA5 vs PBLHmod) and validity changes
#           08_climatology_<site><period>.csv        monthly (and annual, from the monthly sums) climatologies: target
#                                                    area, 80 % area / distance, Jaccard
#           08_target_area_<site><period>.csv        target-area fraction (p10/median/p90, % of the full footprint) by
#                                                    stability class, difference to PBLHmod [pp], 90 % region at the edge
#           08_rebmann_<site><period>.csv            ICOS target-area criteria (Rebmann et al. 2018) per case
#           08_anchor_<site><period>.csv             +-20 % h vs Kljun et al. (2015) Table 5
#   figures 08_coverage_<site>, 08_dx80_by_zmh_<site><period>, 08_dx80_diurnal_<site><period>,
#           08_clim_contours_<site><period>, 08_topaffected_<site><period>
# ERA5 tower cell only in figures; era5_bilinear in the tables.
#
# Usage: Rscript scripts/08_footprint_compare.R SITE [--start=YYYY-MM-DD --end=YYYY-MM-DD]

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
site_id <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code.")

cfg <- get_cfg()
fp_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint", site_id)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
t0 <- if (!is.null(opt("start"))) as.POSIXct(opt("start"), tz = "UTC") else NULL
t1 <- if (!is.null(opt("end"))) as.POSIXct(opt("end"), tz = "UTC") + 86400 else NULL
per <- if (is.null(t0)) "" else paste0("_", format(t0, "%Y-%m"), if (format(t0, "%Y-%m") != format(t1 - 1, "%Y-%m")) paste0("_", format(t1 - 1, "%Y-%m")) else "")
label <- if (per == "") site_id else paste0(site_id, ", ", sub("^_", "", per), " (test, one month)")
inp <- readRDS(file.path(fp_dir, paste0(site_id, "_footprint_inputs", per, ".rds")))
tab <- inp$tab
sub_cols <- inp$subsets
mask <- readRDS(file.path(fp_dir, paste0(site_id, "_target_area_mask.rds")))
grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 3)
cli::cli_h1("Footprint comparison: {label}")

# ---- Chunks ----
cases <- sub("^h_", "", grep("^h_", names(tab), value = TRUE))
cases <- setdiff(cases, "oper")
months <- sort(unique(format(tab$timestamp_start_utc, "%Y-%m")))
chunk <- function(cs, mo, what) file.path(fp_dir, "chunks", paste0(site_id, "_", cs, "_", mo, "_", what))
hh <- list()
clim <- list()
for (cs in cases) {
  f_hh <- chunk(cs, months, "hh.csv.gz")
  have <- file.exists(f_hh) & file.exists(chunk(cs, months, "meta.rds"))
  if (!all(have)) cli::cli_warn("{cs}: chunks missing for {paste(months[!have], collapse = ', ')}.")
  if (!any(have)) next
  x <- dplyr::bind_rows(lapply(f_hh[have], readr::read_csv, show_col_types = FALSE, progress = FALSE))
  hh[[cs]] <- x[match(as.numeric(tab$timestamp_start_utc), as.numeric(x$timestamp_start_utc)), ]
  for (mo in months[have]) {
    for (cl in readRDS(chunk(cs, mo, "clim.rds"))) clim[[length(clim) + 1]] <- c(list(case = cs, month = mo), cl)
  }
}
cases <- names(hh)
cli::cli_inform("Cases with chunks: {paste(cases, collapse = ', ')}.")
if (!"pblhmod" %in% cases) cli::cli_abort("No PBLHmod chunks for {label}: run scripts/08b_footprint_run.R first.")
sets <- c("kljun", "kljun_norsl", "etc", "none")
sets <- sets[vapply(sets, function(s) any(sub_cols[[paste("valid", s, "pblhmod", sep = "__")]]), logical(1))]
subsets <- c("valid", "valid_qc", "valid_qc_used", "valid_qc_h")
qty <- c("x_peak", "x_50", "x_70", "x_80", "x_90")

# ---- Strata ----
t_mid <- tab$timestamp_start_utc + 900
strata <- list(
  all = rep("all", nrow(tab)),
  regime = as.character(tab$regime),
  zl_class = as.character(zl_class(tab$zm_m, tab$L)),
  zmh_class = as.character(zmh_class(tab$zm_m, tab$h_pblhmod)),
  season = as.character(season_of(t_mid)),
  hour = sprintf("%02d", hour_bin_centred(t_mid)),
  daynight = ifelse(tab$day %in% TRUE, "day", ifelse(tab$day %in% FALSE, "night", NA))
)
valid_of <- function(u, s, cs) {
  v <- sub_cols[[paste(u, s, cs, sep = "__")]]
  if (is.null(v)) rep(FALSE, nrow(tab)) else v
}
ffp_ok <- function(cs) hh[[cs]]$ffp_ok %in% TRUE

# Paired summary of case `cs` against reference `ref` on rows valid in both (and in `also`).
pair_summary <- function(cs, ref, s, u, also = NULL) {
  common <- valid_of(u, s, cs) & valid_of(u, s, ref) & ffp_ok(cs) & ffp_ok(ref)
  for (a in also) common <- common & valid_of(u, s, a) & ffp_ok(a)
  dta <- 100 * (hh[[cs]]$f_target_area - hh[[ref]]$f_target_area)
  out <- list()
  for (q in qty) {
    lr <- fp_log_ratio(hh[[cs]][[q]], hh[[ref]][[q]])
    el <- fp_elasticity(hh[[cs]][[q]], hh[[ref]][[q]], hh[[cs]]$h, hh[[ref]]$h)
    for (sn in names(strata)) {
      g <- strata[[sn]]
      for (lv in sort(unique(g[common & !is.na(g)]))) {
        i <- common & g %in% lv
        out[[length(out) + 1]] <- dplyr::bind_cols(
          dplyr::tibble(case = cs, reference = ref, set = s, subset = u, quantity = q, stratum = sn, level = lv),
          fp_summary(lr[i], el[i], if (q == "x_80") dta[i] else NULL))
      }
    }
  }
  dplyr::bind_rows(out)
}

# ---- Half-hourly differences ----
hm <- list()
for (s in sets) for (u in subsets) {
  for (cs in setdiff(cases, "pblhmod")) hm[[length(hm) + 1]] <- pair_summary(cs, "pblhmod", s, u)
  if ("measured" %in% cases) {
    # FI-Hyy: PBLHmod and ERA5 vs the measured MLH on rows with an exact MLH (exact values over-represent shallow
    # layers: spec 07)
    for (cs in intersect(c("pblhmod", "era5"), cases)) hm[[length(hm) + 1]] <- pair_summary(cs, "measured", s, u, also = setdiff(c("pblhmod", "era5"), cs))
  }
}
hm <- dplyr::bind_rows(hm)
readr::write_csv(hm, file.path(tables_dir, paste0("08_halfhour_metrics_", site_id, per, ".csv")))
cli::cli_h2("ln(x_80 case / x_80 ref), all rows, subset valid")
print_tbl(hm[hm$quantity == "x_80" & hm$stratum == "all" & hm$subset == "valid",
             c("case", "reference", "set", "n", "median_lr", "p25_lr", "p75_lr", "share_gt10", "share_gt25", "median_elasticity", "median_dta_pp", "p90_abs_dta_pp")])
cli::cli_h2("PBLHmod vs ERA5 by regime, x_80, set none: valid vs valid_qc_used")
print_tbl(hm[hm$quantity == "x_80" & hm$stratum == "regime" & hm$case == "era5" & hm$reference == "pblhmod" & hm$set == "none" & hm$subset %in% c("valid", "valid_qc_used"),
             c("subset", "level", "n", "median_lr", "p25_lr", "p75_lr", "share_gt10", "share_gt25", "median_dta_pp")])

# ---- Which half-hours matter (ERA5 vs PBLHmod) ----
ta <- list()
if ("era5" %in% cases) for (s in sets) {
  common <- valid_of("valid", s, "era5") & valid_of("valid", s, "pblhmod") & ffp_ok("era5") & ffp_ok("pblhmod")
  alr <- abs(fp_log_ratio(hh$era5$x_80, hh$pblhmod$x_80))
  thr <- stats::quantile(alr[common], c(0.75, 0.9), na.rm = TRUE)
  grp <- list(all = common, top25 = common & alr >= thr[1], top10 = common & alr >= thr[2],
              valid_pblhmod_only = valid_of("valid", s, "pblhmod") & !valid_of("valid", s, "era5"),
              valid_era5_only = valid_of("valid", s, "era5") & !valid_of("valid", s, "pblhmod"))
  for (gn in names(grp)) {
    i <- grp[[gn]] %in% TRUE
    both_qc <- valid_of("valid_qc", s, "pblhmod") | valid_of("valid_qc", s, "era5")
    both_used <- valid_of("valid_qc_used", s, "pblhmod") | valid_of("valid_qc_used", s, "era5")
    both_h <- valid_of("valid_qc_h", s, "pblhmod") | valid_of("valid_qc_h", s, "era5")
    reg <- table(factor(tab$regime[i], levels = c("stable", "unstable_diag", "convective")), useNA = "no")
    ta[[length(ta) + 1]] <- dplyr::tibble(set = s, group = gn, n = sum(i),
      abs_lr_threshold = if (gn == "top25") thr[[1]] else if (gn == "top10") thr[[2]] else NA_real_,
      share_stable = reg[["stable"]] / sum(i), share_unstable_diag = reg[["unstable_diag"]] / sum(i), share_convective = reg[["convective"]] / sum(i),
      share_night = mean(strata$daynight[i] %in% "night"),
      share_removed_by_qc = mean(!both_qc[i]), share_removed_by_qc_or_used = mean(!both_used[i]),
      share_removed_by_h_qc = mean(!both_h[i]))
  }
}
ta <- dplyr::bind_rows(ta)
readr::write_csv(ta, file.path(tables_dir, paste0("08_topaffected_", site_id, per, ".csv")))
cli::cli_h2("Top-affected half-hours (|ln x_80 ratio|, ERA5 vs PBLHmod) and validity changes")
print_tbl(ta)

# ---- Climatologies ----
# Monthly: from the chunks. Annual (years with more than one month in the period): from the monthly summed fields
# (sum of sums / sum of counts, smoothed at the end; ffp_clim_from_sums()), as calc_footprint_FFP_climatology() would
# give for all half-hours of the year.
for (yr in unique(substr(months, 1, 4))) {
  mo_y <- months[substr(months, 1, 4) == yr]
  if (length(mo_y) < 2) next
  for (cs in cases) {
    f_s <- chunk(cs, mo_y, "sums.rds")
    if (!all(file.exists(f_s))) { cli::cli_warn("{cs} {yr}: monthly sums missing, no annual climatology."); next }
    parts <- lapply(f_s, ffp_sums_read)
    gk <- unique(unlist(lapply(parts, function(p) names(p$sums))))
    for (key0 in unique(sub("\\|(day|night)\\|", "|", gk))) {
      p <- strsplit(key0, "\\|")[[1]]
      for (dn in c("day", "night", "all")) {
        k <- paste(p[1], p[2], if (dn == "all") c("day", "night") else dn, p[3], sep = "|")
        res <- ffp_clim_from_sums(parts, k, grid, mask)
        if (is.null(res)) next
        clim[[length(clim) + 1]] <- list(case = cs, month = yr, set = p[1], subset = p[2], daynight = dn, rowset = p[3],
                                         n = res$n, metrics = res$metrics, region80 = res$region)
      }
    }
    rm(parts)
  }
}
cl_tab <- dplyr::bind_rows(lapply(clim, function(c0) dplyr::bind_cols(
  dplyr::tibble(case = c0$case, period = c0$month, period_type = if (nchar(c0$month) == 4) "year" else "month",
                set = c0$set, subset = c0$subset, daynight = c0$daynight,
                rowset = c0$rowset, n = c0$n, share80_in_target_area = region_share_inside(c0$region80, mask)),
  c0$metrics)))
key <- function(c0) paste(c0$month, c0$set, c0$subset, c0$daynight, c0$rowset, sep = "|")
ref_regions <- list()
for (c0 in clim) if (c0$case == "pblhmod") ref_regions[[key(c0)]] <- c0$region80
cl_tab$jaccard80_vs_pblhmod <- vapply(clim, function(c0) {
  r0 <- ref_regions[[key(c0)]]
  if (c0$case == "pblhmod" || is.null(r0) || c0$rowset == "own") NA_real_ else jaccard_index(c0$region80, r0)
}, numeric(1))
readr::write_csv(cl_tab, file.path(tables_dir, paste0("08_climatology_", site_id, per, ".csv")))
cli::cli_h2("Monthly climatologies, common rows (PBLHmod, ERA5{if ('measured' %in% cases) ', measured'}), set none")
print_tbl(cl_tab[cl_tab$set == "none" & cl_tab$rowset != "own" & cl_tab$case %in% c("pblhmod", "era5", "measured"),
                 c("case", "period", "subset", "daynight", "rowset", "n", "f_target_area", "share80_in_target_area", "area_80", "xmax_80", "jaccard80_vs_pblhmod")])

# ---- Target-area fraction distribution ----
# f_target_area = integral of the (smoothed) footprint over the target area, relative to the FULL footprint (the FFP
# function integrates to 1 over the plane; f_domain = part inside the +-1000 m domain). edge_90 = the 90 % source
# region reaches the domain edge (includes f_domain < 0.9), the analogue of the code's flag_err = 2.
q3 <- function(v) if (any(!is.na(v))) unname(stats::quantile(v, c(0.1, 0.5, 0.9), na.rm = TRUE)) else rep(NA_real_, 3)
tg <- list()
for (cs in cases) for (s in sets) for (u in subsets) {
  own <- valid_of(u, s, cs) & ffp_ok(cs)
  common <- own & valid_of(u, s, "pblhmod") & ffp_ok("pblhmod")
  dta <- 100 * (hh[[cs]]$f_target_area - hh$pblhmod$f_target_area)
  for (sn in c("all", "zl_class", "daynight")) {
    g <- strata[[sn]]
    for (lv in sort(unique(g[own & !is.na(g)]))) {
      i <- own & g %in% lv
      j <- common & g %in% lv
      qa <- q3(100 * hh[[cs]]$f_target_area[i])
      qd <- if (cs == "pblhmod") rep(NA_real_, 3) else q3(dta[j])
      tg[[length(tg) + 1]] <- dplyr::tibble(case = cs, set = s, subset = u, stratum = sn, level = lv, n = sum(i),
        p10_ta_pct = qa[1], median_ta_pct = qa[2], p90_ta_pct = qa[3],
        median_f_domain_pct = 100 * stats::median(hh[[cs]]$f_domain[i]),
        share_edge_90 = mean(hh[[cs]]$edge_90[i] %in% TRUE), share_f_domain_lt90 = mean(hh[[cs]]$f_domain[i] < 0.9),
        n_common = if (cs == "pblhmod") NA_integer_ else sum(j),
        p10_dta_pp = qd[1], median_dta_pp = qd[2], p90_dta_pp = qd[3])
    }
  }
}
tg <- dplyr::bind_rows(tg)
readr::write_csv(tg, file.path(tables_dir, paste0("08_target_area_", site_id, per, ".csv")))
cli::cli_h2("Target-area fraction [% of the full footprint] by stability class, set none, subset valid")
print_tbl(tg[tg$set == "none" & tg$subset == "valid" & tg$stratum %in% c("all", "zl_class") & tg$case %in% c("pblhmod", "era5", "measured"),
             c("case", "level", "n", "p10_ta_pct", "median_ta_pct", "p90_ta_pct", "median_f_domain_pct", "share_edge_90",
               "p10_dta_pp", "median_dta_pp", "p90_dta_pp")])

# ---- ICOS target-area criteria (Rebmann et al. 2018), per half-hour, on rows common with PBLHmod ----
# A site-design criterion (Rebmann2018): >= 90 % of the flux from the area of interest also under stable conditions for
# >= 50 % of the time (p. 475); cumulative footprint >= 80 % for an acceptable measurement (p. 487).
rb <- list()
for (cs in setdiff(cases, "pblhmod")) for (s in sets) for (u in subsets) {
  common <- valid_of(u, s, cs) & valid_of(u, s, "pblhmod") & ffp_ok(cs) & ffp_ok("pblhmod")
  for (gn in c("all", "night", "stable")) {
    i <- common & switch(gn, all = TRUE, night = strata$daynight %in% "night", stable = strata$zl_class %in% "stable")
    a <- rebmann_criteria(hh[[cs]]$f_target_area[i])
    b <- rebmann_criteria(hh$pblhmod$f_target_area[i])
    rb[[length(rb) + 1]] <- dplyr::tibble(case = cs, set = s, subset = u, rows = gn, n = a$n,
      share_ta90_case = a$share_ta90, share_ta90_pblhmod = b$share_ta90, share_cif80_case = a$share_cif80,
      share_cif80_pblhmod = b$share_cif80, pass_ta90_case = a$pass_ta90_50pct, pass_ta90_pblhmod = b$pass_ta90_50pct,
      verdict_differs = a$pass_ta90_50pct != b$pass_ta90_50pct)
  }
}
rb <- dplyr::bind_rows(rb)
readr::write_csv(rb, file.path(tables_dir, paste0("08_rebmann_", site_id, per, ".csv")))
cli::cli_h2("Rebmann et al. (2018) target-area criteria, subset valid, set none")
print_tbl(rb[rb$subset == "valid" & rb$set == "none", ])

# ---- Anchor: +-20 % h vs Kljun et al. (2015) Table 5 ----
an <- list()
for (cs in intersect(c("pblhmod_x0.8", "pblhmod_x1.2"), cases)) for (s in sets) {
  common <- valid_of("valid", s, cs) & valid_of("valid", s, "pblhmod") & ffp_ok(cs) & ffp_ok("pblhmod")
  pc <- fp_peak_change(hh[[cs]]$x_peak, hh$pblhmod$x_peak)
  for (zc in c("unstable", "neutral", "stable")) {
    i <- common & strata$zl_class %in% zc
    an[[length(an) + 1]] <- dplyr::tibble(case = cs, set = s, zl_class = zc, n = sum(i),
      median_dxmax_pct = stats::median(pc$dxmax_pct[i]), p90_abs_dxmax_pct = unname(stats::quantile(abs(pc$dxmax_pct[i]), 0.9, na.rm = TRUE)),
      median_dfymax_pct = stats::median(pc$dfymax_pct[i]))
  }
}
an <- dplyr::left_join(dplyr::bind_rows(an), kljun2015_table5_h20(), by = "zl_class")
readr::write_csv(an, file.path(tables_dir, paste0("08_anchor_", site_id, per, ".csv")))
cli::cli_h2("Anchor: h x 0.8 / x 1.2 vs Kljun et al. (2015) Table 5 (|change| for h +-20 %)")
print_tbl(an[an$set == "none", ])

# ---- Figures ----
pal <- c(pblhmod = "#0072B2", era5 = "#D55E00", era5_bilinear = "#E69F00", measured = "#009E73", bound = "#CC79A7",
         "pblhmod_x0.8" = "#56B4E9", "pblhmod_x1.2" = "#000000")
fig_sub <- function(x) ggplot2::labs(subtitle = paste0(label, if (nzchar(x)) paste0("; ", x) else ""))

# Coverage by month (full-period input table)
full_file <- file.path(fp_dir, paste0(site_id, "_footprint_inputs.rds"))
if (file.exists(full_file)) {
  fi <- readRDS(full_file)
  mo <- format(fi$tab$timestamp_start_utc, "%Y-%m")
  cv <- list()
  for (cs in intersect(c("pblhmod", "era5", "measured"), sub("^h_", "", grep("^h_", names(fi$tab), value = TRUE)))) for (s in setdiff(sets, "kljun")) for (u in subsets) {
    v <- fi$subsets[[paste(u, s, cs, sep = "__")]]
    cv[[length(cv) + 1]] <- dplyr::tibble(case = cs, set = s, subset = u, month = as.Date(paste0(names(tapply(v, mo, sum)), "-01")),
                                          share = as.numeric(tapply(v, mo, mean)))
  }
  cv <- dplyr::bind_rows(cv)
  p <- ggplot2::ggplot(cv, ggplot2::aes(month, 100 * share, colour = case)) + ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::facet_grid(set ~ subset) + ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = NULL, y = "Valid half-hours [% of month]", colour = NULL, title = "Footprint coverage per PBLH case",
                  subtitle = paste0(site_id, ", full period")) + theme_pub()
  save_figure(p, paste0("08_coverage_", site_id), 18, 14)
}

if ("era5" %in% cases) {
  common <- valid_of("valid", "none", "era5") & valid_of("valid", "none", "pblhmod") & ffp_ok("era5") & ffp_ok("pblhmod")
  dd <- dplyr::tibble(zmh = tab$zm_m / tab$h_pblhmod, lr = fp_log_ratio(hh$era5$x_80, hh$pblhmod$x_80),
                      regime = tab$regime, qc_used = valid_of("valid_qc_used", "none", "era5") & valid_of("valid_qc_used", "none", "pblhmod"),
                      qc_h = valid_of("valid_qc_h", "none", "era5") & valid_of("valid_qc_h", "none", "pblhmod"),
                      hour = hour_bin_centred(t_mid), season = season_of(t_mid))[common, ]
  p <- ggplot2::ggplot(dd, ggplot2::aes(zmh, lr, colour = regime)) + ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_point(size = 0.7, alpha = 0.6) + ggplot2::scale_x_log10() +
    ggplot2::labs(x = "zm / h (PBLHmod)", y = "ln(x80 ERA5 / x80 PBLHmod)", colour = NULL,
                  title = "Footprint 80 % distance: ERA5 vs PBLHmod h") + fig_sub("set none, subset valid") + theme_pub()
  save_figure(p, paste0("08_dx80_by_zmh_", site_id, per))

  dl <- dplyr::bind_rows(
    dplyr::summarise(dplyr::group_by(dd, season, hour), med = stats::median(lr), n = dplyr::n(), .groups = "drop") |> dplyr::mutate(subset = "valid"),
    dplyr::summarise(dplyr::group_by(dd[dd$qc_used, ], season, hour), med = stats::median(lr), n = dplyr::n(), .groups = "drop") |> dplyr::mutate(subset = "valid_qc_used"),
    dplyr::summarise(dplyr::group_by(dd[dd$qc_h, ], season, hour), med = stats::median(lr), n = dplyr::n(), .groups = "drop") |> dplyr::mutate(subset = "valid_qc_h"))
  p <- ggplot2::ggplot(dl, ggplot2::aes(hour, med, colour = subset)) + ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) + ggplot2::facet_wrap(~season) +
    ggplot2::scale_x_continuous(breaks = seq(0, 21, 3)) +
    ggplot2::labs(x = "Hour (UTC, centred bins)", y = "Median ln(x80 ERA5 / x80 PBLHmod)", colour = NULL,
                  title = "Diurnal cycle of the footprint difference") + fig_sub("set none") + theme_pub()
  save_figure(p, paste0("08_dx80_diurnal_", site_id, per))

  tp <- ta[ta$set == "none" & ta$group %in% c("all", "top25", "top10"), ]
  wc <- c("share_stable", "share_unstable_diag", "share_convective", "share_night", "share_removed_by_qc", "share_removed_by_qc_or_used", "share_removed_by_h_qc")
  tl <- dplyr::bind_rows(lapply(wc, function(w) dplyr::tibble(group = tp$group, what = w, share = tp[[w]])))
  tl$what <- factor(sub("^share_", "", tl$what), levels = c("stable", "unstable_diag", "convective", "night", "removed_by_qc", "removed_by_qc_or_used", "removed_by_h_qc"))
  tl$group <- factor(tl$group, levels = c("all", "top25", "top10"))
  p <- ggplot2::ggplot(tl, ggplot2::aes(what, 100 * share, fill = group)) + ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = c(all = "grey60", top25 = "#E69F00", top10 = "#D55E00")) +
    ggplot2::labs(x = NULL, y = "Share of half-hours [%]", fill = NULL,
                  title = "Half-hours most affected by the PBLH source (|ln x80 ratio|, ERA5 vs PBLHmod)") +
    fig_sub("set none, subset valid") + theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  save_figure(p, paste0("08_topaffected_", site_id, per))

  # 80 % source regions of the monthly climatologies (common rows), over the target area; 5 m resolution for drawing
  ds <- 5
  iy <- seq(1, length(grid$y), by = ds)
  ix <- seq(1, length(grid$x), by = ds)
  as_df <- function(m, what) dplyr::tibble(x = rep(grid$x[ix], each = length(iy)), y = rep(grid$y[iy], times = length(ix)), z = c(m[iy, ix] * 1), what = what)
  cc <- list(as_df(mask, "target area"))
  for (c0 in clim) {
    if (!(c0$case %in% c("pblhmod", "era5") && c0$set == "none" && c0$rowset == "common_pe" && c0$subset %in% c("valid", "valid_qc_used", "valid_qc_h") && nchar(c0$month) == 7 && c0$daynight %in% c("day", "night"))) next
    m <- matrix(FALSE, length(grid$y), length(grid$x))
    m[c0$region80] <- TRUE
    cc[[length(cc) + 1]] <- dplyr::mutate(as_df(m, c0$case), subset = c0$subset, daynight = c0$daynight)
  }
  ta_df <- cc[[1]]
  rg <- dplyr::bind_rows(cc[-1])
  if (nrow(rg)) {
    panels <- unique(rg[c("subset", "daynight")])
    ta_df <- dplyr::bind_rows(lapply(seq_len(nrow(panels)), function(k) dplyr::mutate(ta_df, subset = panels$subset[k], daynight = panels$daynight[k])))
    p <- ggplot2::ggplot() + ggplot2::geom_contour(data = ta_df, ggplot2::aes(x, y, z = z), breaks = 0.5, colour = "grey40", linetype = "dashed") +
      ggplot2::geom_contour(data = rg, ggplot2::aes(x, y, z = z, colour = what), breaks = 0.5) +
      ggplot2::facet_grid(daynight ~ subset) + ggplot2::coord_equal() + ggplot2::scale_colour_manual(values = pal) +
      ggplot2::labs(x = "x [m, east]", y = "y [m, north]", colour = NULL,
                    title = "80 % source region, monthly climatology", caption = "Dashed: ICOS target area") +
      fig_sub("set none, common rows") + theme_pub()
    save_figure(p, paste0("08_clim_contours_", site_id, per), 20, 14)
  }
}
cli::cli_inform("Tables and figures written for {label}.")
