# Spec 09 A4: network collection of the per-site outputs (steps 01, 11a, 02, 04, 08a; scripts/run_site.R).
# Reads only the small per-site tables, so it runs in seconds on the login node or as a short job.
# Output (<output_dir>/tables/):
#   09_network_inputs.csv          files, reference period, ERA5 period and cell distance per site (01)
#   09_network_gamma_monthly.csv   gamma climatology per site x month (11a) and 09_network_gamma_summary.csv
#                                  (median, p25, p75 over the sites per month and overall; sites with only 2 fit levels)
#   09_network_regimes.csv         PBLHmod, PBLHmod-clim (and gamma_med, gamma_site) per regime: share, quantiles (02)
#   09_network_filters.csv         u* filter and L check counts, validity, gamma per case (02)
#   09_network_era5_metrics.csv    PBLHmod and PBLHmod-clim vs ERA5 (tower cell, bilinear): pairs 1, 2, 6; side-result
#                                  constants gamma_med, gamma_site: pairs 7, 8 (04; amendment 3)
#   09_network_era5_cells.csv      ERA5 cell distance and orography
#   09_network_footprint_coverage.csv  footprint coverage per case x set (08a, full period, stratum all)
#   09_network_nee_qc.csv          FLUXNET NEE QC source (VUT / CUT / none) and status per site and year (08a; amendment 2)
#   09_network_etc_multiplier.csv  ETC z0 multiplier per site: m_loop, m_used, m_source loop / override, approved
#                                  reference year (08a; spec 08 amendment 2)
#   08_reference_years.csv         reference-year proposal (spec 09 A5), with the approved year (approved_year): review,
#                                  then commit config/reference_years.csv
#                                  (exceptions: config/reference_year_exceptions.csv, rule coverage; amendment 2)
#   09_network_run_status.csv      last status per site x step; 09_network_flags.csv sites with problems and the
#                                  stations not run (config/sites_not_run.csv: flag "not_run", spec 09 amendment 1)
#   09_network_gamma_site.csv      gamma_site per site: median of the site's daily ERA5 gamma (11a; amendment 3)
#   09_network_record_start.csv    EC record start = ICOS labelling date (first valid pblhmod_default value, local
#                                  date), last valid value, n_valid (G 2026-11-12; spec 14; the only input read
#                                  outside the per-site tables: time and value columns of the step-04 standard file)
# Derived: <derived_dir>/gamma/network_gamma.csv: gamma_med = median of the site medians of the daily ERA5 gamma (pooled
#   median kept for reference), read by step 02
#   (case gamma_med; spec 09 amendment 3); written only when every site has its 11a table and the value changed.
# Figure: 09_network_logratio_<n>sites.png, ln(PBLHmod / ERA5) by site and regime: PBLHmod median and p25-p75 bar,
#   PBLHmod-clim median as open circles on top (spec 09 amendment 2); 09_network_logratio_supp_<n>sites.png, the same
#   with the side-result constants gamma_med and gamma_site (pairs 7, 8; amendment 3).
#
# Spec 15 (G 2026-11-12), with a '#' header (definitions, signs, units; read with comment = "#"):
#   09_network_fp_target_area.csv  target-area fraction difference [pp] between cases (era5 - pblhmod, pblhmod_clim -
#                                  pblhmod, era5 - pblhmod_clim), per site and network, regime x day/night
#   09_network_fp_fetch.csv        ln(x_80 ratio) between cases and x_80 per case [m], same strata
#   09_network_fp_topaffected.csv  the 10 % rows with the largest |TA difference|: shares and what removes them
#   09_network_fp_icos_criteria.csv  Rebmann et al. (2018) criteria and the ICOS labelling rule (ETC) per site and case
#   09_network_fp_exceptions.csv   strata with |p10| or |p90| of the TA difference > 5 pp
#   09s_gamma_site.csv, 09s_gamma_network.csv  gamma per site (quantiles) and gamma_med with a site bootstrap
#   Figures 09_fp_ta_diff, 09_fp_fetch_ratio, 09_fp_topaffected, 09_fp_icos_criteria (set etc, subset valid_qc_used);
#   supplement 09s_fp_kljun_norsl_* and 09s_fp_valid_*; 09s_gamma_site.
#   The footprint part reads the footprint input tables and hh chunks of the reference years (not only per-site tables).
#
# Usage: Rscript scripts/09_collect.R [--sites-file=config/run_sites.txt] [--gamma-only] [--no-footprint]
#   --gamma-only: only the network gamma (first pass of the cluster order, after step 11a at every site).
#   --no-footprint: skip the spec-15 footprint summary (the rest runs in seconds).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
cfg <- get_cfg()
ids <- trimws(readLines(opt("sites-file", here::here("config", "run_sites.txt")), warn = FALSE))
sites <- ids[nzchar(ids) & !startsWith(ids, "#")]
# Stations not run (spec 09 amendment 1): listed in the flags as "not run (<reason>)", not as missing steps.
not_run <- read_sites_not_run(if (is.null(cfg$sites_not_run_file)) NULL else resolve_path(cfg$sites_not_run_file))
if (length(intersect(sites, not_run$site_id))) {
  cli::cli_alert_warning("In the site list but marked as not run (dropped): {.val {intersect(sites, not_run$site_id)}}.")
}
sites <- setdiff(sites, not_run$site_id)
tdir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
derived <- resolve_path(cfg$derived_dir, must_exist = FALSE)
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
per_site <- function(fmt, dir = tdir) stats::setNames(file.path(dir, sprintf(fmt, sites)), sites)
write_tab <- function(x, name) {
  readr::write_csv(x, file.path(tdir, name))
  cli::cli_alert_success("{name}: {nrow(x)} rows")
}
cli::cli_h1("Network collection: {length(sites)} sites")
if (nrow(not_run)) cli::cli_inform("Not run: {paste0(not_run$site_id, ' (', not_run$reason, ')', collapse = ', ')}.")

# --- network gamma constants (spec 09 amendment 3): gamma_med for step 02, gamma_site per site for the record ---
# Needs the 11a daily tables of every site; written only when all are present and only if the value changed (the
# file is an input of step 02: a rewrite would rerun it).
g_daily_p <- per_site("%s_daily.csv", file.path(derived, "gamma"))
g_miss <- names(g_daily_p)[!file.exists(g_daily_p)]
if (length(g_miss) == length(g_daily_p)) {
  cli::cli_alert_warning("No 11a daily gamma table: network gamma not computed.")
} else {
  g_daily <- lapply(g_daily_p[file.exists(g_daily_p)], function(p) readr::read_csv(p, show_col_types = FALSE, progress = FALSE))
  g_const <- dplyr::bind_rows(lapply(names(g_daily), function(s) dplyr::mutate(gamma_site_constant(g_daily[[s]]), site_id = s, .before = 1)))
  gnet <- gamma_network_median(g_daily)
  write_tab(dplyr::rename(g_const, gamma_site = gamma), "09_network_gamma_site.csv")
  cli::cli_inform("Network gamma: median of the site medians at {gnet$n_sites} sites = {signif(gnet$gamma_med, 4)} K m-1 (reference: median of the {gnet$n_days} pooled daily values {signif(gnet$gamma_med_pooled, 4)}).")
  if (length(g_miss)) {
    cli::cli_alert_warning("11a daily gamma missing at {.val {g_miss}}: {.path {gamma_network_path(cfg)}} not written (gamma_med needs every site).")
  } else if (write_gamma_network(gnet, gamma_network_path(cfg))) {
    cli::cli_alert_success("{.path {gamma_network_path(cfg)}} written: step 02 reruns at every site (case gamma_med).")
  } else {
    cli::cli_inform("{.path {gamma_network_path(cfg)}} unchanged.")
  }
}
if ("--gamma-only" %in% args) quit(save = "no", status = 0)

# --- run status ---
status <- dplyr::bind_rows(lapply(sites, function(s) dplyr::bind_rows(lapply(names(pipeline_steps()), function(st) {
  x <- read_run_status(run_status_path(s, st, cfg))
  if (is.null(x)) return(NULL)
  dplyr::tibble(site_id = s, step = st, state = x$state, start = x$start, end = x$end, exit_code = x$exit_code,
                git_commit = x$git_commit)
}))))
if (nrow(status) == 0) status <- dplyr::tibble(site_id = character(), step = character(), state = character())
write_tab(status, "09_network_run_status.csv")

# --- 01 inputs (long -> wide) ---
inp <- bind_site_tables(per_site("01_inputs_%s.csv"))
if (nrow(inp)) {
  inp_w <- stats::reshape(as.data.frame(inp), idvar = "site_id", timevar = "item", direction = "wide")
  names(inp_w) <- sub("^value\\.", "", names(inp_w))
  write_tab(dplyr::as_tibble(inp_w), "09_network_inputs.csv")
}

# --- 11a gamma ---
g_month <- bind_site_tables(per_site("%s_monthly.csv", file.path(derived, "gamma")))
g_summ <- bind_site_tables(per_site("11a_gamma_summary_%s.csv"))
if (nrow(g_month)) {
  write_tab(g_month, "09_network_gamma_monthly.csv")
  write_tab(network_gamma_summary(g_month, g_summ), "09_network_gamma_summary.csv")
}

# --- 02 regimes and filters ---
reg <- bind_site_tables(per_site("02_regimes_%s.csv"))
if (nrow(reg)) write_tab(reg[reg$case %in% c("default", "clim", "gamma_med", "gamma_site"), ], "09_network_regimes.csv")
cases_tbl <- bind_site_tables(per_site("02_cases_%s.csv"))
if (nrow(cases_tbl)) write_tab(cases_tbl, "09_network_filters.csv")

# --- 04 metrics vs ERA5, cells ---
met <- bind_site_tables(per_site("04_metrics_%s.csv"))
pairs_main <- c("1_pblhmod_default_vs_era5_center", "2_pblhmod_default_vs_era5_bilinear", "6_pblhmod_clim_vs_era5_center",
                "7_pblhmod_gamma_med_vs_era5_center", "8_pblhmod_gamma_site_vs_era5_center")
if (nrow(met)) {
  met <- met[met$pair %in% pairs_main & (met$stratum == "all" | startsWith(met$stratum, "regime=")), ]
  write_tab(met, "09_network_era5_metrics.csv")
}
cm_path <- file.path(derived, "era5", "era5_cell_metadata.csv")
cells <- if (file.exists(cm_path)) readr::read_csv(cm_path, show_col_types = FALSE, progress = FALSE) else dplyr::tibble()
if (nrow(cells)) {
  cells <- cells[cells$site_id %in% sites, ]
  write_tab(cells, "09_network_era5_cells.csv")
}

# --- 08a coverage and reference years ---
cov <- bind_site_tables(per_site("08_coverage_%s.csv"))
if (nrow(cov)) write_tab(cov[cov$stratum == "all", ], "09_network_footprint_coverage.csv")
nq <- bind_site_tables(per_site("08_nee_qc_%s.csv"))
if (nrow(nq)) write_tab(nq, "09_network_nee_qc.csv")
yc <- bind_site_tables(per_site("08_year_counts_%s.csv"))
# Approved reference years (config/reference_years.csv), if committed; the step-08 outputs of a year are not read here
ry_path <- resolve_path(cfg$reference_years_file, must_exist = FALSE)
ry_ok <- if (file.exists(ry_path)) read_reference_years(ry_path) else NULL
# ETC z0 multiplier per site: loop value, value used, loop / override (spec 08 amendment 2)
em <- network_etc_multiplier(bind_site_tables(per_site("08_etc_multiplier_%s.csv")), ry_ok)
if (nrow(em)) {
  write_tab(em, "09_network_etc_multiplier.csv")
  for (i in which(em$m_source == "override")) {
    cli::cli_inform("{em$site_id[i]}: ETC multiplier override m_used = {em$m_used[i]} (m_loop = {em$m_loop[i]}).")
  }
}
# Exceptions to the reference-year rule (spec 09 amendment 2): used only when the standard rule finds no year
rex <- read_reference_year_exceptions(if (is.null(cfg$reference_year_exceptions_file)) NULL else
  resolve_path(cfg$reference_year_exceptions_file, must_exist = FALSE))
if (nrow(yc)) {
  ry <- dplyr::bind_rows(lapply(split(yc, yc$site_id), function(x) {
    e <- rex[rex$site_id == x$site_id[1], ]
    p <- if (nrow(e)) propose_reference_year(x, rule = e$rule, rule_reason = e$reason) else propose_reference_year(x)
    dplyr::bind_cols(dplyr::tibble(site_id = x$site_id[1]), p)
  }))
  # approved year next to the proposal (config/reference_years.csv; e.g. DE-Tha 2025, spec 08 amendment 2)
  ry$approved_year <- if (is.null(ry_ok)) NA_integer_ else ry_ok$year[match(ry$site_id, ry_ok$site_id)]
  write_tab(ry, "08_reference_years.csv")
  cli::cli_inform("Reference-year proposal: review {.path {file.path(tdir, '08_reference_years.csv')}}, then commit the approved list as {.file config/reference_years.csv} (site_id, year, reason).")
}

# --- EC record start = ICOS labelling date (G 2026-11-12): first valid pblhmod_default value (standard file of step
# 04; only its time and value columns are read), for the supplement ---
st_sites <- read_sites(resolve_path(cfg$sites_file))
std_root <- resolve_path(cfg$standard_dir, must_exist = FALSE)
rs <- network_record_start(stats::setNames(file.path(std_root, sites, paste0(sites, "__pblhmod_default.csv.gz")), sites),
                           stats::setNames(st_sites$utc_offset_h[match(sites, st_sites$site_id)], sites))
write_tab(rs, "09_network_record_start.csv")
if (any(is.na(rs$n_valid))) cli::cli_alert_warning("No pblhmod_default standard file at {.val {rs$site_id[is.na(rs$n_valid)]}}.")

# --- spec 15: network footprint summary (reference year per site; per-half-hour chunk metrics) ---
# Reads every site's footprint input table and hh chunks of the reference year (a few minutes); --no-footprint
# skips it. Tables with a '#' header (read with comment = "#").
fp_out <- NULL
if (!"--no-footprint" %in% args) {
  if (is.null(ry_ok)) {
    cli::cli_alert_warning("No approved reference years: network footprint summary skipped.")
  } else {
    cli::cli_h2("Network footprint summary (spec 15)")
    fp_dir <- file.path(derived, "footprint")
    fp_sites <- intersect(sites, ry_ok$site_id)
    rows <- list()
    # fix-15: a missing month is expected when it is not in chunk_list.csv (no valid rows), a gap when it is listed
    cl_path <- file.path(fp_dir, "chunk_list.csv")
    chunk_list <- if (file.exists(cl_path)) readr::read_csv(cl_path, col_types = readr::cols(.default = "c"), progress = FALSE)
    n_gap <- 0L
    for (s in fp_sites) {
      y <- ry_ok$year[ry_ok$site_id == s]
      x <- tryCatch(fp_site_rows(s, y, fp_dir), error = function(e) e)
      if (inherits(x, "error")) {
        cli::cli_alert_danger("{s} {y}: {conditionMessage(x)}")
        next
      }
      if (length(attr(x, "missing"))) {
        cm <- fp_classify_missing(attr(x, "missing"), s, chunk_list)
        n_gap <- n_gap + sum(cm$kind == "gap")
        msg <- fp_missing_messages(cm)
        for (k in names(msg)) {
          f <- switch(k, i = cli::cli_alert_info, x = cli::cli_alert_danger, cli::cli_alert_warning)
          f(paste0(s, " ", y, ": ", msg[[k]]))
        }
      }
      rows[[s]] <- x
    }
    rows <- dplyr::bind_rows(rows)
    if (n_gap > 0) cli::cli_alert_danger("{n_gap} listed case-month{?s} without chunks (real gaps): rerun those chunks (08c --pending) before using the summary.")
    if (nrow(rows)) {
      st08 <- fp_topaffected_step08(bind_site_tables(stats::setNames(
        file.path(tdir, paste0("08_topaffected_", fp_sites, "_", ry_ok$year[match(fp_sites, ry_ok$site_id)], ".csv")), fp_sites)))
      fp_out <- fp_network_summary(rows, step08_top = st08)
      hd <- fp_table_headers()
      for (k in names(fp_out)) {
        write_csv_header(fp_out[[k]], file.path(tdir, paste0("09_network_fp_", k, ".csv")), hd[[k]])
        cli::cli_alert_success("09_network_fp_{k}.csv: {nrow(fp_out[[k]])} rows")
      }
      # figures: main set etc / valid_qc_used; supplement kljun_norsl / valid_qc_used and etc / valid
      igbp <- st_sites$igbp[match(fp_sites, st_sites$site_id)]
      cm <- fp_common(rows, "etc", "valid_qc_used")
      mh <- tapply(rows$h__pblhmod[cm], rows$site_id[cm], stats::median)[fp_sites]
      site_order <- fp_sites[order(igbp, mh)] # land cover, then median PBLHmod h (main rows)
      ex_sites <- unique(fp_out$exceptions$site_id)
      variants <- list(list(prefix = "09_fp_", set = "etc", subset = "valid_qc_used"),
                       list(prefix = "09s_fp_kljun_norsl_", set = "kljun_norsl", subset = "valid_qc_used"),
                       list(prefix = "09s_fp_valid_", set = "etc", subset = "valid"))
      # fix-15: short subtitle, details in a wrapped caption (titles were clipped); heights per site row
      lab_of <- function(v) paste0("Validity set ", v$set, ", subset ", v$subset, "; reference year per site")
      cap_of <- function(extra = NULL) paste(c("Sites by land cover (IGBP), then median PBLHmod h.", extra,
                                               fp_exception_caption()), collapse = " ")
      n_ax <- length(site_order) + 1 # sites + network row
      for (v in variants) {
        ta_v <- fp_out$target_area[fp_out$target_area$set == v$set & fp_out$target_area$subset == v$subset, ]
        fe_v <- fp_out$fetch[fp_out$fetch$set == v$set & fp_out$fetch$subset == v$subset, ]
        if (nrow(ta_v)) {
          p <- plot_fp_site_diff(ta_v, "p10_pp", "median_pp", "p90_pp", site_order, ex_sites,
                                 "Target-area fraction difference [pp] (p10-p90, median)", lab_of(v), cap_of())
          save_figure(p, paste0(v$prefix, "ta_diff"), 22, fp_fig_height(n_ax, 2))
          p <- plot_fp_site_diff(fe_v, "p10_lnratio", "median_lnratio", "p90_lnratio", site_order, ex_sites,
                                 "ln(x_80 ratio) [-] (p10-p90, median)", lab_of(v), cap_of())
          save_figure(p, paste0(v$prefix, "fetch_ratio"), 22, fp_fig_height(n_ax, 2))
        }
        ic_v <- fp_out$icos_criteria[fp_out$icos_criteria$set == v$set & fp_out$icos_criteria$subset == v$subset, ]
        if (nrow(ic_v)) {
          p <- plot_fp_icos_criteria(ic_v, site_order, paste0(lab_of(v), "; labelling rule on valid_qc"), ex_sites,
                                     cap_of("Cases dodged vertically within each site row (PBLHmod on top)."))
          save_figure(p, paste0(v$prefix, "icos_criteria"), 22, fp_fig_height(length(site_order), 1, 0.45))
        }
      }
      for (s in c("etc", "kljun_norsl")) {
        tp <- fp_out$topaffected[fp_out$topaffected$selection == "abs_ta_diff (spec 15)" & fp_out$topaffected$set == s &
                                   fp_out$topaffected$pair == "era5_vs_pblhmod", ]
        if (nrow(tp)) {
          p <- plot_fp_topaffected(tp, site_order, paste0("Top 10 % |TA difference| ERA5 - PBLHmod; set ", s, ", subset valid"),
                                   paste("Sites by land cover (IGBP), then median PBLHmod h; network bar = pooled",
                                         "shares (sums to 100 %); median of the site shares in 09_network_fp_topaffected."))
          save_figure(p, if (s == "etc") "09_fp_topaffected" else "09s_fp_kljun_norsl_topaffected", 22,
                      fp_fig_height(n_ax, 1, 0.38))
        }
      }
      cli::cli_alert_success("Figures 09_fp_* and 09s_fp_* written.")
    }
  }
}

# --- spec 15 C: gamma supplementary tables (needs the 11a daily tables of the sites) ---
if (exists("g_daily") && length(g_daily)) {
  gst <- gamma_site_table(g_daily, stats::setNames(st_sites$igbp, st_sites$site_id))
  gnt <- gamma_network_table(g_daily)
  write_csv_header(gst, file.path(tdir, "09s_gamma_site.csv"), c(
    "Daily morning gamma from ERA5 (scripts/11a_gamma_era5.R, 500-3500 m above the site) per site [K m-1].",
    "gamma_site = median of the valid daily values (side result, spec 09 amendment 3); p25, p75, p10, p90 of the daily values;",
    "n_days = valid days, n_months = year-months with a valid day, months_used = calendar months present; igbp = land cover."))
  write_csv_header(gnt, file.path(tdir, "09s_gamma_network.csv"), c(
    "gamma_med = median of the site medians [K m-1] (recommended operational EC-only value, decision 3b, G 2026-11-11);",
    "gamma_med_lo95 / hi95 = 2.5 / 97.5 % of the medians of n_boot resamples of the sites with replacement (seed given);",
    "gamma_med_pooled = median of all site-days; site_median_min / max = range of the site medians."))
  save_figure(plot_gamma_site(gst, gnt), "09s_gamma_site", 16, 5 + 0.3 * nrow(gst))
  cli::cli_alert_success("09s_gamma_site.csv, 09s_gamma_network.csv (gamma_med {signif(gnt$gamma_med, 4)} [{signif(gnt$gamma_med_lo95, 3)}, {signif(gnt$gamma_med_hi95, 3)}]).")
}

# --- flags ---
fl <- network_flags(sites, status, if (nrow(cases_tbl)) cases_tbl else dplyr::tibble(site_id = character(), case = character(), gt5km = integer()),
                    if (nrow(cov)) cov else dplyr::tibble(site_id = character(), case = character(), set = character(), stratum = character(),
                                                          n_rows = integer(), n_valid_qc_used = integer()),
                    if (nrow(cells)) cells else dplyr::tibble(site_id = character(), cell_dist_km = numeric(), dalt_m = numeric()),
                    not_run = not_run, nee_qc = if (nrow(nq)) nq else NULL)
write_tab(fl, "09_network_flags.csv")
if (nrow(fl)) print(as.data.frame(fl), row.names = FALSE)

# --- figure: median ln(PBLHmod / ERA5) by site and regime ---
# PBLHmod: median and p25-p75 bar; PBLHmod-clim: open circles on top (spec 09 amendment 2)
if (nrow(met)) {
  d <- network_logratio_table(met)
  no_iqr <- unique(d$site_id[d$source == "pblhmod_default" & is.na(d$p25)])
  if (length(no_iqr)) {
    cli::cli_alert_warning("No p25/p75 log ratio (step-04 metrics older than spec 09 amendment 2; re-run 04) at {.val {no_iqr}}: no bar there.")
  }
  p <- plot_network_logratio(d)
  pth <- save_figure(p, paste0("09_network_logratio_", length(unique(d$site_id)), "sites"), 22,
                     6 + 0.35 * length(unique(d$site_id)))
  cli::cli_alert_success("{.path {pth[['png']]}}")
}
# Supplementary version (spec 09 amendment 3): also the side-result constants gamma_med and gamma_site (pairs 7, 8)
if (nrow(met)) {
  ds <- network_logratio_table(met, pairs = network_logratio_pairs(supp = TRUE))
  for (s in c("pblhmod_gamma_med", "pblhmod_gamma_site")) {
    miss <- setdiff(unique(d$site_id), ds$site_id[ds$source == s])
    if (length(miss)) cli::cli_alert_warning("{s}: no step-04 pair at {length(miss)} site{?s} ({.val {miss}}; gamma_med needs the second pass, see spec 09 amendment 3).")
  }
  p <- plot_network_logratio(ds)
  pth <- save_figure(p, paste0("09_network_logratio_supp_", length(unique(ds$site_id)), "sites"), 22,
                     6 + 0.35 * length(unique(ds$site_id)))
  cli::cli_alert_success("{.path {pth[['png']]}}")
}
