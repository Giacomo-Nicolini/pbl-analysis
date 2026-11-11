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
# Derived: <derived_dir>/gamma/network_gamma.csv: gamma_med = median of the site medians of the daily ERA5 gamma (pooled
#   median kept for reference), read by step 02
#   (case gamma_med; spec 09 amendment 3); written only when every site has its 11a table and the value changed.
# Figure: 09_network_logratio_<n>sites.png, ln(PBLHmod / ERA5) by site and regime: PBLHmod median and p25-p75 bar,
#   PBLHmod-clim median as open circles on top (spec 09 amendment 2); 09_network_logratio_supp_<n>sites.png, the same
#   with the side-result constants gamma_med and gamma_site (pairs 7, 8; amendment 3).
#
# Usage: Rscript scripts/09_collect.R [--sites-file=config/run_sites.txt] [--gamma-only]
#   --gamma-only: only the network gamma (first pass of the cluster order, after step 11a at every site).

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
