# Spec 09 A4: network collection of the per-site outputs (steps 01, 11a, 02, 04, 08a; scripts/run_site.R).
# Reads only the small per-site tables, so it runs in seconds on the login node or as a short job.
# Output (<output_dir>/tables/):
#   09_network_inputs.csv          files, reference period, ERA5 period and cell distance per site (01)
#   09_network_gamma_monthly.csv   gamma climatology per site x month (11a) and 09_network_gamma_summary.csv
#                                  (median, p25, p75 over the sites per month and overall; sites with only 2 fit levels)
#   09_network_regimes.csv         PBLHmod and PBLHmod-clim per regime: share and quantiles (02)
#   09_network_filters.csv         u* filter and L check counts, validity, gamma per case (02)
#   09_network_era5_metrics.csv    PBLHmod and PBLHmod-clim vs ERA5 (tower cell, bilinear): pairs 1, 2, 6 (04)
#   09_network_era5_cells.csv      ERA5 cell distance and orography
#   09_network_footprint_coverage.csv  footprint coverage per case x set (08a, full period, stratum all)
#   09_network_nee_qc.csv          FLUXNET NEE QC source (VUT / CUT / none) and status per site and year (08a; amendment 2)
#   08_reference_years.csv         reference-year proposal (spec 09 A5): review, then commit config/reference_years.csv
#                                  (exceptions: config/reference_year_exceptions.csv, rule coverage; amendment 2)
#   09_network_run_status.csv      last status per site x step; 09_network_flags.csv sites with problems and the
#                                  stations not run (config/sites_not_run.csv: flag "not_run", spec 09 amendment 1)
# Figure: 09_network_logratio_<n>sites.png, ln(PBLHmod / ERA5) by site and regime: PBLHmod median and p25-p75 bar,
#   PBLHmod-clim median as open circles on top (spec 09 amendment 2).
#
# Usage: Rscript scripts/09_collect.R [--sites-file=config/run_sites.txt]

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
if (nrow(reg)) write_tab(reg[reg$case %in% c("default", "clim"), ], "09_network_regimes.csv")
cases_tbl <- bind_site_tables(per_site("02_cases_%s.csv"))
if (nrow(cases_tbl)) write_tab(cases_tbl, "09_network_filters.csv")

# --- 04 metrics vs ERA5, cells ---
met <- bind_site_tables(per_site("04_metrics_%s.csv"))
pairs_main <- c("1_pblhmod_default_vs_era5_center", "2_pblhmod_default_vs_era5_bilinear", "6_pblhmod_clim_vs_era5_center")
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
# Exceptions to the reference-year rule (spec 09 amendment 2): used only when the standard rule finds no year
rex <- read_reference_year_exceptions(if (is.null(cfg$reference_year_exceptions_file)) NULL else
  resolve_path(cfg$reference_year_exceptions_file, must_exist = FALSE))
if (nrow(yc)) {
  ry <- dplyr::bind_rows(lapply(split(yc, yc$site_id), function(x) {
    e <- rex[rex$site_id == x$site_id[1], ]
    p <- if (nrow(e)) propose_reference_year(x, rule = e$rule, rule_reason = e$reason) else propose_reference_year(x)
    dplyr::bind_cols(dplyr::tibble(site_id = x$site_id[1]), p)
  }))
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
