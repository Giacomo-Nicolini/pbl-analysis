# Spec 19: why the PBLH source changes the footprint at CH-Dav and IT-Ren and not at other mountain sites.
# Analysis only: no PBLHmod or FFP run. Reads, per site, the step-08a footprint input table and the step-08b chunk
# metrics of the reference year (spec-15 rows `common_pe`, set etc, subsets valid_qc_used and valid), W_SIGMA of the
# L2 FLUXES file and the ONEFlux VUT u* thresholds (BIF file); the ERA5 cell table of 09_collect.R; the existing
# measured comparison tables (section 3).
#
# Outputs (<output_dir>/tables, `#` header; <output_dir>/figures):
#   19_site_descriptors.csv     section 1, long: site x variable x subset x stratum (terrain, tower, PBLH, turbulence,
#                               night-time flows, footprint, spec-15 status)
#   19_row_level.csv            section 2: R2 of the TA difference and of ln(x_80 ratio) at night (selected sites)
#   19_row_level_rows.csv    night rows of the selected sites (figure 19s_row_level only)
#   19_night_flows.csv          wind-direction histograms per site x day/night (figure 19s_night_flows)
#   19_night_closeness.csv      section 3, from 07_metrics_FI-Hyy.csv and 10_metrics_<site>.csv
#   19_spec15_check.csv         n and p10 of the TA difference here vs the spec-15 table (row-selection check)
#   figures 19s_alpine_drivers, 19s_row_level, 19s_night_flows
# Sites whose chunks are absent (local test) use the row basis flags_only (valid by the 08a flags, footprint NA).
#
# Usage: Rscript scripts/19a_site_descriptors.R [--sites-file=config/run_sites.txt | --sites=A,B]
#          [--row-sites=CH-Dav,IT-Ren,...] [--flow-sites=...] [--label-sites=CZ-BK1,IT-MBo,FR-Lqu]
#          [--cells=FILE] [--spec15=FILE] [--measured-tables=DIR]
#        Rscript scripts/19a_site_descriptors.R --from-tables=DIR   (figures only, from the collected 19_* tables)

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
known <- c("--sites-file", "--sites", "--row-sites", "--flow-sites", "--label-sites", "--cells", "--spec15",
           "--measured-tables", "--from-tables")
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), known)
if (length(bad)) cli::cli_abort("Unknown argument{?s}: {.val {bad}}.")
split <- function(x) trimws(strsplit(x, ",")[[1]])

cfg <- get_cfg()
out_dir <- resolve_path(cfg$output_dir, must_exist = FALSE)
tables_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
threshold_pp <- 5
min_n <- 50
row_sites <- split(opt("row-sites", "CH-Dav,IT-Ren,DE-Tha,CZ-BK1,FR-Fon,BE-Bra,IT-MBo"))
flow_sites <- split(opt("flow-sites", "CH-Dav,IT-Ren,DE-Tha,CZ-BK1,IT-MBo,BE-Bra"))
label_sites <- split(opt("label-sites", "CZ-BK1,IT-MBo,FR-Lqu"))
hd <- alpine_table_headers(threshold_pp, min_n)
rd <- function(p) readr::read_csv(p, comment = "#", show_col_types = FALSE, progress = FALSE, guess_max = 1e5)

draw_figures <- function(desc, nr, wh, fig_dir) {
  basis <- paste(sort(unique(desc$row_basis)), collapse = ", ")
  ft <- alpine_drivers_table(desc, label_sites = label_sites, min_n = min_n)
  if (any(is.finite(ft$x) & is.finite(ft$y))) {
    cap <- fp_wrap(paste0("Rows: spec-15 common_pe of the reference year (set etc, subset valid_qc_used; row basis: ",
                          basis, "). x: median ln(PBLHmod / ERA5) on the neutral and stable night rows; y: p10 of the",
                          " target-area difference ERA5 - PBLHmod in the panel's stratum. Open symbols: fewer than ",
                          min_n, " rows (not flagged as exception, spec 15). Dashed line: -", threshold_pp, " pp."))
    save_figure(plot_alpine_drivers(ft, threshold_pp, cap), "19s_alpine_drivers", 18, 12, dir = fig_dir)
  } else cli::cli_alert_warning("19s_alpine_drivers: no footprint values (row basis {basis}), not drawn.")
  if (nrow(nr) && any(is.finite(nr$ta_diff_pp))) {
    cap <- fp_wrap("Night rows (zm/L neutral and stable) of the reference year, spec-15 common_pe (etc, valid_qc_used).")
    save_figure(plot_alpine_row_level(nr, cap), "19s_row_level", 18, 14, dir = fig_dir)
  } else cli::cli_alert_warning("19s_row_level: no footprint values, not drawn.")
  wh <- wh[wh$site_id %in% flow_sites, ]
  if (nrow(wh)) {
    cap <- fp_wrap(paste0("Wind direction (L2 WD, direction the wind blows from), 30-degree bins; rows of the reference",
                          " year (row basis: ", basis, ", subset valid_qc_used, all regimes); day = SW_IN_POT > 0."))
    save_figure(plot_alpine_flows(wh, cap), "19s_night_flows", 18, 13, dir = fig_dir)
  }
}

from <- opt("from-tables")
if (!is.null(from)) {
  desc <- rd(file.path(from, "19_site_descriptors.csv"))
  nr <- if (file.exists(p <- file.path(from, "19_row_level_rows.csv"))) rd(p) else dplyr::tibble()
  wh <- rd(file.path(from, "19_night_flows.csv"))
  draw_figures(desc, nr, wh, fig_dir)
  cli::cli_alert_success("Figures redrawn from {.path {from}} into {.path {fig_dir}}.")
  quit(save = "no")
}

ids <- if (!is.null(opt("sites"))) split(opt("sites")) else {
  x <- trimws(readLines(opt("sites-file", here::here("config", "run_sites.txt")), warn = FALSE))
  x[nzchar(x) & !startsWith(x, "#")]
}
sites <- read_sites(resolve_path(cfg$sites_file))
ref <- read_reference_years(resolve_path(cfg$reference_years_file))
fp_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint")
cells_path <- opt("cells", file.path(tables_dir, "09_network_era5_cells.csv"))
cells_all <- if (file.exists(cells_path)) rd(cells_path) else NULL
if (is.null(cells_all)) cli::cli_alert_info("{.path {cells_path}} not found: ERA5 cell table computed per site.")

desc_l <- list(); nr_l <- list(); wh_l <- list(); log_l <- list()
for (s in ids) {
  cli::cli_h2(s)
  res <- tryCatch({
    yr <- ref$year[ref$site_id == s]
    if (length(yr) != 1) cli::cli_abort("no reference year in {.path {cfg$reference_years_file}}.")
    meta <- sites[sites$site_id == s, ]
    r <- fp_site_rows(s, yr, fp_dir, sets = "etc", subsets = c("valid", "valid_qc_used"), tab_cols = alpine_tab_cols())
    miss <- attr(r, "missing")
    months <- unique(format(r$timestamp_start_utc, "%Y-%m", tz = "UTC"))
    no_chunks <- vapply(c("pblhmod", "pblhmod_clim", "era5"), function(cs) {
      sum(startsWith(miss, paste0(cs, " "))) == length(months)
    }, logical(1))
    require_ffp <- !any(no_chunks)
    basis_lab <- if (require_ffp) "common_pe" else "flags_only"
    if (!require_ffp) {
      cli::cli_alert_warning("no chunks for {paste(names(no_chunks)[no_chunks], collapse = ', ')} in {yr}: row basis flags_only (footprint values NA).")
      # partial chunks (e.g. a one-month PC test) must not enter the flags_only descriptors
      for (cs in c("pblhmod", "pblhmod_clim", "era5")) r[[paste0("ta__", cs)]] <- r[[paste0("x80__", cs)]] <- NA_real_
    }
    # sigma_w from L2 FLUXES
    fx <- read_icos_fluxes(find_icos_files(s, cfg)$fluxes, meta$utc_offset_h)
    w <- fx$W_SIGMA[match(as.numeric(r$timestamp_start_utc), as.numeric(fx$timestamp_start_utc))]
    # ONEFlux VUT u* threshold of the year (local standard time), BIF file
    yl <- as.integer(format(r$timestamp_start_utc + 3600 * meta$utc_offset_h, "%Y", tz = "UTC"))
    thr_src <- "VUT"
    thr <- tryCatch({
      hh <- find_icos_files(s, cfg)$fluxnet
      th <- oneflux_ustar_ref(read_oneflux_bif(find_oneflux_bif(s, cfg, hh_path = hh), "GRP_FLUX_REF", hh_path = hh), "VUT")
      v <- th$ustar_threshold[match(yl, th$year)]
      if (all(is.na(v))) thr_src <- paste0("VUT: no threshold for ", paste(unique(yl), collapse = "/"))
      v
    }, error = function(e) {
      thr_src <<- paste0("none (", conditionMessage(e), ")")
      rep(NA_real_, nrow(r))
    })
    r <- alpine_row_vars(r, w, thr)
    bases <- list(valid_qc_used = alpine_basis(r, "etc", "valid_qc_used", require_ffp),
                  valid = alpine_basis(r, "etc", "valid", require_ffp))
    cells <- if (!is.null(cells_all)) cells_all[cells_all$site_id == s, ] else
      tryCatch(era5_cell_metadata(s, cfg, sites), error = function(e) NULL)
    d <- alpine_site_descriptors(r, bases, basis_lab, list(cells = cells, igbp = meta$igbp, ustar_thr_source = thr_src),
                                 threshold_pp, min_n)
    cli::cli_inform("{yr}: {sum(bases$valid_qc_used)} rows ({basis_lab}, valid_qc_used), {sum(bases$valid)} (valid); W_SIGMA on {sum(is.finite(w[bases$valid_qc_used]))}; u* threshold: {thr_src}.")
    list(desc = d, nr = if (s %in% row_sites) alpine_night_rows(r, bases$valid_qc_used) else NULL,
         wh = alpine_wd_hist(r, bases$valid_qc_used),
         log = dplyr::tibble(site_id = s, year = yr, row_basis = basis_lab, n_valid_qc_used = sum(bases$valid_qc_used),
                             n_valid = sum(bases$valid), status = "ok", message = NA_character_))
  }, error = function(e) {
    cli::cli_alert_danger("skipped: {conditionMessage(e)}")
    list(log = dplyr::tibble(site_id = s, year = NA_integer_, row_basis = NA_character_, n_valid_qc_used = NA_integer_,
                             n_valid = NA_integer_, status = "skipped", message = conditionMessage(e)))
  })
  desc_l[[s]] <- res$desc; nr_l[[s]] <- res$nr; wh_l[[s]] <- res$wh; log_l[[s]] <- res$log
}
desc <- dplyr::bind_rows(desc_l)
nr <- dplyr::bind_rows(nr_l)
wh <- dplyr::bind_rows(wh_l)
lg <- dplyr::bind_rows(log_l)
if (nrow(desc) == 0) cli::cli_abort("No site processed.")

write_csv_header(desc, file.path(tables_dir, "19_site_descriptors.csv"), hd$descriptors)
rl <- if (nrow(nr)) alpine_row_level(nr) else dplyr::tibble()
write_csv_header(rl, file.path(tables_dir, "19_row_level.csv"), hd$row_level)
if (nrow(nr)) write_csv_header(nr, file.path(tables_dir, "19_row_level_rows.csv"), hd$rows)
write_csv_header(wh, file.path(tables_dir, "19_night_flows.csv"), hd$flows)
cl <- alpine_night_closeness(opt("measured-tables", tables_dir))
write_csv_header(cl, file.path(tables_dir, "19_night_closeness.csv"), hd$closeness)

# Row-selection check against the spec-15 table (same n and p10 per site x regime x day/night)
sp <- opt("spec15", file.path(tables_dir, "09_network_fp_target_area.csv"))
if (file.exists(sp) && any(desc$row_basis == "common_pe")) {
  t15 <- rd(sp)
  t15 <- t15[t15$scope == "site" & t15$set == "etc" & t15$subset == "valid_qc_used" & t15$pair == "era5_vs_pblhmod" &
               t15$site_id %in% desc$site_id[desc$row_basis == "common_pe"], ]
  one <- function(v) {
    x <- desc[desc$subset %in% "valid_qc_used" & desc$variable == v & desc$regime != "neutral_stable",
              c("site_id", "regime", "daynight", "value")]
    names(x)[4] <- v
    x
  }
  mine <- dplyr::inner_join(one("ta_diff_n"), one("ta_diff_p10"), by = c("site_id", "regime", "daynight"))
  ck <- dplyr::inner_join(dplyr::mutate(t15[c("site_id", "regime", "daynight", "n", "p10_pp")], regime = as.character(regime),
                                        daynight = as.character(daynight)), mine, by = c("site_id", "regime", "daynight"))
  ck$n_equal <- ck$n == ck$ta_diff_n
  ck$p10_diff_pp <- ck$ta_diff_p10 - ck$p10_pp
  write_csv_header(ck, file.path(tables_dir, "19_spec15_check.csv"),
                   c("Spec 19: rows here vs the spec-15 target-area table (era5_vs_pblhmod, etc, valid_qc_used, site scope).",
                     "n, p10_pp = spec 15; ta_diff_n, ta_diff_p10 = spec 19; n_equal, p10_diff_pp = check."))
  cli::cli_inform("Spec-15 check: {sum(ck$n_equal)} of {nrow(ck)} cells with equal n; max |p10 difference| {signif(max(abs(ck$p10_diff_pp), na.rm = TRUE), 3)} pp.")
} else cli::cli_alert_info("Spec-15 check skipped (no common_pe site or no {.path {sp}}).")

draw_figures(desc, nr, wh, fig_dir)
cli::cli_h2("Sites")
print(as.data.frame(lg), row.names = FALSE)
cli::cli_alert_success("Spec 19 tables in {.path {tables_dir}}, figures in {.path {fig_dir}}.")
