# Spec 18 collection: network tables and figures from the per-site tables of scripts/18a_pblh_ec_site.R.
#   Reads <DIR>/18_<part>_<site>.csv (DIR = --from-tables, default <output_dir>/tables).
#   Writes <output_dir>/tables/18_variables.csv, 18_sl_share.csv (site + network rows), 18_sl_share_flags.csv,
#   18_ustar_cross.csv, 18_nee_residual.csv, 18_storage_h.csv, 18_storage_nights.csv, 18_panofsky.csv,
#   18_panofsky_rank.csv, 18_panofsky_bins.csv, 18_scf_h.csv; figures 18_sl_share, 18_nee_residual, 18_storage_h,
#   18_panofsky, 18s_scf_h.
#
# Usage: Rscript scripts/18b_pblh_ec_collect.R [--from-tables=DIR] [--label=TEXT] [--panofsky-sites=A,B,...]
#   --label: appended to every caption (e.g. "test: 4 local sites").
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), c("--from-tables", "--label", "--panofsky-sites"))
if (length(bad)) cli::cli_abort("Unknown argument{?s}: {.val {bad}}.")
cfg <- get_cfg()
out_dir <- resolve_path(cfg$output_dir, must_exist = FALSE)
tables_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
from <- opt("from-tables", tables_dir)
label <- opt("label")
tag <- if (is.null(label)) "" else paste0(" [", label, "]")
pan_sites <- trimws(strsplit(opt("panofsky-sites", "FR-Gri,CH-Dav,FI-Hyy,BE-Bra,DE-Tha,IT-Ren"), ",")[[1]])
hd <- pec_table_headers()

tb <- lapply(stats::setNames(pec_parts(), pec_parts()), function(p) pec_bind_part(from, p))
if (!nrow(tb$variables)) cli::cli_abort("No 18_variables_<site>.csv in {.path {from}}.")
cli::cli_inform("{nrow(tb$variables)} site{?s}: {paste(tb$variables$site_id, collapse = ', ')}.")

sl <- tb$sl_share
sl <- dplyr::bind_rows(sl, pec_network_median(sl, c("source", "row_basis", "daynight", "zl", "season"),
                                              c("share_gt_0.1", "share_gt_0.05", "share_gt_0.2", "zm_h_median", "zm_h_p90")))
flags <- pec_stable_night_flags(sl)
uc <- tb$ustar_cross
uc <- dplyr::bind_rows(uc, pec_network_median(uc, c("source", "filter_def"),
                                              c("share_above_removed", "share_pass_above", "share_above", "share_filtered")))
pan <- tb$panofsky
rk <- pec_panofsky_rank(pan)
pb <- tb$panofsky_bins
scf <- tb$scf_h
if (nrow(scf)) scf <- dplyr::bind_rows(scf, pec_network_median(scf, c("source", "x_var", "bin", "scf_var"), c("median", "x_mid")))
vars <- dplyr::left_join(tb$variables, dplyr::rename(flags, stable_night_gt25_sources = sources,
                                                     stable_night_max_share = max_share), by = "site_id")

hdr <- function(h) c(h, paste0("Collected by scripts/18b_pblh_ec_collect.R", tag, "; site_id = network: median over sites (n_sites with a value)."))
write_csv_header(vars, file.path(tables_dir, "18_variables.csv"), hdr(c(hd$variables, "stable_night_gt25_sources: sources with a stable-night share zm/h > 0.1 above 25 % (common rows).")))
write_csv_header(sl, file.path(tables_dir, "18_sl_share.csv"), hdr(hd$sl_share))
write_csv_header(flags, file.path(tables_dir, "18_sl_share_flags.csv"),
                 hdr("Spec 18 section 1: sites whose stable-night share of half-hours with zm/h > 0.1 exceeds 25 % for any source (common rows); sources ;-separated, max_share [-]."))
write_csv_header(uc, file.path(tables_dir, "18_ustar_cross.csv"), hdr(hd$ustar_cross))
write_csv_header(tb$nee_residual, file.path(tables_dir, "18_nee_residual.csv"), hdr(hd$nee_residual))
write_csv_header(tb$storage_h, file.path(tables_dir, "18_storage_h.csv"), hdr(hd$storage_h))
write_csv_header(tb$storage_nights, file.path(tables_dir, "18_storage_nights.csv"), hdr(hd$storage_nights))
write_csv_header(pan, file.path(tables_dir, "18_panofsky.csv"), hdr(hd$panofsky))
write_csv_header(rk, file.path(tables_dir, "18_panofsky_rank.csv"),
                 hdr("Spec 18 section 5: per row basis and source, sites ranked (n_sites), sites where the source ranks first (n_first), share_first [-]."))
write_csv_header(pb, file.path(tables_dir, "18_panofsky_bins.csv"), hdr(hd$panofsky_bins))
write_csv_header(scf, file.path(tables_dir, "18_scf_h.csv"), hdr(hd$scf_h))

cli::cli_h2("Section 1: sites with stable-night share > 25 % (zm/h > 0.1, common rows)")
print(as.data.frame(flags), row.names = FALSE)
cli::cli_h2("Section 5: ranks")
print(as.data.frame(rk), row.names = FALSE)

# figures
rows_txt <- "Rows: every half-hour from the labelling date, exclusions removed; zm = EC height - 2/3 hc; night = SW_IN_POT = 0."
cap1 <- fp_wrap(paste0("Stable nights (zm/L > 0.1), rows with h of all three model sources. h sources: ERA5 (tower cell, interpolated to the half-hour mid), PBLHmod, PBLHmod-clim. ",
                       pec_circ_note(), " Dashed line: 25 %. ", rows_txt, tag), 130)
save_figure(plot_pec_sl_share(sl, vars, cap1), "18_sl_share", 17, max(8, 3 + 0.45 * nrow(vars)), dir = fig_dir)
ord <- flags$site_id
x1 <- sl[sl$row_basis == "common" & sl$daynight == "night" & sl$zl == "stable" & sl$season == "all" &
           sl$source == "era5_blh_center" & sl$site_id != "network", ]
ord <- x1$site_id[order(-x1$share_gt_0.1)]
if (nrow(tb$nee_residual) && any(is.finite(tb$nee_residual$diff))) {
  okc <- unique(tb$nee_residual$site_id[is.finite(tb$nee_residual$diff) & tb$nee_residual$flux == "NEE" & tb$nee_residual$measure == "abs"])
  no_cell <- setdiff(tb$variables$site_id, okc)
  no_cell <- if (length(no_cell)) paste(no_cell, collapse = ", ") else "none"
  cap3 <- fp_wrap(paste0("Night, u*-passed (USTAR >= ONEFlux VUT threshold), measured (NEE_VUT_REF_QC = 0) half-hours: median of NEE_VUT_REF - RECO_NT_VUT_REF above the surface layer (zm/h > 0.1) minus within it, per u* / u*_thr bin (pooled over TA_F terciles; 'all' = all cells), ",
                         "cells with >= 30 rows each side; 90 % night-block bootstrap interval. h source: ERA5 (main) and PBLHmod. RECO_NT is fitted on the same night data, so a bias is diluted, not created. ",
                         pec_circ_note(), " Sites ordered by the ERA5 stable-night share above the surface layer. Sites without a valid cell (too few u*-passed rows above the surface layer): ", no_cell, ".", tag), 130)
  n_s <- length(unique(tb$nee_residual$site_id[is.finite(tb$nee_residual$diff)]))
  save_figure(plot_pec_nee_residual(tb$nee_residual, ord, caption = cap3), "18_nee_residual", 19,
              max(10, 5 + 3.2 * ceiling(n_s / 4)), dir = fig_dir)
} else cli::cli_alert_warning("18_nee_residual: no valid cell, not drawn.")
if (nrow(tb$storage_nights)) {
  cap4 <- fp_wrap(paste0("Sites with a CO2 profile (SC_FLAG 0): sum of positive SC per night (sunset to sunrise, >= 80 % of rows valid) against the night median h; lines: median per h bin. h source: ERA5 and PBLHmod. ",
                         pec_circ_note(), tag), 130)
  n_s <- length(unique(tb$storage_nights$site_id))
  save_figure(plot_pec_storage_h(tb$storage_nights, tb$storage_h, cap4), "18_storage_h", 19,
              max(10, 5 + 3.5 * ceiling(n_s / 4)), dir = fig_dir)
}
if (nrow(pb)) {
  cap5 <- fp_wrap(paste0("Day, 0.05 <= -zm/L <= 2, zm/h < 0.1 for the source: binned median σv/u* against -h/L (bins with >= 20 rows); dashed: Panofsky et al. (1977), (12 + 0.5 (-h/L))^(1/3). h sources: ERA5, PBLHmod, PBLHmod-clim and the measured MLH (day). ",
                         "A site-specific offset (terrain, meandering) is expected.", tag), 130)
  save_figure(plot_pec_panofsky(pb, pan_sites, cap5), "18_panofsky", 19, 13, dir = fig_dir)
}
if (nrow(scf)) {
  cap6 <- fp_wrap(paste0("Descriptive. Total spectral correction factors (EddyPro *_SCF_STAT) of unstable half-hours (zm/L < -0.1), binned medians per site (grey) and the median over sites (line), against zm/h and -h/L with h = ERA5. ",
                         "The low-frequency part of these factors is analytic and contains no h: any relation reflects their shared dependence on wind speed, zm and L, not a measured loss.", tag), 130)
  save_figure(plot_pec_scf(scf, "era5_blh_center", cap6), "18s_scf_h", 17, 16, dir = fig_dir)
}
cli::cli_alert_success("Spec 18 tables in {.path {tables_dir}}, figures in {.path {fig_dir}}.")
