# Spec 09 (footprint run 2026-11-05): rebuild the target-area masks from all layers of the KML files, before any
# footprint chunk runs. The earlier builder read only the first KML layer: it failed at CZ-Lnz (first layer empty on
# the cluster's GDAL) and may have missed polygons elsewhere. For each site: the cached mask (if any) is compared with
# the new one (target_area_mask(): every layer, polygons only, union); with --apply the new mask is written
# atomically. The site's chunks in the chunk list are counted (done = _meta.rds present).
# With --apply, a site whose mask changed (or had no mask) also has its finished chunks invalidated: their _meta.rds
# files are deleted, so 08b recomputes them (the chunk hash includes only the number of mask cells, so a changed mask
# with the same cell count would otherwise be skipped). The per-half-hour target-area fraction needs the full
# footprint fields, which are not stored: a changed mask means a full rerun of the site's chunks.
# Output: <output_dir>/tables/08_target_area_masks.csv (_check.csv without --apply; per site: layers, polygons, cells
#         and area old / new, first layer and its polygon area (what the earlier builder read), same,
#         chunks, done, invalidated) and <derived_dir>/footprint/chunk_list_rerun.csv (rows of the chunk list without
#         _meta.rds after the invalidation: failed, never run or invalidated chunks; --apply only).
#
# Usage: Rscript scripts/08t_target_area_masks.R [SITE ...] [--apply] [--chunk-list=PATH]
#   Without sites: every site of config/reference_years.csv. Without --apply: report only (the table), nothing written
#   or deleted. Run once without, check the table, then with --apply (it compares against the cached masks).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
apply <- "--apply" %in% args
cfg <- get_cfg()
site_ids <- grep("^--", args, value = TRUE, invert = TRUE)
ry <- read_reference_years(resolve_path(cfg$reference_years_file, must_exist = FALSE),
                           if (length(site_ids)) site_ids else NULL)
site_ids <- ry$site_id
sites <- read_sites(resolve_path(cfg$sites_file))
derived <- resolve_path(cfg$derived_dir, must_exist = FALSE)
fp_root <- file.path(derived, "footprint")
ta_dir <- resolve_path(cfg$target_areas_dir)
extra_sites <- footprint_extra_sites(resolve_path(cfg$footprint_extra_sites_file))
cl_file <- opt("chunk-list", file.path(fp_root, "chunk_list.csv"))
cl <- if (file.exists(cl_file)) readr::read_csv(cl_file, col_types = readr::cols(.default = "c"), progress = FALSE) else NULL
if (is.null(cl)) cli::cli_alert_warning("No chunk list {.path {cl_file}}: masks only, no chunk counts.")
meta_path <- function(d) {
  dir <- ifelse(d$case == "oper" & !d$site_id %in% extra_sites, file.path(fp_root, d$site_id, "oper_internal", "chunks"),
                file.path(fp_root, d$site_id, "chunks"))
  file.path(dir, paste0(d$site_id, "_", d$case, "_", d$part, "_meta.rds"))
}

grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)
rows <- list()
for (s in site_ids) {
  meta <- sites[sites$site_id == s, ]
  if (nrow(meta) != 1 || is.na(meta$lat) || is.na(meta$lon)) cli::cli_abort("{s}: no coordinates in sites.csv.")
  kml <- file.path(ta_dir, paste0(s, "_TA.kml"))
  if (!file.exists(kml)) cli::cli_abort("{s}: {.path {kml}} not found.")
  mask_file <- file.path(fp_root, s, paste0(s, "_target_area_mask.rds"))
  old <- if (file.exists(mask_file)) readRDS(mask_file) else NULL
  new <- target_area_mask(kml, meta$lat, meta$lon, grid)
  same <- !is.null(old) && identical(as.vector(old), as.vector(new))
  # what the earlier builder read: the first layer reported by this GDAL build (polygons only; 0 if it has none)
  first <- sf::st_layers(kml)$name[1]
  g1 <- sf::st_zm(sf::st_geometry(sf::st_read(kml, layer = first, quiet = TRUE)))
  g1 <- g1[sf::st_geometry_type(g1) %in% c("POLYGON", "MULTIPOLYGON")]
  crs <- sprintf("+proj=aeqd +lat_0=%.8f +lon_0=%.8f +x_0=0 +y_0=0 +datum=WGS84 +units=m", meta$lat, meta$lon)
  area_first <- if (length(g1)) sum(as.numeric(sf::st_area(sf::st_transform(g1, crs)))) else 0
  if (apply) write_rds_atomic(new, mask_file)
  d <- if (is.null(cl)) cl else cl[cl$site_id == s, ]
  mp <- if (is.null(d)) character() else meta_path(d)
  done <- file.exists(mp)
  n_inv <- 0L
  if (!same && any(done)) {
    if (apply) {
      n_inv <- sum(file.remove(mp[done]))
      cli::cli_alert_warning("{s}: mask changed, {n_inv} finished chunk{?s} invalidated (meta deleted).")
    } else {
      cli::cli_alert_warning("{s}: mask changed, {sum(done)} finished chunk{?s} to invalidate: rerun with --apply.")
    }
  }
  rows[[s]] <- dplyr::tibble(site_id = s, n_layers = attr(new, "n_layers"),
                             layers = paste(attr(new, "layers"), collapse = ";"), n_polygons = attr(new, "n_polygons"),
                             old_mask = !is.null(old), old_n_layers_attr = !is.null(attr(old, "n_layers")),
                             cells_old = if (is.null(old)) NA_integer_ else sum(old), cells_new = sum(new),
                             area_old_m2 = if (is.null(old)) NA_real_ else round(attr(old, "area_m2")),
                             area_new_m2 = round(attr(new, "area_m2")), first_layer = first,
                             area_first_layer_m2 = round(area_first), same = same, n_chunks = length(mp),
                             n_done = sum(done), n_invalidated = n_inv)
  cli::cli_inform("{s}: layers {attr(new, 'n_layers')} ({paste(attr(new, 'layers'), collapse = ', ')}), {attr(new, 'n_polygons')} polygon{?s}, {sum(new)} cells (old {if (is.null(old)) 'none' else sum(old)}), {if (same) 'same' else 'CHANGED'}; chunks {sum(done)}/{length(mp)} done.")
}
tab <- dplyr::bind_rows(rows)
tdir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(tab, file.path(tdir, paste0("08_target_area_masks", if (apply) "" else "_check", ".csv")))
print(as.data.frame(tab[, c("site_id", "n_layers", "n_polygons", "cells_old", "cells_new", "same", "n_chunks", "n_done",
                            "n_invalidated")]), row.names = FALSE)
if (apply && !is.null(cl)) {
  todo <- cl[!file.exists(meta_path(cl)), ]
  out <- file.path(fp_root, "chunk_list_rerun.csv")
  readr::write_csv(todo, out)
  cli::cli_alert_success("{nrow(todo)} of {nrow(cl)} chunks without _meta.rds: {.path {out}}.")
}
