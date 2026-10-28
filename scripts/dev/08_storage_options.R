# Spec 08 -> spec 09: storage options for the monthly summed footprint fields (_sums.rds of 08b).
# Decides the storage design of the cluster run (spec 09): which groups keep a summed grid and in which form, with a
# target of <= ~20 GB per site and year (the test month stores ~0.8 GB per site x case x month as float32 + gzip,
# i.e. ~2.7 TB for 40 sites x 12 months x 7 cases; G 2026-10-27).
# Runs on: the 08b test-month chunks (2022-09, chunk hash v3, rerun of 2026-10-27 with 2 workers) of FI-Hyy and
# BE-Bra, cases pblhmod and era5 (sizes and metrics); all cases of the month for the group counts.
#
# Options (each a lossless or lossy encoding of one group's summed field, computed at 1 m in every case):
#   full_gzip  float32, gzip (current 08b format)          full_xz  float32, xz
#   c5         5 m block sums (tower cell at a block centre); reconstructed piecewise constant on the 1 m grid
#   sp3, sp4   sparse: the cells holding 1 - 1e-3 / 1 - 1e-4 of the domain mass (delta-coded int32 index + float32)
#   c5_sp4     5 m, then sparse 1 - 1e-4
# The "summed grids for the etc set only" option is a group selection, evaluated in the projection.
# Each lossy option is checked on 16 groups (sets etc, none x subsets valid, valid_qc_used x day, night; row set
# common_pe): climatology as ffp_clim_from_sums() (sum / n, ffp_smooth(), ffp_field_metrics()) vs the full field.
# Tolerances (G 2026-10-27): |d area_80| <= 1 %, |d xmax_80| <= max(5 m, 1 %), |d f_target_area| <= 0.1 pp,
# Jaccard of the 80 % region vs full >= 0.98, |d Jaccard ERA5 vs PBLHmod| <= 0.01.
#
# Output: <output_dir>/tables/08_storage_sizes.csv    per site, case, group, option: count, non-zero share, MB
#         <output_dir>/tables/08_storage_metrics.csv  per evaluated group and option: climatology metrics, Jaccard
#         <output_dir>/tables/08_storage_summary.csv  per option: max deviations, pass/fail, GB per site-year (all
#                                                     sets / etc only)
# Usage: Rscript scripts/dev/08_storage_options.R [--sites=FI-Hyy,BE-Bra] [--month=2022-09]

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
sites <- strsplit(opt("sites", "FI-Hyy,BE-Bra"), ",")[[1]]
month <- opt("month", "2022-09")
eval_cases <- c("pblhmod", "era5")
cfg <- get_cfg()
fp_root <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint")
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)
nn <- length(grid$x)
tol <- list(area_80_rel = 0.01, xmax_80_m = 5, xmax_80_rel = 0.01, f_target_area = 0.001, jaccard_full = 0.98,
            d_jaccard_pe = 0.01)

# ---- Encodings ----
f32 <- function(v) writeBin(as.vector(v), raw(), size = 4)
gz <- function(r) length(memCompress(r, "gzip"))
xz <- function(r) length(memCompress(r, "xz"))
# Block sums over k x k cells; the grid is padded with (k - 1) / 2 zero cells per side so that the tower cell
# (centre of the 2001 x 2001 grid) is the centre of a block.
coarsen <- function(m, k = 5) {
  p <- (k - 1) / 2
  big <- matrix(0, nn + 2 * p, nn + 2 * p)
  big[(p + 1):(p + nn), (p + 1):(p + nn)] <- m
  nb <- nrow(big) / k
  colSums(aperm(array(big, c(k, nb, k, nb)), c(1, 3, 2, 4)), dims = 2)
}
upsample <- function(b, k = 5) {
  p <- (k - 1) / 2
  kronecker(b / k^2, matrix(1, k, k))[(p + 1):(p + nn), (p + 1):(p + nn)]
}
sparsify <- function(m, delta) {
  v <- c(m)
  o <- order(v, decreasing = TRUE)
  cs <- cumsum(v[o])
  keep <- sort(o[seq_len(which(cs >= (1 - delta) * cs[length(cs)])[1])])
  list(idx = keep, val = v[keep], dim = dim(m))
}
desparse <- function(s) { m <- matrix(0, s$dim[1], s$dim[2]); m[s$idx] <- s$val; m }
sp_bytes <- function(s) gz(writeBin(diff(c(0L, as.integer(s$idx))), raw(), size = 4)) + gz(f32(s$val))
options <- list(
  full_gzip = list(size = function(m) gz(f32(m)), recon = function(m) m),
  full_xz = list(size = function(m) xz(f32(m)), recon = NULL),
  c5 = list(size = function(m) gz(f32(coarsen(m))), recon = function(m) upsample(coarsen(m))),
  sp3 = list(size = function(m) sp_bytes(sparsify(m, 1e-3)), recon = function(m) desparse(sparsify(m, 1e-3))),
  sp4 = list(size = function(m) sp_bytes(sparsify(m, 1e-4)), recon = function(m) desparse(sparsify(m, 1e-4))),
  c5_sp4 = list(size = function(m) sp_bytes(sparsify(coarsen(m), 1e-4)),
                recon = function(m) upsample(desparse(sparsify(coarsen(m), 1e-4)))))
eval_keys <- paste(rep(c("etc", "none"), each = 4), rep(rep(c("valid", "valid_qc_used"), each = 2), 2),
                   rep(c("day", "night"), 4), "common_pe", sep = "|")

# ---- Sizes (all groups) and metrics (evaluated groups) ----
sizes <- list(); mets <- list(); regions <- list(); ngroups <- list()
for (site in sites) {
  chunk_dir <- file.path(fp_root, site, "chunks")
  mask <- readRDS(file.path(fp_root, site, paste0(site, "_target_area_mask.rds")))
  for (path in list.files(chunk_dir, pattern = paste0("_", month, "_sums\\.rds$"), full.names = TRUE)) {
    cs <- sub(paste0("^", site, "_(.*)_", month, "_sums\\.rds$"), "\\1", basename(path))
    x <- readRDS(path)
    ngroups[[length(ngroups) + 1]] <- dplyr::tibble(site = site, case = cs, file_mb = file.size(path) / 1e6,
      n_groups = length(x$data), n_groups_etc = sum(startsWith(names(x$data), "etc|")))
    cli::cli_inform("{site} {cs}: {length(x$data)} groups, {round(file.size(path) / 1e6)} MB")
    if (!cs %in% eval_cases) { rm(x); invisible(gc()); next }
    for (g in names(x$data)) {
      m <- matrix(readBin(x$data[[g]], "double", n = prod(x$dim), size = 4), x$dim[1], x$dim[2])
      p <- strsplit(g, "\\|")[[1]]
      sz <- vapply(options, function(o) as.numeric(o$size(m)), numeric(1))
      sizes[[length(sizes) + 1]] <- dplyr::tibble(site = site, case = cs, group = g, set = p[1], subset = p[2],
        daynight = p[3], rowset = p[4], n = x$n[[g]], share_nonzero = mean(m > 0), raw_mb = 4 * length(m) / 1e6,
        option = names(sz), mb = sz / 1e6)
      if (g %in% eval_keys) for (on in names(options)) {
        if (is.null(options[[on]]$recon)) next
        res <- ffp_field_metrics(ffp_smooth(options[[on]]$recon(m) / x$n[[g]]), grid, mask, keep_region = 0.8)
        regions[[paste(site, cs, g, on, sep = "#")]] <- res$region
        mets[[length(mets) + 1]] <- dplyr::bind_cols(dplyr::tibble(site = site, case = cs, group = g, option = on,
          share80_in_ta = region_share_inside(res$region, mask)), res$metrics)
      }
      rm(m)
    }
    rm(x); invisible(gc())
  }
}
sizes <- dplyr::bind_rows(sizes)
mets <- dplyr::bind_rows(mets)
ngroups <- dplyr::bind_rows(ngroups)
rg <- function(site, case, group, option) regions[[paste(site, case, group, option, sep = "#")]]
mets$jaccard_vs_full <- mapply(function(s, c0, g, o) jaccard_index(rg(s, c0, g, o), rg(s, c0, g, "full_gzip")),
                               mets$site, mets$case, mets$group, mets$option)
mets$jaccard_era5_pblhmod <- mapply(function(s, c0, g, o) if (c0 != "era5") NA_real_ else
  jaccard_index(rg(s, "era5", g, o), rg(s, "pblhmod", g, o)), mets$site, mets$case, mets$group, mets$option)
ref <- mets[mets$option == "full_gzip", c("site", "case", "group", "area_80", "xmax_80", "f_target_area", "f_domain",
                                          "jaccard_era5_pblhmod")]
mets <- dplyr::left_join(mets, ref, by = c("site", "case", "group"), suffix = c("", "_full"))
mets <- dplyr::mutate(mets,
  d_area_80_rel = area_80 / area_80_full - 1, d_xmax_80_m = xmax_80 - xmax_80_full,
  d_f_target_area_pp = 100 * (f_target_area - f_target_area_full), d_f_domain_pp = 100 * (f_domain - f_domain_full),
  d_jaccard_pe = jaccard_era5_pblhmod - jaccard_era5_pblhmod_full)

# ---- Summary: deviations vs tolerances and projected size per site and year ----
# Size per site-year = 12 x (sum over the cases of the month of n_groups x mean MB per group of the option), the mean
# per group taken over the evaluated cases of the same site (all sets, or etc groups only).
proj <- function(o, etc_only) {
  vapply(sites, function(s) {
    sz <- sizes[sizes$site == s & sizes$option == o & (!etc_only | sizes$set == "etc"), ]
    ng <- ngroups[ngroups$site == s, ]
    12 * sum(if (etc_only) ng$n_groups_etc else ng$n_groups) * mean(sz$mb) / 1e3
  }, numeric(1))
}
summ <- dplyr::bind_rows(lapply(names(options), function(o) {
  mo <- mets[mets$option == o, ]
  pa <- proj(o, FALSE); pe <- proj(o, TRUE)
  ev <- nrow(mo) > 0
  mx <- function(v) if (ev) max(abs(v), na.rm = TRUE) else NA_real_
  dplyr::tibble(option = o,
    max_abs_d_area_80_pct = 100 * mx(mo$d_area_80_rel), max_abs_d_xmax_80_m = mx(mo$d_xmax_80_m),
    max_abs_d_xmax_80_pct = 100 * mx(mo$d_xmax_80_m / mo$xmax_80_full),
    max_abs_d_f_target_area_pp = mx(mo$d_f_target_area_pp), max_abs_d_f_domain_pp = mx(mo$d_f_domain_pp),
    min_jaccard_vs_full = if (ev) min(mo$jaccard_vs_full) else NA_real_, max_abs_d_jaccard_pe = mx(mo$d_jaccard_pe),
    pass = if (!ev) NA else all(abs(mo$d_area_80_rel) <= tol$area_80_rel,
      abs(mo$d_xmax_80_m) <= pmax(tol$xmax_80_m, tol$xmax_80_rel * mo$xmax_80_full),
      abs(mo$d_f_target_area_pp) <= 100 * tol$f_target_area, mo$jaccard_vs_full >= tol$jaccard_full,
      abs(mo$d_jaccard_pe) <= tol$d_jaccard_pe, na.rm = TRUE),
    !!!stats::setNames(as.list(pa), paste0("gb_site_year_all_", sites)),
    !!!stats::setNames(as.list(pe), paste0("gb_site_year_etc_", sites)))
}))
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(sizes, file.path(tables_dir, "08_storage_sizes.csv"))
readr::write_csv(mets, file.path(tables_dir, "08_storage_metrics.csv"))
readr::write_csv(summ, file.path(tables_dir, "08_storage_summary.csv"))
cli::cli_h2("Storage options ({month}; sizes per site and year, from the test month x 12)")
print(as.data.frame(summ), row.names = FALSE, digits = 3)
cli::cli_h2("Groups per case")
print(as.data.frame(ngroups), row.names = FALSE, digits = 3)
