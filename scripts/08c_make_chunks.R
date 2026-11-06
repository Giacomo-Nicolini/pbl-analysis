# Spec 09 A5: footprint chunk list for the cluster array, from the approved reference years only.
# Reads config/reference_years.csv (site_id, year, reason, approved; committed after G's review of the proposal
# output/tables/08_reference_years.csv written by scripts/09_collect.R) and each site's full-period input table
# (scripts/08a_footprint_flags.R). One chunk per site x case x month of the year; a month whose estimated duration
# (half-hours x sec_per_ffp_halfhour / CPUs) exceeds chunk_target_min is split into halves (config.yml).
# Cases: every case of the input table; oper only at the extra sites (supplement) or with --with-oper.
# Output: <derived_dir>/footprint/chunk_list.csv (one row per array task) and
#         <output_dir>/tables/09_chunks_summary.csv (per site: chunks, half-hours, CPU-hours, longest chunk).
# With --pending (spec 13): also <derived_dir>/footprint/chunk_list_pending.csv, the rows of the list whose chunk is
#   new (no _meta.rds) or changed (stored hash differs from footprint_task_hash(), e.g. the common_pme groups of
#   pblhmod / pblhmod_clim / era5 once a measured case is added), with a column `status` (new / changed); submit it
#   with scripts/slurm_submit.R --what=chunks --chunk-list=PATH. Needs the FFP code (ffp_code_dir) for the hash and
#   reads the cached target-area masks. Run it without sites, so that chunk_list.csv stays the full list.
#
# Usage: Rscript scripts/08c_make_chunks.R [SITE ...] [--cpus=2] [--with-oper] [--pending]
#   Without sites: every site of config/reference_years.csv. A requested site without an approved year aborts.

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
cfg <- get_cfg()
cpus <- as.integer(opt("cpus", "2"))
with_oper <- "--with-oper" %in% args
pending <- "--pending" %in% args
site_ids <- grep("^--", args, value = TRUE, invert = TRUE)
ry <- read_reference_years(resolve_path(cfg$reference_years_file, must_exist = FALSE),
                           if (length(site_ids)) site_ids else NULL)
site_ids <- ry$site_id
extra_sites <- footprint_extra_sites(resolve_path(cfg$footprint_extra_sites_file))
derived <- resolve_path(cfg$derived_dir, must_exist = FALSE)
if (pending) {
  sites <- read_sites(resolve_path(cfg$sites_file))
  ffp_md5 <- unname(tools::md5sum(file.path(resolve_path(cfg$ffp_code_dir),
                                            c("calc_footprint_FFP.R", "calc_footprint_FFP_climatology.R"))))
  grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)
  rs <- ffp_run_settings()
  ta_dir <- resolve_path(cfg$target_areas_dir)
}

lst <- list()
for (s in site_ids) {
  inp_file <- file.path(derived, "footprint", s, paste0(s, "_footprint_inputs.rds"))
  if (!file.exists(inp_file)) cli::cli_abort("{s}: {.path {inp_file}} not found (run step 08a).")
  inp <- readRDS(inp_file)
  cases <- sub("^h_", "", grep("^h_", names(inp$tab), value = TRUE))
  if (!(with_oper || s %in% extra_sites)) cases <- setdiff(cases, "oper")
  yr <- ry$year[ry$site_id == s]
  cl <- make_chunk_list(inp, yr, cases, sec_per_hh = cfg$sec_per_ffp_halfhour, cpus = cpus, target_min = cfg$chunk_target_min)
  if (nrow(cl) == 0) cli::cli_abort("{s}: no valid half-hour in {yr}.")
  if (pending) {
    # chunk status as 08b sees it: same task, same hash (R/footprint_chunks.R)
    m <- sites[sites$site_id == s, ]
    mask <- target_area_mask_cached(file.path(derived, "footprint", s, paste0(s, "_target_area_mask.rds")),
                                    file.path(ta_dir, paste0(s, "_TA.kml")), m$lat, m$lon, grid)
    cl$status <- vapply(seq_len(nrow(cl)), function(j) {
      tk <- footprint_task(inp, cl$case[j], cl$part[j], as.POSIXct(cl$start[j], tz = "UTC"),
                           as.POSIXct(cl$end[j], tz = "UTC"))
      d <- if (cl$case[j] == "oper" && !s %in% extra_sites) file.path("oper_internal", "chunks") else "chunks"
      mf <- file.path(derived, "footprint", s, d, paste0(s, "_", cl$case[j], "_", cl$part[j], "_meta.rds"))
      if (!file.exists(mf)) return("new")
      h <- footprint_task_hash(tk, ffp_md5, rs$r_levels, grid, mask, rs$block)
      if (identical(readRDS(mf)$hash, h)) "done" else "changed"
    }, "")
  }
  lst[[s]] <- cl
  cli::cli_inform("{s} {yr}: {nrow(cl)} chunks ({sum(grepl('[ab]$', cl$part))} half-months), {sum(cl$n_halfhours)} FFP half-hours, {round(sum(cl$cpu_hours))} CPU-hours; cases {paste(unique(cl$case), collapse = ', ')}{if (pending) paste0('; pending ', sum(cl$status != 'done')) else ''}.")
}
cl_all <- dplyr::bind_rows(lst)
cl_all <- cl_all[order(-cl_all$est_min), ] # longest first: the array starts them early
out <- file.path(derived, "footprint", "chunk_list.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(cl_all[setdiff(names(cl_all), "status")], out)
summ <- dplyr::summarise(cl_all, .by = site_id, n_chunks = dplyr::n(), n_halfhours = sum(n_halfhours),
                         cpu_hours = round(sum(cpu_hours), 1), max_est_min = max(est_min))
tdir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(summ, file.path(tdir, "09_chunks_summary.csv"))
print(as.data.frame(summ), row.names = FALSE)
cli::cli_alert_success("{nrow(cl_all)} chunks, {round(sum(cl_all$cpu_hours))} CPU-hours ({cpus} CPUs per chunk, {cfg$sec_per_ffp_halfhour} s per half-hour): {.path {out}}")
if (pending) {
  pd <- cl_all[cl_all$status %in% c("new", "changed"), ]
  out_p <- file.path(derived, "footprint", "chunk_list_pending.csv")
  readr::write_csv(pd, out_p)
  if (nrow(pd) > 0) print(as.data.frame(dplyr::count(pd, site_id, case, status)), row.names = FALSE)
  cli::cli_alert_success("{nrow(pd)} pending chunks ({sum(pd$status == 'new')} new, {sum(pd$status == 'changed')} changed), {round(sum(pd$cpu_hours))} CPU-hours: {.path {out_p}}")
}
