# Spec 16 (phase 5b): footprint uncertainty from the FFP inputs. Step 16b: FFP runs of one chunk (site x design x
# part) of <derived_dir>/footprint_mc/chunk_list.csv (16a), restartable by inputs hash as 08b.
#   Design A: per sampled row the base run and each input at -1 / +1 sigma (17 runs; L alone, u* fixed), plus hc and z0
#   at -1 / +1 sigma with the etc roughness-sublayer test off (validity rsl_off, spec 16 amendment 1; 21 runs).
#   Design B: per subsampled row the joint draws of 16a for each h option (L linked to u*: L' = L0 (u*'/u*0)^3 / eps_L).
#   Every run: perturbed inputs (mc_perturb()), validity under the etc set recomputed (mc_validity(); ETC z0 multiplier
#   of 08a kept); FFP only for valid runs, once per distinct input set, with the 08b wrapper and settings
#   (ffp_halfhour_run(): Kljun calc_footprint_FFP_climatology() for one half-hour, umean only, z0 = NaN, rslayer = 1,
#   domain +-1000 m, dx = 1 m, smoothed as the code; ffp_ci_paper() for x_peak, x_50 ... x_90).
# Output: <derived_dir>/footprint_mc/<site>/<site>_<part>_runs.csv.gz (one row per run: perturbed inputs, validity,
#   first failing condition, TA fraction, distances), _meta.rds (hash, written last: a chunk without it is
#   recomputed).
#
# Usage (cluster array): Rscript scripts/16b_mc_run.R --chunk-list=PATH --chunk-row=N [--workers=2]
# Usage (PC): Rscript scripts/16b_mc_run.R SITE [--design=A,B] [--workers=2] (every chunk of the site in the list)
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), c("--chunk-list", "--chunk-row", "--workers", "--design"))
if (length(bad)) cli::cli_abort("Unknown option{?s} {.val {bad}}.")
workers <- as.integer(opt("workers", "2"))
cfg <- get_cfg()
derived <- resolve_path(cfg$derived_dir, must_exist = FALSE)
mc_dir <- file.path(derived, "footprint_mc")
cl_file <- opt("chunk-list", file.path(mc_dir, "chunk_list.csv"))
if (!file.exists(cl_file)) cli::cli_abort("Chunk list {.path {cl_file}} not found: run scripts/16a_mc_sample.R first.")
cl_tab <- readr::read_csv(cl_file, col_types = readr::cols(first_row = "i", last_row = "i", .default = "c"), progress = FALSE)
if (!is.null(opt("chunk-row"))) {
  row <- as.integer(opt("chunk-row"))
  if (is.na(row) || row < 1 || row > nrow(cl_tab)) cli::cli_abort("--chunk-row must be in 1..{nrow(cl_tab)}.")
  todo <- cl_tab[row, ]
} else {
  site <- grep("^--", args, value = TRUE, invert = TRUE)
  if (length(site) != 1) cli::cli_abort("Pass one station code, or --chunk-list and --chunk-row.")
  des <- strsplit(opt("design", "A,B"), ",")[[1]]
  todo <- cl_tab[cl_tab$site_id == site & cl_tab$design %in% des, ]
  if (!nrow(todo)) cli::cli_abort("No chunk of {site} (design {paste(des, collapse = ', ')}) in {.path {cl_file}}.")
}
sites <- read_sites(resolve_path(cfg$sites_file))
ffp_dir <- resolve_path(cfg$ffp_code_dir)
ffp_md5 <- unname(tools::md5sum(file.path(ffp_dir, c("calc_footprint_FFP.R", "calc_footprint_FFP_climatology.R"))))
r_levels <- ffp_run_settings()$r_levels
grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)

for (j in seq_len(nrow(todo))) {
  ck <- todo[j, ]
  site_id <- ck$site_id
  draws <- readRDS(file.path(mc_dir, paste0(site_id, "_draws.rds")))
  meta <- sites[sites$site_id == site_id, ]
  fp_dir <- file.path(derived, "footprint", site_id)
  mask <- target_area_mask_cached(file.path(fp_dir, paste0(site_id, "_target_area_mask.rds")),
                                  file.path(resolve_path(cfg$target_areas_dir), paste0(site_id, "_TA.kml")),
                                  meta$lat, meta$lon, grid)
  smp <- draws$sample
  rows <- if (ck$design == "A") smp else smp[smp$in_b, ]
  rows <- rows[ck$first_row:ck$last_row, ]
  runs <- mc_run_inputs(rows, ck$design, draws = list(z = draws$z_b, h_options = draws$h_options), spec = draws$spec,
                        m_etc = draws$m_etc)
  hash <- mc_chunk_hash(runs, ffp_md5, r_levels, grid, mask)
  out_dir <- file.path(mc_dir, site_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  b <- file.path(out_dir, paste0(site_id, "_", ck$part))
  mf <- paste0(b, "_meta.rds")
  if (file.exists(mf) && identical(readRDS(mf)$hash, hash)) {
    cli::cli_inform("{site_id} {ck$part}: done (same inputs hash), skipped.")
    next
  }
  t0 <- Sys.time()
  res <- mc_run_ffp(runs, ffp_dir, grid, mask, r_levels, workers = workers)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  readr::write_csv(res, paste0(b, "_runs.csv.gz"))
  n_ffp <- attr(res, "n_ffp")
  saveRDS(list(hash = hash, n_runs = nrow(res), n_ffp = n_ffp, seconds = secs, workers = workers, seed = draws$seed),
          mf)
  cli::cli_inform(sprintf("%s %s %s: %d rows, %d runs, %d invalid under their set, %d FFP calls (distinct valid inputs), %.0f s, %d worker(s), %.2f s per FFP call per core",
                          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), site_id, ck$part, nrow(rows), nrow(res),
                          sum(!res$valid_run), n_ffp, secs, workers, if (n_ffp) secs * workers / n_ffp else NA))
}
