# Spec 08 (phase 5a) / spec 09 A5: footprint runs (input table from scripts/08a_footprint_flags.R).
# One chunk = site x PBLH case x month, or half-month (YYYY-MMa / YYYY-MMb) when the chunk list splits it
# (R/footprint_chunks.R). Restartable: a chunk with the same inputs hash is skipped. Within a chunk the half-hours are
# processed in one-day blocks by --workers processes (fork on Linux, PSOCK on Windows); results do not depend on the
# number of workers.
#   A  crosswind-integrated distances x_peak, x_50, x_70, x_80, x_90 per half-hour, from the paper's equations
#      (ffp_ci_paper(); verified against Kljun's code by scripts/08v_ffp_verify.R);
#   B  per half-hour 2D footprint on the fixed grid (domain +-1000 m, dx = 1 m) from Kljun's code
#      (calc_footprint_FFP_climatology() with one half-hour, umean only, rslayer = 1, unsmoothed), smoothed as the code
#      does (ffp_smooth()); fraction inside the ICOS target area, source areas and distances for r = 50/70/80/90 %;
#   C  summed fields per validity set x subset (valid, valid_qc, valid_qc_used, valid_qc_h) x day/night x row set
#      (own; common_pe: also valid for pblhmod and era5; common_pme: and measured), stored as 5 m x 5 m block sums
#      (storage c5, float32; spec 09 A5), and the climatologies (day, night, all) from them: sum, upsample, smooth,
#      normalise, as for the annual climatologies of scripts/08_footprint_compare.R.
# FFP computed only for rows valid in at least one requested set (never for rows failing the code checks).
# oper chunks go to <...>/oper_internal/chunks/ except at the extra sites (supplement).
# Output: <derived_dir>/footprint/<site>/chunks/<site>_<case>_<part>_hh.csv.gz, _clim.rds, _sums.rds, _meta.rds
#         (meta written last: a chunk without it is recomputed); log in <derived_dir>/footprint/<site>/logs/.
#
# Usage (site mode, PC tests): Rscript scripts/08b_footprint_run.R SITE [--start=YYYY-MM-DD --end=YYYY-MM-DD]
#          [--case=a,b] [--set=a,b] [--workers=4] [--with-oper] [--limit=N (testing: N half-hours per chunk)]
#   Uses the input table of the same period (08a with --start/--end), one chunk per case x month.
# Usage (chunk mode, cluster array): Rscript scripts/08b_footprint_run.R --chunk-list=PATH --chunk-row=N [--workers=6]
#   Runs row N of the chunk list of scripts/08c_make_chunks.R (input table of the full period).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
with_oper <- "--with-oper" %in% args
workers <- as.integer(opt("workers", "4"))
limit <- as.integer(opt("limit", NA)) # testing only: at most this many half-hours per chunk
sets <- strsplit(opt("set", "kljun,kljun_norsl,etc,none"), ",")[[1]]
r_levels <- c(0.5, 0.7, 0.8, 0.9)
block <- 5L # storage c5 (G 2026-10-28)
cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
extra_sites <- footprint_extra_sites(resolve_path(cfg$footprint_extra_sites_file))

if (!is.null(opt("chunk-list"))) {
  cl_tab <- readr::read_csv(opt("chunk-list"), col_types = readr::cols(.default = "c"), progress = FALSE)
  row <- as.integer(opt("chunk-row"))
  if (is.na(row) || row < 1 || row > nrow(cl_tab)) cli::cli_abort("--chunk-row must be in 1..{nrow(cl_tab)}.")
  ck <- cl_tab[row, ]
  site_id <- ck$site_id
  per <- ""
  todo <- dplyr::tibble(case = ck$case, part = ck$part, start = as.POSIXct(ck$start, tz = "UTC"),
                        end = as.POSIXct(ck$end, tz = "UTC"))
} else {
  site_id <- grep("^--", args, value = TRUE, invert = TRUE)
  if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code, or --chunk-list and --chunk-row.")
  t0 <- if (!is.null(opt("start"))) as.POSIXct(opt("start"), tz = "UTC") else NULL
  t1 <- if (!is.null(opt("end"))) as.POSIXct(opt("end"), tz = "UTC") + 86400 else NULL
  per <- if (is.null(t0)) "" else paste0("_", format(t0, "%Y-%m"), if (format(t0, "%Y-%m") != format(t1 - 1, "%Y-%m")) paste0("_", format(t1 - 1, "%Y-%m")) else "")
  todo <- NULL
}
meta <- sites[sites$site_id == site_id, ]
extra <- site_id %in% extra_sites
fp_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint", site_id)
ffp_dir <- resolve_path(cfg$ffp_code_dir)
inp_file <- file.path(fp_dir, paste0(site_id, "_footprint_inputs", per, ".rds"))
if (!file.exists(inp_file)) cli::cli_abort("Input table {.path {inp_file}} not found: run scripts/08a_footprint_flags.R with the same period.")
inp <- readRDS(inp_file)
tab <- inp$tab
cases_all <- sub("^h_", "", grep("^h_", names(tab), value = TRUE))
if (!with_oper && !extra) cases_all <- setdiff(cases_all, "oper")
log_dir <- file.path(fp_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("08b_", site_id, per, "_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid(), ".log"))
logm <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(...))
  cat(msg, "\n", file = log_file, append = TRUE)
  message(msg)
}
chunk_dir_of <- function(case) {
  d <- if (case == "oper" && !extra) file.path(fp_dir, "oper_internal", "chunks") else file.path(fp_dir, "chunks")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Target-area mask on the grid (computed once per site, cached)
grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)
mask_file <- file.path(fp_dir, paste0(site_id, "_target_area_mask.rds"))
kml <- file.path(resolve_path(cfg$target_areas_dir), paste0(site_id, "_TA.kml"))
if (!file.exists(mask_file)) saveRDS(target_area_mask(kml, meta$lat, meta$lon, grid), mask_file)
mask <- readRDS(mask_file)
logm(site_id, ": target area ", round(attr(mask, "area_m2")), " m2 (", sum(mask), " grid cells)")
ffp_md5 <- unname(tools::md5sum(file.path(ffp_dir, c("calc_footprint_FFP.R", "calc_footprint_FFP_climatology.R"))))

# Tasks: one chunk-list row, or case x month of the period
if (is.null(todo)) {
  cases <- if (is.null(opt("case"))) cases_all else intersect(strsplit(opt("case"), ",")[[1]], cases_all)
  months <- sort(unique(format(tab$timestamp_start_utc, "%Y-%m")))
  todo <- dplyr::bind_rows(lapply(cases, function(k) dplyr::bind_rows(lapply(months, function(mo) {
    dplyr::mutate(chunk_parts(mo), case = k, .before = 1)
  }))))
} else if (!todo$case %in% cases_all) {
  cli::cli_abort("Case {.val {todo$case}} not in the input table of {site_id} (or oper without --with-oper).")
}
tasks <- list()
for (j in seq_len(nrow(todo))) {
  tk <- footprint_task(inp, todo$case[j], todo$part[j], todo$start[j], todo$end[j], sets, limit)
  if (is.null(tk)) next
  tk$hash <- digest::digest(list(tk$rows, tk$grp, ffp_md5, r_levels, grid$domain, grid$dx, sum(mask), block, "v4"))
  tk$chunk_dir <- chunk_dir_of(tk$case)
  tasks[[length(tasks) + 1]] <- tk
}
done <- vapply(tasks, function(tk) {
  mf <- file.path(tk$chunk_dir, paste0(tk$site_id, "_", tk$case, "_", tk$part, "_meta.rds"))
  file.exists(mf) && identical(readRDS(mf)$hash, tk$hash)
}, logical(1))
logm(site_id, per, ": ", length(tasks), " chunks (", sum(done), " already done), ", sum(vapply(tasks, function(tk) nrow(tk$rows), 1L)),
     " half-hours; cases ", paste(unique(todo$case), collapse = ", "), "; sets ", paste(sets, collapse = ", "), "; workers ", workers)
for (tk in tasks[!done]) {
  logm(run_footprint_chunk(tk, ffp_dir, grid, mask, r_levels, tk$chunk_dir, workers = workers, block = block))
}
logm(site_id, per, ": done.")
