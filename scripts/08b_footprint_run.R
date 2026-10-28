# Spec 08 (phase 5a): footprint runs for one station and period (input table from scripts/08a_footprint_flags.R).
# Per PBLH case and month (one chunk each; restartable: a chunk with the same inputs hash is skipped; parallel over
# chunks with --workers, PSOCK on Windows, forked on Linux):
#   A  crosswind-integrated distances x_peak, x_50, x_70, x_80, x_90 per half-hour, from the paper's equations
#      (ffp_ci_paper(); verified against Kljun's code by scripts/08v_ffp_verify.R);
#   B  per half-hour 2D footprint on the fixed grid (domain +-1000 m, dx = 1 m) from Kljun's code
#      (calc_footprint_FFP_climatology() with one half-hour, umean only, rslayer = 1, unsmoothed), smoothed as the code
#      does (ffp_smooth()); fraction inside the ICOS target area, source areas and distances for r = 50/70/80/90 %;
#   C  monthly climatologies from the summed B fields (the code's climatology is their mean, then smoothed), per
#      validity set x subset (valid, valid_qc, valid_qc_used, valid_qc_h) x day/night (all = day + night), on three row sets:
#      own (rows valid for the case), common_pe (valid for pblhmod and era5), common_pme (and measured, FI-Hyy).
#   Monthly sums: the unsmoothed summed fields and counts per group (day and night separately), float32, so that
#      climatologies over several months (e.g. annual) are built by summing months and smoothing at the end
#      (ffp_clim_from_sums()). About 14 MB per group, up to 96 groups per chunk.
# FFP computed only for rows valid in at least one requested set (never for rows failing the code checks).
# Output: <derived_dir>/footprint/<site>/chunks/<site>_<case>_<YYYY-MM>_hh.csv.gz, _clim.rds, _sums.rds, _meta.rds
#         (meta written last: a chunk without it is recomputed); log in <derived_dir>/footprint/<site>/logs/.
#
# Usage: Rscript scripts/08b_footprint_run.R SITE [--start=YYYY-MM-DD --end=YYYY-MM-DD] [--case=a,b] [--set=a,b]
#          [--workers=4] [--with-oper] [--limit=N (testing: N half-hours per chunk)]

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
with_oper <- "--with-oper" %in% args
site_id <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code.")
workers <- as.integer(opt("workers", "4"))
limit <- as.integer(opt("limit", NA)) # testing only: at most this many half-hours per chunk
sets <- strsplit(opt("set", "kljun,kljun_norsl,etc,none"), ",")[[1]]
r_levels <- c(0.5, 0.7, 0.8, 0.9)

cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
meta <- sites[sites$site_id == site_id, ]
fp_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint", site_id)
ffp_dir <- resolve_path(cfg$ffp_code_dir)
t0 <- if (!is.null(opt("start"))) as.POSIXct(opt("start"), tz = "UTC") else NULL
t1 <- if (!is.null(opt("end"))) as.POSIXct(opt("end"), tz = "UTC") + 86400 else NULL
per <- if (is.null(t0)) "" else paste0("_", format(t0, "%Y-%m"), if (format(t0, "%Y-%m") != format(t1 - 1, "%Y-%m")) paste0("_", format(t1 - 1, "%Y-%m")) else "")
inp_file <- file.path(fp_dir, paste0(site_id, "_footprint_inputs", per, ".rds"))
if (!file.exists(inp_file)) cli::cli_abort("Input table {.path {inp_file}} not found: run scripts/08a_footprint_flags.R with the same period.")
inp <- readRDS(inp_file)
tab <- inp$tab
cases_all <- sub("^h_", "", grep("^h_", names(tab), value = TRUE))
if (!with_oper) cases_all <- setdiff(cases_all, "oper")
cases <- if (is.null(opt("case"))) cases_all else intersect(strsplit(opt("case"), ",")[[1]], cases_all)
chunk_dir <- if (with_oper) file.path(fp_dir, "oper_internal", "chunks") else file.path(fp_dir, "chunks")
log_dir <- file.path(fp_dir, "logs")
for (d0 in c(chunk_dir, log_dir)) dir.create(d0, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("08b_", site_id, per, "_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".log"))
logm <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(...))
  cat(msg, "\n", file = log_file, append = TRUE)
  message(msg)
}

# Target-area mask on the grid (computed once per site, cached)
grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)
mask_file <- file.path(fp_dir, paste0(site_id, "_target_area_mask.rds"))
kml <- file.path(resolve_path(cfg$target_areas_dir), paste0(site_id, "_TA.kml"))
if (!file.exists(mask_file)) saveRDS(target_area_mask(kml, meta$lat, meta$lon, grid), mask_file)
mask <- readRDS(mask_file)
logm(site_id, ": target area ", round(attr(mask, "area_m2")), " m2 (", sum(mask), " grid cells)")

# Row sets for the climatologies
sub_cols <- inp$subsets
rowset_mask <- function(rs, u, s) {
  get <- function(k) { v <- sub_cols[[paste(u, s, k, sep = "__")]]; if (is.null(v)) rep(FALSE, nrow(tab)) else v }
  switch(rs, common_pe = get("pblhmod") & get("era5"), common_pme = get("pblhmod") & get("era5") & get("measured"))
}
has_measured <- "h_measured" %in% names(tab)
ffp_md5 <- unname(tools::md5sum(file.path(ffp_dir, c("calc_footprint_FFP.R", "calc_footprint_FFP_climatology.R"))))

# Tasks: case x month
month <- format(tab$timestamp_start_utc, "%Y-%m")
tasks <- list()
for (k in cases) {
  run <- Reduce(`|`, lapply(sets, function(s) sub_cols[[paste("valid", s, k, sep = "__")]]))
  for (mo in sort(unique(month[run]))) {
    idx <- which(run & month == mo)
    if (!is.na(limit)) idx <- utils::head(idx, limit)
    rows <- tab[idx, c("timestamp_start_utc", "zm_m", "umean", "L", "sigmav", "ustar", "wind_dir", "day", paste0("h_", k))]
    names(rows)[ncol(rows)] <- "h"
    grp <- list()
    for (s in sets) for (u in c("valid", "valid_qc", "valid_qc_used", "valid_qc_h")) {
      own <- sub_cols[[paste(u, s, k, sep = "__")]][idx]
      for (dn in c("day", "night")) {
        dmask <- if (dn == "day") rows$day %in% TRUE else rows$day %in% FALSE
        grp[[paste(s, u, dn, "own", sep = "|")]] <- own & dmask
        if (k %in% c("pblhmod", "era5", "era5_bilinear")) grp[[paste(s, u, dn, "common_pe", sep = "|")]] <- rowset_mask("common_pe", u, s)[idx] & dmask
        if (has_measured && k %in% c("pblhmod", "era5", "measured")) grp[[paste(s, u, dn, "common_pme", sep = "|")]] <- rowset_mask("common_pme", u, s)[idx] & dmask
      }
    }
    grp <- grp[vapply(grp, any, logical(1))]
    hash <- digest::digest(list(rows, grp, ffp_md5, r_levels, grid$domain, grid$dx, sum(mask), "v3"))
    tasks[[length(tasks) + 1]] <- list(site_id = site_id, case = k, month = mo, rows = rows, grp = grp, hash = hash)
  }
}
base <- function(tk) file.path(chunk_dir, paste0(tk$site_id, "_", tk$case, "_", tk$month))
done <- vapply(tasks, function(tk) {
  mf <- paste0(base(tk), "_meta.rds")
  file.exists(mf) && identical(readRDS(mf)$hash, tk$hash)
}, logical(1))
logm(site_id, per, ": ", length(tasks), " chunks (", sum(done), " already done), ", sum(vapply(tasks, function(tk) nrow(tk$rows), 1L)),
     " half-hours; cases ", paste(cases, collapse = ", "), "; sets ", paste(sets, collapse = ", "), "; workers ", workers)
tasks <- tasks[!done]

run_chunk <- function(tk, ffp_dir, grid, mask, r_levels, chunk_dir, log_file) {
  t_start <- Sys.time()
  env <- load_ffp_code(ffp_dir)
  rows <- tk$rows
  a <- ffp_ci_paper(rows$zm_m, rows$h, rows$umean, rows$ustar, r_levels)
  n <- nrow(rows)
  bm <- vector("list", n)
  sums <- list()
  cnt <- list()
  for (i in seq_len(n)) {
    fld <- function() ffp_field_halfhour(env, grid, rows$zm_m[i], rows$h[i], rows$umean[i], rows$L[i], rows$sigmav[i], rows$ustar[i], rows$wind_dir[i])
    f <- fld()
    if (is.null(f)) {
      bm[[i]] <- dplyr::tibble(ffp_ok = FALSE)
      next
    }
    # Guard (2026-10-27): one non-reproducible NA f_domain in the first test run. A field with non-finite cells is
    # recomputed once; cells still non-finite are set to 0 and counted, so they never enter the sums.
    n_retry <- 0L
    if (any(!is.finite(f))) { f <- fld(); n_retry <- 1L }
    n_nonfinite <- sum(!is.finite(f))
    f[!is.finite(f)] <- 0
    met <- function() ffp_field_metrics(ffp_smooth(f), grid, mask, r_levels, keep_region = NULL)$metrics
    m <- met()
    if (!all(is.finite(c(m$f_domain, m$f_target_area)))) { m <- met(); n_retry <- n_retry + 1L }
    bm[[i]] <- dplyr::bind_cols(dplyr::tibble(ffp_ok = TRUE, n_nonfinite = n_nonfinite, n_retry = n_retry), m)
    for (g in names(tk$grp)) {
      if (!tk$grp[[g]][i]) next
      if (is.null(sums[[g]])) { sums[[g]] <- f; cnt[[g]] <- 1L } else { sums[[g]] <- sums[[g]] + f; cnt[[g]] <- cnt[[g]] + 1L }
    }
  }
  hh <- dplyr::bind_cols(dplyr::tibble(timestamp_start_utc = rows$timestamp_start_utc, case = tk$case, h = rows$h), a,
                         dplyr::bind_rows(bm))
  # Climatologies: day, night, and all = (day sum + night sum) / (n day + n night)
  keys <- unique(sub("\\|(day|night)\\|", "|", names(sums)))
  clim <- list()
  for (key in keys) {
    p <- strsplit(key, "\\|")[[1]]
    parts <- list(day = paste(p[1], p[2], "day", p[3], sep = "|"), night = paste(p[1], p[2], "night", p[3], sep = "|"))
    for (dn in c("day", "night", "all")) {
      gk <- if (dn == "all") unlist(parts) else parts[[dn]]
      gk <- gk[gk %in% names(sums)]
      if (length(gk) == 0) next
      ntot <- sum(unlist(cnt[gk]))
      fm <- Reduce(`+`, sums[gk]) / ntot
      res <- ffp_field_metrics(ffp_smooth(fm), grid, mask, r_levels, keep_region = 0.8)
      clim[[length(clim) + 1]] <- list(set = p[1], subset = p[2], daynight = dn, rowset = p[3], n = ntot,
                                       metrics = res$metrics, region80 = res$region)
    }
  }
  b <- file.path(chunk_dir, paste0(tk$site_id, "_", tk$case, "_", tk$month))
  readr::write_csv(hh, paste0(b, "_hh.csv.gz"))
  saveRDS(clim, paste0(b, "_clim.rds"))
  ffp_sums_write(sums, cnt, paste0(b, "_sums.rds"))
  saveRDS(list(hash = tk$hash, n = n, seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))), paste0(b, "_meta.rds"))
  msg <- sprintf("%s %s %s %s: %d half-hours, %d climatologies, %.0f s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tk$site_id,
                 tk$case, tk$month, n, length(clim), as.numeric(difftime(Sys.time(), t_start, units = "secs")))
  cat(msg, "\n", file = log_file, append = TRUE)
  msg
}

if (length(tasks) > 0) {
  if (workers > 1) {
    cl <- if (.Platform$OS.type == "windows") parallel::makePSOCKcluster(workers) else parallel::makeForkCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    rfiles <- list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)
    parallel::clusterExport(cl, "rfiles", envir = environment())
    invisible(parallel::clusterEvalQ(cl, for (f in rfiles) source(f)))
    # largest chunks first
    tasks <- tasks[order(-vapply(tasks, function(tk) nrow(tk$rows), 1L))]
    out <- parallel::parLapplyLB(cl, tasks, run_chunk, ffp_dir = ffp_dir, grid = grid, mask = mask, r_levels = r_levels,
                                 chunk_dir = chunk_dir, log_file = log_file)
  } else {
    out <- lapply(tasks, run_chunk, ffp_dir = ffp_dir, grid = grid, mask = mask, r_levels = r_levels,
                  chunk_dir = chunk_dir, log_file = log_file)
  }
  for (m in out) message(m)
}
logm(site_id, per, ": done.")
