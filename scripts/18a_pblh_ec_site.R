# Spec 18 per site: where the PBLH source matters for EC fluxes (sections 1-6). Analysis only: no PBLHmod or FFP
# run, no existing output changed. Call directly, never through run_site.R (the R/ changes of specs 17-19 change
# every run_site hash).
#   Reads per site the L2 FLUXES and FLUXNET HH files, the ONEFlux BIF (VUT u* thresholds), the standard files
#   era5_blh_center, pblhmod_default, pblhmod_clim and the primary measured source (config/measured_sites.csv).
#   Writes <output_dir>/tables/18_<part>_<site>.csv (`#` header) for part in variables, sl_share, ustar_cross,
#   nee_residual, storage_h, storage_nights, panofsky, panofsky_bins, scf_h. A part whose variables are missing is
#   skipped (empty table; reason in 18_variables part_<k>). Restartable: a site whose nine tables exist is skipped
#   unless --force.
#
# Usage: Rscript scripts/18a_pblh_ec_site.R [--sites=A,B | --sites-file=config/run_sites.txt | --site-index=N]
#          [--out=DIR] [--n-boot=1000] [--force]
#   --site-index=N: the N-th site of --sites-file (SLURM array task).
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
known <- c("--sites", "--sites-file", "--site-index", "--out", "--n-boot", "--force")
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), known)
if (length(bad)) cli::cli_abort("Unknown argument{?s}: {.val {bad}}.")
force <- "--force" %in% args
n_boot <- as.integer(opt("n-boot", "1000"))

cfg <- get_cfg()
out_dir <- opt("out", file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ids <- if (!is.null(opt("sites"))) trimws(strsplit(opt("sites"), ",")[[1]]) else {
  x <- trimws(readLines(opt("sites-file", here::here("config", "run_sites.txt")), warn = FALSE))
  x[nzchar(x) & !startsWith(x, "#")]
}
if (!is.null(opt("site-index"))) {
  k <- as.integer(opt("site-index"))
  if (is.na(k) || k < 1 || k > length(ids)) cli::cli_abort("--site-index={k} outside 1..{length(ids)}.")
  ids <- ids[k]
}

sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
excl <- read_exclusions(resolve_path(cfg$exclusions_file))
ref <- tryCatch(read_reference_years(resolve_path(cfg$reference_years_file)), error = function(e) NULL)
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
hd <- pec_table_headers()
parts <- c("variables", "sl_share", "ustar_cross", "nee_residual", "storage_h", "storage_nights", "panofsky",
           "panofsky_bins", "scf_h")
path_of <- function(p, s) file.path(out_dir, paste0("18_", p, "_", s, ".csv"))

for (s in ids) {
  cli::cli_h1(s)
  if (!force && all(file.exists(path_of(parts, s)))) {
    cli::cli_alert_info("{s}: all spec-18 tables exist, skipped (--force to redo).")
    next
  }
  t0 <- Sys.time()
  x <- pblh_ec_site_rows(s, cfg, sites, sites_tv, excl, ref, ms)
  r <- x$rows
  v <- x$vars
  cli::cli_inform("{nrow(r)} rows from {format(v$labelling_date)} ({v$n_excluded} excluded, {v$n_before_labelling} before the labelling date); measured: {x$meas_src}; storage: {v$storage}; u* threshold: {v$ustar_thr_source}.")
  res <- list()
  run_part <- function(k, need, fun) {
    miss <- pec_missing(need)
    if (length(miss)) {
      v[[paste0("part_", k)]] <<- paste0("skipped: ", paste(miss, collapse = ", "))
      cli::cli_alert_warning("Part {k} skipped: {paste(miss, collapse = ', ')}.")
      return(NULL)
    }
    t1 <- Sys.time()
    out <- fun()
    v[[paste0("part_", k)]] <<- "ok"
    v[[paste0("sec_part_", k)]] <<- round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1)
    out
  }
  has <- function(k) pec_has(r, k)
  res$sl_share <- run_part("1", list(zm = has("zm_m"), h = has("h_era5_blh_center")),
                           function() sl_share_table(r, s, meas_src = if (x$meas_type %in% "turbulence") x$meas_src else NA_character_))
  res$ustar_cross <- run_part("2", list(ustar_threshold = has("ustar_thr"), NEE_VUT_REF_QC = has("nee_qc"),
                                        SW_IN_POT = has("sw_in_pot")),
                              function() ustar_cross_table(r, s))
  res$nee_residual <- run_part("3", list(ustar_threshold = has("ustar_thr"), NEE_VUT_REF = has("nee_vut"),
                                         RECO_NT_VUT_REF = has("reco_nt"), TA_F = has("ta_f_k")),
                               function() nee_residual_table(r, s, fluxes = if (isTRUE(v$storage)) c("NEE", "FC") else "NEE",
                                                             n_boot = n_boot))
  nights <- run_part("4", list(storage_profile = isTRUE(v$storage)), function() storage_nights(r, s))
  res$storage_nights <- nights
  res$storage_h <- if (!is.null(nights)) storage_h_table(r, nights, s) else NULL
  pan <- run_part("5", list(V_SIGMA = has("v_sigma"), USTAR = has("ustar"), MO_LENGTH = has("L")),
                  function() panofsky_table(r, s, x$meas_src))
  res$panofsky <- pan$stats
  res$panofsky_bins <- pan$bins
  res$scf_h <- run_part("6", list(SCF = any(vapply(c("h_scf", "le_scf", "fc_scf", "tau_scf"), has, logical(1)))),
                        function() scf_h_table(r, s))
  v$sec_total <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  for (p in setdiff(parts, "variables")) {
    write_csv_header(if (is.null(res[[p]])) dplyr::tibble(site_id = character()) else res[[p]], path_of(p, s), hd[[p]])
  }
  write_csv_header(v, path_of("variables", s), hd$variables)
  cli::cli_alert_success("{s}: {v$sec_total} s; tables 18_*_{s}.csv in {.path {out_dir}}.")
}
