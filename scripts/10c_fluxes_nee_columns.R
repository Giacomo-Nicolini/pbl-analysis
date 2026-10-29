# NEE columns of the FLUXES and FLUXNET files of every station in config/sites.csv (G 2026-10-29, for spec 09:
# the footprint subsets need FLUXES NEE / NEE_DATA_FLAG and FLUXNET NEE_VUT_REF_QC). Header check, plus the number
# of non-missing values of NEE and NEE_DATA_FLAG (FLUXES; needs reading the two columns).
# Stations without an L2 folder are listed as such.
#
# Output (<output_dir>): tables/10c_fluxes_nee_columns.csv
#
# Usage: Rscript scripts/10c_fluxes_nee_columns.R [--headers-only]
#   Profile: env var R_CONFIG_ACTIVE (default | cluster); run on the cluster for all stations.

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

headers_only <- "--headers-only" %in% commandArgs(trailingOnly = TRUE)
cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
l2 <- resolve_path(cfg$icos_l2_dir)

res <- lapply(sites$site_id, function(s) {
  dir <- file.path(l2, paste0(s, "_output_L2I_2026"))
  site_re <- gsub("([.\\-])", "\\\\\\1", s)
  hit <- function(p) {
    h <- list.files(dir, pattern = p, full.names = TRUE)
    if (length(h) == 1) h else NULL
  }
  fx <- if (dir.exists(dir)) hit(paste0("^ICOSETC_", site_re, "_FLUXES_INTERIM_L2\\.zip$")) else NULL
  fn <- if (dir.exists(dir)) hit(paste0("^ICOSETC_", site_re, "_FLUXNET_FLUXMET_HH_.*\\.zip$")) else NULL
  hx <- if (!is.null(fx)) icos_zip_header(fx) else NULL
  hn <- if (!is.null(fn)) icos_zip_header(fn) else NULL
  out <- dplyr::bind_cols(dplyr::tibble(site_id = s, l2_folder = dir.exists(dir), fluxes_file = !is.null(fx),
                                        fluxnet_file = !is.null(fn)),
                          nee_column_check(hx, hn))
  out$n_rows <- NA_real_
  out$n_NEE <- NA_real_
  out$n_NEE_DATA_FLAG <- NA_real_
  out$n_NEE_DATA_FLAG_0 <- NA_real_
  if (!headers_only && !is.null(fx) && isTRUE(out$has_NEE)) {
    cols <- intersect(c("NEE", "NEE_DATA_FLAG"), hx)
    d <- readr::read_csv(fx, col_select = dplyr::all_of(cols), col_types = readr::cols(.default = "d"),
                         na = c("-9999", "-9999.0", "NA", ""), progress = FALSE, show_col_types = FALSE)
    out$n_rows <- nrow(d)
    out$n_NEE <- sum(!is.na(d$NEE))
    if ("NEE_DATA_FLAG" %in% cols) {
      out$n_NEE_DATA_FLAG <- sum(!is.na(d$NEE_DATA_FLAG))
      out$n_NEE_DATA_FLAG_0 <- sum(d$NEE_DATA_FLAG %in% 0)
    }
  }
  cli::cli_inform("{s}: folder {out$l2_folder}, NEE {out$has_NEE}, NEE_DATA_FLAG {out$has_NEE_DATA_FLAG}, NEE_VUT_REF_QC {out$has_NEE_VUT_REF_QC}")
  out
})
res <- dplyr::bind_rows(res)
readr::write_csv(res, file.path(tables_dir, "10c_fluxes_nee_columns.csv"))
print(as.data.frame(res[res$l2_folder, c("site_id", "has_NEE", "has_NEE_DATA_FLAG", "has_NEE_VUT_REF_QC", "n_rows",
                                         "n_NEE", "n_NEE_DATA_FLAG", "n_NEE_DATA_FLAG_0")]), row.names = FALSE)
cli::cli_inform("{sum(res$l2_folder)} of {nrow(res)} stations with an L2 folder; lacking NEE: {.val {res$site_id[res$fluxes_file & !res$has_NEE]}}; lacking NEE_DATA_FLAG: {.val {res$site_id[res$fluxes_file & !res$has_NEE_DATA_FLAG]}}; lacking NEE_VUT_REF_QC: {.val {res$site_id[res$fluxnet_file & !res$has_NEE_VUT_REF_QC]}}.")
