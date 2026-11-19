# Spec 15: re-render the network footprint figures (09_fp_*, 09s_fp_*) and the exceptions table from collected
# tables, without the footprint chunks (close-out 2026-11-19: minimum n of the exception cells).
#   Reads <DIR>/09_network_fp_{target_area,fetch,topaffected,icos_criteria,exceptions}.csv. Exceptions: the status
#   (exception if n >= --min-n, else below_min_n) is (re)computed from n, so tables written before the minimum work.
#   Site order: <DIR>/09_network_fp_site_order.csv (written by 09_collect.R) or --site-order=A,B,... (top to bottom).
# Output: <output_dir>/tables/09_network_fp_exceptions.csv (with status), figures 09_fp_* and 09s_fp_*.
#
# Usage: Rscript scripts/09r_fp_figures.R --tables=DIR [--site-order=SITE,SITE,...] [--min-n=50]
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
bad <- setdiff(sub("=.*$", "", args), c("--tables", "--site-order", "--min-n"))
if (length(bad)) cli::cli_abort("Unknown argument{?s} {.val {bad}}.")
from <- opt("tables")
if (is.null(from) || !dir.exists(from)) cli::cli_abort("Pass --tables=DIR with the 09_network_fp_*.csv tables.")
min_n <- as.numeric(opt("min-n", "50"))
cfg <- get_cfg()
tdir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)

rd <- function(k) readr::read_csv(file.path(from, paste0("09_network_fp_", k, ".csv")), comment = "#",
                                  show_col_types = FALSE, progress = FALSE)
fp_out <- lapply(stats::setNames(nm = c("target_area", "fetch", "topaffected", "icos_criteria", "exceptions")), rd)
ex <- fp_out$exceptions[setdiff(names(fp_out$exceptions), c("min_n", "status"))]
fp_out$exceptions <- fp_exception_status(ex, min_n)
write_csv_header(fp_out$exceptions, file.path(tdir, "09_network_fp_exceptions.csv"), fp_table_headers(5, min_n)$exceptions)

so_file <- file.path(from, "09_network_fp_site_order.csv")
site_order <- if (!is.null(opt("site-order"))) {
  strsplit(opt("site-order"), ",")[[1]]
} else if (file.exists(so_file)) {
  so <- readr::read_csv(so_file, comment = "#", show_col_types = FALSE, progress = FALSE)
  so$site_id[order(so$order)]
} else cli::cli_abort("No site order: {.path {so_file}} missing and no --site-order.")
miss <- setdiff(unique(fp_out$icos_criteria$site_id), site_order)
if (length(miss)) cli::cli_abort("Sites missing from the site order: {paste(miss, collapse = ', ')}.")

exs <- fp_exception_sites(fp_out$exceptions)
lo <- fp_out$exceptions[fp_out$exceptions$status == "below_min_n", ]
cli::cli_inform("Exception sites (n >= {min_n}): {if (length(exs)) paste(exs, collapse = ', ') else 'none'}; cells below the minimum: {if (nrow(lo)) paste0(lo$site_id, ' ', lo$pair, ' ', lo$regime, '_', lo$daynight, ' n ', lo$n, collapse = '; ') else 'none'}.")
f <- fp_network_figures(fp_out, site_order)
cli::cli_alert_success("09_network_fp_exceptions.csv and {length(f)} figure{?s} written: {paste(f, collapse = ', ')}.")
