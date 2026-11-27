# Spec 17 supplement (G 2026-11-27): redraw 17s_gamma_ratio_vs_height from the table 17_gamma_ratio_vs_height.csv
# written by 17b (part D), without rerunning 17a / 17b. Same figure and caption as 17b.
#   Reads <DIR>/17_gamma_ratio_vs_height.csv; writes the figure to --figures=DIR (default <output_dir>/figures).
#
# Usage: Rscript scripts/17r_gamma_ratio_height_figure.R --tables=DIR [--figures=DIR]
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
bad <- setdiff(sub("=.*$", "", args), c("--tables", "--figures"))
if (length(bad)) cli::cli_abort("Unknown argument{?s} {.val {bad}}.")
from <- opt("tables")
if (is.null(from) || !dir.exists(from)) cli::cli_abort("Pass --tables=DIR with 17_gamma_ratio_vs_height.csv.")
fig_dir <- opt("figures")
if (is.null(fig_dir)) fig_dir <- file.path(resolve_path(get_cfg()$output_dir, must_exist = FALSE), "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

hv <- readr::read_csv(file.path(from, "17_gamma_ratio_vs_height.csv"), comment = "#", show_col_types = FALSE,
                      progress = FALSE)
st <- gamma_ratio_height_stats(hv)
print(st$rho)
print(st$groups)
f <- save_figure(plot_gamma_ratio_vs_height(hv, gamma_ratio_height_caption(hv)), "17s_gamma_ratio_vs_height", 17, 13,
                 dir = fig_dir)
cli::cli_alert_success("Written: {paste(f, collapse = ', ')}.")
