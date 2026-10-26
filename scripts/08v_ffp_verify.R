# Spec 08: verification of the footprint shortcuts against Kljun's FFP code (local only, cfg$ffp_code_dir).
#   A  ffp_ci_paper() (paper's equations) vs calc_footprint_FFP() (peak = argmax of f_ci, x_r by cumulative f_ci,
#      as ETC) on --n random half-hours of the full-period input table (rows valid in the `none` set, PBLHmod h).
#   C  mean of the per-half-hour fields (ffp_field_halfhour(), as summed in scripts/08b_footprint_run.R) vs
#      calc_footprint_FFP_climatology() on --nclim half-hours of the test month, smooth_data = 0; with
#      smooth_data = 1 also ffp_smooth() vs the code's EBImage::filter2() smoothing, if EBImage is installed.
# Output: <output_dir>/tables/08v_ffp_verify_<site>.csv
#
# Usage: Rscript scripts/08v_ffp_verify.R SITE [--n=1000] [--nclim=10] [--month=2022-09] [--seed=1]

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
site_id <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code.")
n_a <- as.integer(opt("n", "1000"))
n_c <- as.integer(opt("nclim", "10"))
mon <- opt("month", "2022-09")
set.seed(as.integer(opt("seed", "1")))

cfg <- get_cfg()
env <- load_ffp_code(resolve_path(cfg$ffp_code_dir))
if (is.null(env)) cli::cli_abort("Kljun FFP code needed for the verification.")
fp_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "footprint", site_id)
rd_inp <- function(per) {
  p <- file.path(fp_dir, paste0(site_id, "_footprint_inputs", per, ".rds"))
  if (!file.exists(p)) cli::cli_abort("{.path {p}} not found: run scripts/08a_footprint_flags.R first.")
  readRDS(p)
}
out <- list()

# ---- A: crosswind-integrated distances, paper vs code ----
inp <- rd_inp("")
tab <- inp$tab
ok <- which(inp$subsets[["valid__none__pblhmod"]])
idx <- sort(sample(ok, min(n_a, length(ok))))
cli::cli_inform("{site_id} A: {length(idx)} half-hours of {length(ok)} valid (none set, PBLHmod).")
t_a <- Sys.time()
pap <- ffp_ci_paper(tab$zm_m[idx], tab$h_pblhmod[idx], tab$umean[idx], tab$ustar[idx])
cod <- dplyr::bind_rows(lapply(idx, function(i) {
  ffp_ci_code(env, tab$zm_m[i], tab$h_pblhmod[i], tab$umean[i], tab$L[i], tab$sigmav[i], tab$ustar[i], tab$wind_dir[i])
}))
sec_a <- as.numeric(difftime(Sys.time(), t_a, units = "secs"))
both <- cod$flag_err == 0
for (v in c("x_peak", "x_50", "x_70", "x_80", "x_90")) {
  rel <- abs(pap[[v]][both] - cod[[v]][both]) / cod[[v]][both]
  out[[length(out) + 1]] <- dplyr::tibble(part = "A", quantity = v, n = sum(both), n_code_error = sum(!both),
    max_abs_diff = max(abs(pap[[v]][both] - cod[[v]][both])), max_rel_diff = max(rel), median_rel_diff = stats::median(rel))
}
cli::cli_inform("A: {round(sec_a)} s for the code calls ({round(sec_a / length(idx), 3)} s per half-hour).")

# ---- C: climatology from the summed per-half-hour fields vs the climatology function ----
inp_m <- rd_inp(paste0("_", mon))
tm <- inp_m$tab
okm <- which(inp_m$subsets[["valid__none__pblhmod"]])
ic <- sort(sample(okm, min(n_c, length(okm))))
grid <- ffp_grid(c(-1000, 1000, -1000, 1000), 1)
cli::cli_inform("{site_id} C: {length(ic)} half-hours of {mon}.")
t_c <- Sys.time()
fields <- lapply(ic, function(i) ffp_field_halfhour(env, grid, tm$zm_m[i], tm$h_pblhmod[i], tm$umean[i], tm$L[i],
                                                    tm$sigmav[i], tm$ustar[i], tm$wind_dir[i]))
keep <- !vapply(fields, is.null, logical(1))
mean_ours <- Reduce(`+`, fields[keep]) / sum(keep)
sec_b <- as.numeric(difftime(Sys.time(), t_c, units = "secs")) / length(ic)
clim_code <- function(smooth) {
  junk <- utils::capture.output(f <- without_ebimage_notice(env$calc_footprint_FFP_climatology(
    zm = tm$zm_m[ic], z0 = NaN, umean = tm$umean[ic], h = tm$h_pblhmod[ic], ol = tm$L[ic], sigmav = tm$sigmav[ic],
    ustar = tm$ustar[ic], wind_dir = tm$wind_dir[ic], domain = grid$domain, dx = grid$dx, dy = grid$dx, r = NaN,
    rslayer = 1, smooth_data = smooth, pulse = 0)))
  f
}
cc <- clim_code(0)
d0 <- abs(cc$fclim_2d - mean_ours)
out[[length(out) + 1]] <- dplyr::tibble(part = "C", quantity = "fclim_2d (smooth_data = 0)", n = cc$n,
  n_code_error = sum(!keep), max_abs_diff = max(d0), max_rel_diff = max(d0) / max(cc$fclim_2d), median_rel_diff = NA_real_)
m_ours <- ffp_field_metrics(ffp_smooth(mean_ours), grid)$metrics
if (requireNamespace("EBImage", quietly = TRUE)) {
  cs <- clim_code(1)
  d1 <- abs(cs$fclim_2d - ffp_smooth(mean_ours))
  out[[length(out) + 1]] <- dplyr::tibble(part = "C", quantity = "fclim_2d (smooth_data = 1, EBImage)", n = cs$n,
    n_code_error = sum(!keep), max_abs_diff = max(d1), max_rel_diff = max(d1) / max(cs$fclim_2d), median_rel_diff = NA_real_)
} else {
  cli::cli_inform("EBImage not installed: smoothing verified only against the FFT circular convolution (tests).")
}
cli::cli_inform("C: {round(sec_b, 1)} s per half-hour field (B), climatology area_80 {round(m_ours$area_80)} m2.")

res <- dplyr::bind_rows(out)
res$seconds_per_halfhour_code_A <- sec_a / length(idx)
res$seconds_per_halfhour_field_B <- sec_b
print(as.data.frame(res), digits = 4)
tdir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(res, file.path(tdir, paste0("08v_ffp_verify_", site_id, ".csv")))
cli::cli_inform("Written {.path {file.path(tdir, paste0('08v_ffp_verify_', site_id, '.csv'))}}.")
