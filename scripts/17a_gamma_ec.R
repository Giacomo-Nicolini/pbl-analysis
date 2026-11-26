# Spec 17 section 1: daily and monthly gamma from the station's own data (EC-derived gamma, gamma_EC). No ERA5.
# Inputs per site: the `default` PBLHmod output of step 02 (<derived_dir>/pblhmod/<site>_pblhmod_default.csv.gz:
# forcing wtheta, T_K, T_source, p_Pa, ustar, L, the convective window, h_diag_m for h0, flags) and the FLUXNET file
# (SW_IN_F, SW_IN_POT for the clear-sky test; H, LE, NETRAD, G for the energy-balance ratio of part D).
# Method A (primary): B&G integration as run_pblh_model() with the mixed-layer temperature; Method B (closed-form
# encroachment, sensitivity). Day selection by flags (nothing deleted); climatology = site x month median of the used
# days (gamma_ec_climatology(), layout of the ERA5 monthly table, read by run_pblh_model() for ec_clim / ec_daily).
# Parameters: config.yml gamma_ec_* (gamma_ec_options_cfg()).
#
# Outputs: <derived_dir>/gamma_ec/<site>_daily.csv, <site>_monthly.csv, <site>_ebr.csv (part D: daytime JJA energy-
#          balance ratio of the reference year; all years where the site has no reference year);
#          <output_dir>/tables/17a_gamma_ec_summary_<site>.csv (flag counts, used days, gamma_A vs gamma_B).
#
# Usage: Rscript scripts/17a_gamma_ec.R SITE [--label=test] [--start-offset-h=H]
#   --label: written into the summary table (e.g. "test" for the PC runs).
#   --start-offset-h: fit start H hours after the first positive heat flux (overrides config gamma_ec_start_offset_h;
#     D1 pre-test, G 2026-11-26).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), c("--label", "--start-offset-h"))
if (length(bad)) cli::cli_abort("Unknown argument{?s}: {.val {bad}}.")
site_id <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code.")
label <- opt("label", "")

cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
meta <- sites[sites$site_id == site_id, ]
if (nrow(meta) != 1) cli::cli_abort("{.val {site_id}} not in the site table.")
derived_dir <- resolve_path(cfg$derived_dir, must_exist = FALSE)
out_dir <- file.path(derived_dir, "gamma_ec")
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
for (d0 in c(out_dir, tables_dir)) dir.create(d0, recursive = TRUE, showWarnings = FALSE)
opts <- gamma_ec_options_cfg(cfg)
if (!is.null(opt("start-offset-h"))) { # D1 pre-test (G 2026-11-26): override of gamma_ec_start_offset_h
  opts <- do.call(gamma_ec_options, utils::modifyList(opts[setdiff(names(opts), c("gamma_min", "gamma_max"))],
                                                      list(start_offset_h = as.numeric(opt("start-offset-h")))))
}
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

cli::cli_h1("EC-derived gamma: {site_id}")
p_def <- file.path(derived_dir, "pblhmod", paste0(site_id, "_pblhmod_default.csv.gz"))
if (!file.exists(p_def)) cli::cli_abort("{.path {p_def}} not found: run step 02 first.")
pm <- readr::read_csv(p_def, show_col_types = FALSE, guess_max = 1e5, progress = FALSE)
pm$timestamp_start_utc <- as_utc_time(pm$timestamp_start_utc)
pm_opts <- pblh_options_cfg(cfg)
fn <- read_icos_fluxnet(find_icos_files(site_id, cfg)$fluxnet, meta$utc_offset_h)
i <- match(as.numeric(pm$timestamp_start_utc), as.numeric(fn$timestamp_start_utc))
col <- function(k) if (k %in% names(fn)) fn[[k]][i] else rep(NA_real_, nrow(pm))
sw_in <- col("fn_SW_IN_F")
sw_pot <- col("fn_SW_IN_POT")
cli::cli_inform(c("{nrow(pm)} half-hours; FLUXNET SW_IN_F on {sum(!is.na(sw_in))}, SW_IN_POT on {sum(!is.na(sw_pot))}.",
                  "T source of the model: {paste(names(table(pm$T_source)), table(pm$T_source), collapse = ', ')}.",
                  "Options: {paste(names(opts), unlist(opts), sep = ' = ', collapse = '; ')}."))
tm <- system.time(daily <- gamma_ec_daily(pm, sw_in, sw_pot, meta, opts, n_sub = pm_opts$n_sub))
cli::cli_inform("Fit time: {round(tm[['elapsed']], 1)} s for {nrow(daily)} days.")
clim <- gamma_ec_climatology(daily, opts$min_days)
w1 <- write_csv_if_changed(daily, file.path(out_dir, paste0(site_id, "_daily.csv")))
w2 <- write_csv_if_changed(dplyr::mutate(clim, site_id = site_id, .before = 1), file.path(out_dir, paste0(site_id, "_monthly.csv")))
# written only when the content changes: step 02 hashes these tables (two passes 02 -> 17a -> 02 converge)
cli::cli_inform("gamma_EC tables: {if (w1 || w2) 'written (content changed)' else 'unchanged (files kept)'}.")

# Part D: daytime summer energy-balance ratio on the reference year (diagnostic, no correction)
ry <- tryCatch(read_reference_years(resolve_path(cfg$reference_years_file), site_id)$year, error = function(e) NA)
yr_loc <- as.integer(format(fn$timestamp_start_utc + 3600 * meta$utc_offset_h, "%Y", tz = "UTC"))
sel <- if (length(ry) == 1 && !is.na(ry)) yr_loc == as.integer(ry) else rep(TRUE, nrow(fn))
eb <- dplyr::mutate(energy_balance_ratio(fn[sel, ]), site_id = site_id,
                     period = if (length(ry) == 1 && !is.na(ry)) paste0("reference year ", ry) else "all years",
                     .before = 1)
write_csv_if_changed(eb, file.path(out_dir, paste0(site_id, "_ebr.csv")))

# Summary
fc <- gamma_ec_flag_counts(daily)
u <- daily[daily$used, ]
rat <- u$gamma_B / u$gamma_A
summ <- dplyr::tibble(site_id = site_id, label = label, n_days = nrow(daily), n_fitted = sum(!is.na(daily$gamma_A)),
                      n_used = nrow(u), gamma_A_median = stats::median(u$gamma_A), gamma_A_p25 = q_fin(u$gamma_A, 0.25),
                      gamma_A_p75 = q_fin(u$gamma_A, 0.75), ratio_B_A_median = med_fin(rat),
                      r2_A_median = med_fin(u$r2_A), n_tsonic_rows = sum(daily$n_tsonic),
                      n_bound_low = sum(daily$gamma_A_bound %in% "low"), n_bound_high = sum(daily$gamma_A_bound %in% "high"),
                      share_clear_unknown = mean(grepl("clear_unknown", daily$flags)), months_with_value = sum(!is.na(clim$median)),
                      ebr = eb$ebr, ebr_no_g = eb$ebr_no_g)
summ <- dplyr::bind_cols(summ, as.data.frame(stats::setNames(as.list(fc$n_days), paste0("n_", fc$flag)))[, -(1:2)])
readr::write_csv(summ, file.path(tables_dir, paste0("17a_gamma_ec_summary_", site_id, ".csv")))

cli::cli_h2("Days and flags (a day can carry several flags)")
print_tbl(fc)
cli::cli_h2("gamma_A vs gamma_B on the used days")
cli::cli_inform("Used days: {nrow(u)}; median gamma_A {signif(summ$gamma_A_median, 3)} K m-1 (p25 {signif(summ$gamma_A_p25, 3)}, p75 {signif(summ$gamma_A_p75, 3)}); median gamma_B / gamma_A {signif(summ$ratio_B_A_median, 3)}; median R2_A {signif(summ$r2_A_median, 3)}.")
cli::cli_h2("Monthly climatology [K m-1] (NA: fewer than {opts$min_days} used days)")
print_tbl(clim)
cli::cli_h2("Energy-balance ratio (part D)")
print_tbl(eb)
cli::cli_alert_success("Written {.path {out_dir}}/{site_id}_{{daily,monthly,ebr}}.csv and 17a_gamma_ec_summary_{site_id}.csv")
