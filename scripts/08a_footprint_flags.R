# Spec 08 (phase 5a): footprint input table, validity flags and coverage for one station, over the requested period
# (default: the full reference period). No FFP run: the FFP code's own checks are reproduced as conditions.
# PBLH cases: pblhmod, era5 (tower cell, interpolated to the half-hour midpoint), era5_bilinear, pblhmod_x0.8,
# pblhmod_x1.2, measured (exact MLH, qc 0/2) and bound (qc 3 lower bounds) where a Cloudnet MLH exists; oper with
# --with-oper (internal use only: outputs under output/*/oper_internal/).
# Validity sets: kljun (Kljun et al. 2015, Eq. 27), kljun_norsl (Eq. 27 without the roughness condition), etc (ETC
# pipeline, stability from zm/L), none; all include the code checks and zm < h. Subsets: valid, valid_qc (NEE QC good
# in FLUXES), valid_qc_used (and NEE not gap-filled in FLUXNET).
# Output: <derived_dir>/footprint/<site>/<site>_footprint_inputs[_<period>].rds (input table + flags, for 08b)
#         <output_dir>/tables/08_coverage_<site>[_<period>].csv, 08_etc_multiplier_<site>.csv,
#         08_ustar_check_<site>.csv (ONEFlux VUT thresholds vs FLUXNET NEE QC)
#
# Usage: Rscript scripts/08a_footprint_flags.R SITE [--start=YYYY-MM-DD] [--end=YYYY-MM-DD] [--with-oper]

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
with_oper <- "--with-oper" %in% args
site_id <- grep("^--", args, value = TRUE, invert = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one station code.")

cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
meta <- sites[sites$site_id == site_id, ]
derived_dir <- resolve_path(cfg$derived_dir, must_exist = FALSE)
std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
fp_dir <- file.path(derived_dir, "footprint", site_id)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
if (with_oper) tables_dir <- file.path(tables_dir, "oper_internal")
for (d0 in c(fp_dir, tables_dir)) dir.create(d0, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(x), row.names = FALSE, digits = 4)

cli::cli_h1("Footprint flags and coverage: {site_id}")
d <- read_icos_site(site_id, cfg, sites, sites_tv, exclusions)
pm <- readr::read_csv(file.path(derived_dir, "pblhmod", paste0(site_id, "_pblhmod_default.csv.gz")),
                      show_col_types = FALSE, guess_max = 1e5, progress = FALSE)
rd <- function(src) {
  p <- file.path(std_dir, paste0(site_id, "__", src, ".csv.gz"))
  if (file.exists(p)) suppressMessages(read_standard_pblh(p)) else NULL
}
cases <- footprint_h_cases(d$timestamp_start_utc, rd("pblhmod_default"), rd("era5_blh_center"), rd("era5_blh_bilinear"),
                           measured = rd("cloudnet_eps_mlh"), oper = if (with_oper) rd("icos_oper") else NULL)
tab <- footprint_input_table(d, pm, cases)

# ETC multiplier: once per site, full period, rows with zm and z0 (independent of h and of the requested period)
m_etc <- etc_z0_multiplier(tab$zm_m, tab$z0_m)
mt <- dplyr::tibble(site_id = site_id, m_etc = m_etc, n_rows = sum(!is.na(tab$zm_m) & !is.na(tab$z0_m)),
                    zm_over_z0_median = stats::median(tab$zm_m / tab$z0_m, na.rm = TRUE),
                    n_rejected = sum(tab$zm_m < m_etc * tab$z0_m, na.rm = TRUE))
readr::write_csv(mt, file.path(tables_dir, paste0("08_etc_multiplier_", site_id, ".csv")))
cli::cli_inform("ETC roughness multiplier m = {m_etc} (zm/z0 median {round(mt$zm_over_z0_median, 2)}; {mt$n_rejected} of {mt$n_rows} rows rejected).")

# Period
t0 <- as.POSIXct(opt("start", format(min(tab$timestamp_start_utc), "%Y-%m-%d")), tz = "UTC")
t1 <- as.POSIXct(opt("end", format(max(tab$timestamp_start_utc), "%Y-%m-%d")), tz = "UTC") + 86400
per <- if (is.null(opt("start")) && is.null(opt("end"))) "" else paste0("_", format(t0, "%Y-%m"),
                                                                       if (format(t0, "%Y-%m") != format(t1 - 1, "%Y-%m")) paste0("_", format(t1 - 1, "%Y-%m")) else "")
tab <- tab[tab$timestamp_start_utc >= t0 & tab$timestamp_start_utc < t1, ]
cli::cli_inform("Period {format(t0, '%Y-%m-%d')} to {format(t1 - 1, '%Y-%m-%d')}: {nrow(tab)} half-hours; {sum(tab$excluded)} in exclusion periods.")

fl <- footprint_flags_all(tab, m_etc)
cov <- footprint_coverage(tab, fl)
readr::write_csv(cov, file.path(tables_dir, paste0("08_coverage_", site_id, per, ".csv")))
saveRDS(list(site_id = site_id, tab = tab, flags = fl$flags, subsets = fl$subsets, m_etc = m_etc, period = per,
             t0 = t0, t1 = t1), file.path(fp_dir, paste0(site_id, "_footprint_inputs", per, ".rds")))

cli::cli_h2("Coverage (all strata), per case and set")
print_tbl(cov[cov$stratum == "all", c("case", "set", "n_rows", "n_with_h", "n_inputs", "n_valid", "n_valid_qc",
                                       "n_valid_qc_used", "n_h_lt_zm", "n_h_le_10")])
cli::cli_h2("Rejections per condition (rows with h and inputs), pblhmod and era5")
print_tbl(cov[cov$stratum == "all" & cov$case %in% c("pblhmod", "era5") & cov$set == "none",
              c("case", grep("^fail_", names(cov), value = TRUE))])

cli::cli_h2("ONEFlux u* thresholds vs FLUXNET NEE QC")
uc <- tryCatch({
  fr <- read_oneflux_bif(find_oneflux_bif(site_id, cfg), "GRP_FLUX_REF", hh_path = find_icos_files(site_id, cfg)$fluxnet)
  thr <- oneflux_ustar_ref(fr, "VUT")
  print_tbl(thr)
  x <- oneflux_ustar_check(d, thr)
  readr::write_csv(dplyr::bind_cols(dplyr::tibble(site_id = site_id), x), file.path(tables_dir, paste0("08_ustar_check_", site_id, ".csv")))
  print_tbl(x[is.na(x$year) | x$period == "all", ])
}, error = function(e) cli::cli_alert_danger("u* check not possible: {conditionMessage(e)}"))
