# Spec 11 A: free-atmosphere gamma from the ERA5 pressure levels (tower cell), per day (morning profile: the
# 6-hourly profile closest to and not after sunrise) and as a site x month climatology.
# Fit layer: default 500-3500 m above the site (gamma_options(), G 2026-10-30) and the sensitivity 1000-3500 m
# (case gamma_era5clim_zmin1000). Also reported: the same fit above the ERA5 cell orography (gamma_orog).
# Output: <derived_dir>/gamma/<site>_daily.csv, <site>_monthly.csv (default layer) and <site>_zmin1000_daily.csv,
# <site>_zmin1000_monthly.csv; tables/11a_gamma_summary_<sites>.csv (levels used, flags, R2, by season).
#
# Usage: Rscript scripts/11a_gamma_era5.R [SITE ...]   (default: BE-Bra FI-Hyy FR-Gri CH-Dav)
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- get_cfg()
site_ids <- commandArgs(trailingOnly = TRUE)
if (length(site_ids) == 0) site_ids <- c("BE-Bra", "FI-Hyy", "FR-Gri", "CH-Dav")
sites <- read_sites(resolve_path(cfg$sites_file))
g_dir <- file.path(resolve_path(cfg$derived_dir, must_exist = FALSE), "gamma")
dir.create(g_dir, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(dplyr::mutate(x, dplyr::across(dplyr::where(is.double), ~ signif(.x, 3)))),
                               row.names = FALSE)
layers <- list(default = gamma_options(), zmin1000 = gamma_options(zmin_agl = 1000))

summ <- list()
for (site_id in site_ids) {
  cli::cli_h1("ERA5 gamma: {site_id}")
  meta <- sites[sites$site_id == site_id, ]
  pl <- read_era5_pressure_levels(find_era5_files(site_id, cfg, "pressure_level", sites = sites), cells = "center")
  st <- read_era5_static(site_id, cfg, sites)
  orog <- st$orog_m[st$cell_row == 0 & st$cell_col == 0]
  cli::cli_inform("Site altitude {meta$altitude_m} m; ERA5 tower-cell orography {round(orog, 1)} m (difference {round(meta$altitude_m - orog, 1)} m).")
  lv <- dplyr::summarise(dplyr::group_by(pl, p_hPa), z_agl_p05 = stats::quantile(zg_m - meta$altitude_m, 0.05),
                         z_agl_p50 = stats::median(zg_m - meta$altitude_m),
                         z_agl_p95 = stats::quantile(zg_m - meta$altitude_m, 0.95), .groups = "drop")
  print_tbl(lv[lv$p_hPa >= 500, ])
  days <- seq(as.Date(min(pl$timestamp_utc)), as.Date(max(pl$timestamp_utc)), by = "day")
  for (ly in names(layers)) {
    o <- layers[[ly]]
    prof <- gamma_profiles(pl, meta$altitude_m, orog, o)
    daily <- gamma_daily(prof, days, meta$lat, meta$lon, o)
    clim <- gamma_climatology(daily)
    sfx <- if (ly == "default") "" else paste0("_", ly)
    readr::write_csv(daily, file.path(g_dir, paste0(site_id, sfx, "_daily.csv")))
    readr::write_csv(dplyr::mutate(clim, site_id = site_id, .before = 1), file.path(g_dir, paste0(site_id, sfx, "_monthly.csv")))
    cli::cli_h2("Layer {ly}: {o$zmin_agl}-{o$zmax_agl} m above the site; {nrow(daily)} days")
    cli::cli_inform("Flags: {paste(names(table(daily$flag)), table(daily$flag), sep = ' ', collapse = ', ')}; profile hour (UTC): {paste(names(table(format(daily$profile_utc, '%H'))), table(format(daily$profile_utc, '%H')), sep = ' ', collapse = ', ')}.")
    se <- as.character(season_of(as.POSIXct(daily$date, tz = "UTC")))
    by_lev <- dplyr::summarise(dplyr::group_by(dplyr::tibble(levels = daily$levels, n = daily$n_levels, g = daily$gamma,
                                                             raw = daily$gamma_raw, r2 = daily$r2), levels, n),
                               days = dplyr::n(), share = NA_real_, gamma_median = stats::median(g, na.rm = TRUE),
                               raw_median = stats::median(raw, na.rm = TRUE), r2_median = stats::median(r2, na.rm = TRUE),
                               n_na = sum(is.na(g)), .groups = "drop")
    by_lev$share <- by_lev$days / sum(by_lev$days)
    cli::cli_inform("By levels used:")
    print_tbl(by_lev)
    by_se <- dplyr::bind_rows(lapply(c("all", "DJF", "MAM", "JJA", "SON"), function(s) {
      i <- if (s == "all") rep(TRUE, nrow(daily)) else se == s
      g <- daily$gamma[i]
      dplyr::tibble(site_id = site_id, layer = ly, season = s, n_days = sum(i), n_valid = sum(!is.na(g)),
                    share_2_levels = mean(daily$n_levels[i] == 2), share_implausible = mean(daily$flag[i] == "implausible"),
                    share_too_few = mean(daily$flag[i] == "too_few_levels"),
                    gamma_p25 = stats::quantile(g, 0.25, na.rm = TRUE, names = FALSE),
                    gamma_median = stats::median(g, na.rm = TRUE),
                    gamma_p75 = stats::quantile(g, 0.75, na.rm = TRUE, names = FALSE),
                    r2_median = stats::median(daily$r2[i], na.rm = TRUE),
                    gamma_orog_median = stats::median(daily$gamma_orog[i], na.rm = TRUE),
                    ratio_to_0.01 = stats::median(g, na.rm = TRUE) / 0.01)
    }))
    summ[[length(summ) + 1]] <- by_se
    cli::cli_inform("By season (gamma [K m-1]):")
    print_tbl(by_se[setdiff(names(by_se), c("site_id", "layer"))])
    cli::cli_inform("Monthly climatology (median [p25, p75], n days):")
    print_tbl(clim)
  }
}
summ <- dplyr::bind_rows(summ)
out <- file.path(tables_dir, paste0("11a_gamma_summary_", paste(site_ids, collapse = "_"), ".csv"))
readr::write_csv(summ, out)
cli::cli_alert_success("{.path {out}}")
