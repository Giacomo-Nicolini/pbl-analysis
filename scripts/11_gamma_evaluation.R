# Spec 11 D: evaluation of the gamma cases at four sites, with the generic spec-10 comparison (measured_time_base()).
#   References: FR-Gri SIRTA CHM15k (aerosol, hourly; primary set = window convective), FI-Hyy Cloudnet epsilon MLH
#   (turbulence, half-hourly, exact values + censored-tail test; sets convective and all regimes), CH-Dav MeteoSwiss
#   STRATfinder (aerosol, half-hourly; primary set = window convective, G 2026-10-30), BE-Bra ERA5 only.
#   Cases: default (gamma 0.01), gamma0025, gamma_era5clim, gamma_era5daily, gamma_era5clim_zmin1000,
#   gamma_era5clim_bgonly (derived PBLHmod outputs of script 02, standardized in memory).
#   Per site x case: metrics vs the measurement (pblh_metrics(), median log ratio with day-block bootstrap interval,
#   symmetric slope), the FI-Hyy censored-tail test, PBLHmod vs ERA5 tower cell (convective and all regimes);
#   regime shares and stable / neutral heights (gamma_era5clim_bgonly vs gamma_era5clim: the effect of gamma in N of the
#   neutral Eq. B3); clear-day growth curves by apparent solar hour (kt >= 0.65, SW_IN on >= 90 % of the daylight
#   half-hours; G 2026-10-30); gamma diagnostics (daily vs climatology).
#
# Output (<output_dir>): tables/11_gamma_cases_metrics.csv (long), 11_gamma_censored_null.csv,
# 11_gamma_regime_heights.csv, 11_gamma_daily_vs_clim.csv, 11_growth_clearday_<site>.csv; figures/
# 11_gamma_by_month_<sites>.png, 11_gamma_cases_<site>.png, 11_growth_clearday_<site>.png.
#
# Usage: Rscript scripts/11_gamma_evaluation.R [SITE ...] [--n-boot=N]   (default: FR-Gri FI-Hyy CH-Dav BE-Bra)
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
nb_arg <- grep("^--n-boot=", args, value = TRUE)
site_ids <- setdiff(args, nb_arg)
if (length(site_ids) == 0) site_ids <- c("FR-Gri", "FI-Hyy", "CH-Dav", "BE-Bra")
opts <- measured_compare_options()
if (length(nb_arg) == 1) opts$n_boot <- as.integer(sub("^--n-boot=", "", nb_arg))
kt_clear <- 0.65 # clear day (spec 10 item 3b, confirmed G 2026-10-30)
kt_cov <- 0.9
cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
pbl_dir <- file.path(resolve_path(cfg$derived_dir), "pblhmod")
g_dir <- file.path(resolve_path(cfg$derived_dir), "gamma")
tables_dir <- file.path(resolve_path(cfg$output_dir, must_exist = FALSE), "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(dplyr::mutate(x, dplyr::across(dplyr::where(is.double), ~ signif(.x, 3)))),
                               row.names = FALSE)
t0 <- Sys.time()

cases <- c("default", "gamma0025", "gamma_era5clim", "gamma_era5daily", "gamma_era5clim_zmin1000",
           "gamma_era5clim_bgonly")
case_pal <- c(default = "#0072B2", gamma0025 = "#56B4E9", gamma_era5clim = "#009E73", gamma_era5daily = "#E69F00",
              gamma_era5clim_zmin1000 = "#CC79A7", gamma_era5clim_bgonly = "#882255")
refs <- list(`FR-Gri` = list(src = "sirta_mlh_chm15k", type = "aerosol", label = "SIRTA CHM15k"),
             `FI-Hyy` = list(src = "cloudnet_eps_mlh", type = "turbulence", label = "Cloudnet epsilon MLH"),
             `CH-Dav` = list(src = "meteoswiss_stratfinder_mlh", type = "aerosol", label = "MeteoSwiss STRATfinder"))
exact_ref <- function(v) v[v$kind %in% "exact", c("timestamp_utc", "pblh_m")]
mets <- c("n", "bias", "median_log_ratio", "mad_log_ratio", "frac_within_factor2", "r_spearman", "fb", "nrmsd")

# Metrics of one pair table (wide, one row).
metric_row <- function(pt, site_id, case, reference, set) {
  m <- pblh_metrics(pt$mod, pt$ref)
  mb <- median_logratio_boot(pt$mod, pt$ref, pt$timestamp_utc, opts)
  cs <- compression_slope(pt$mod, pt$ref, pt$timestamp_utc, opts, diagnostics = FALSE)
  dplyr::bind_cols(dplyr::tibble(site_id = site_id, case = case, reference = reference, stratum_set = set), m,
                   mb[c("n_days", "mlr_lo", "mlr_hi")],
                   cs[c("slope_sym", "slope_sym_lo", "slope_sym_hi", "ba_slope", "ba_slope_lo", "ba_slope_hi")])
}

wide <- list()
nulls <- list()
reg_h <- list()
for (site_id in site_ids) {
  cli::cli_h1("Gamma cases: {site_id}")
  meta <- sites[sites$site_id == site_id, ]
  std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
  rd <- function(src) read_standard_pblh(file.path(std_dir, paste0(site_id, "__", src, ".csv.gz")))
  raw <- list()
  std <- list()
  for (k in cases) {
    p <- file.path(pbl_dir, paste0(site_id, "_pblhmod_", k, ".csv.gz"))
    if (!file.exists(p)) cli::cli_abort("{.path {p}} not found (script 02).")
    raw[[k]] <- readr::read_csv(p, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
    std[[k]] <- suppressMessages(standardize_pblhmod(raw[[k]], k)$data)
  }
  era <- list(era5_blh_center = rd("era5_blh_center"), era5_blh_bilinear = rd("era5_blh_bilinear"))

  # --- regime shares and stable / neutral heights (bgonly vs clim vs default) ---
  for (k in c("default", "gamma_era5clim", "gamma_era5clim_bgonly")) {
    d <- raw[[k]]
    reg_h[[length(reg_h) + 1]] <- dplyr::bind_rows(lapply(c("stable", "unstable_diag", "convective"), function(r) {
      i <- d$regime %in% r & !is.na(d$pblh_m)
      dplyr::tibble(site_id = site_id, case = k, regime = r, share = mean(d$regime %in% r), n = sum(i),
                    share_B2 = mean(d$method[i] %in% "B2"), share_B3 = mean(d$method[i] %in% "B3"),
                    p25_m = stats::quantile(d$pblh_m[i], 0.25, names = FALSE), median_m = stats::median(d$pblh_m[i]),
                    p75_m = stats::quantile(d$pblh_m[i], 0.75, names = FALSE),
                    median_gamma_bg = stats::median(d$gamma_used[i]), median_gamma_n = stats::median(d$gamma_n_used[i]))
    }))
  }

  # --- PBLHmod vs ERA5 tower cell (ERA5 hours), every case ---
  era_c <- dplyr::rename(era$era5_blh_center, timestamp_utc = timestamp_start_utc)
  for (k in cases) {
    pt <- pair_table(match_halfhours_to_hours(std[[k]], era_c$timestamp_utc), era_c, k)
    x <- pt[pt$regime %in% "convective", ]
    wide[[length(wide) + 1]] <- metric_row(x, site_id, k, "era5_blh_center", "convective")
    for (se in levels(x$season)) {
      wide[[length(wide) + 1]] <- metric_row(x[x$season %in% se, ], site_id, k, "era5_blh_center",
                                             paste0("convective:", se))
    }
    wide[[length(wide) + 1]] <- metric_row(pt, site_id, k, "era5_blh_center", "all_regimes")
  }
  cli::cli_inform("PBLHmod vs ERA5 done ({format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}).")

  # --- measured reference ---
  rf <- refs[[site_id]]
  if (is.null(rf)) next
  meas <- rd(rf$src)
  v_rng <- function(d) range(d$timestamp_start_utc[!is.na(d$pblh_m)])
  p_start <- max(v_rng(std$default)[1], v_rng(meas)[1])
  p_end <- min(v_rng(std$default)[2], v_rng(meas)[2])
  meas <- meas[meas$timestamp_start_utc >= p_start & meas$timestamp_start_utc <= p_end, ]
  cli::cli_inform("{rf$label}: comparison period {format(p_start, '%Y-%m-%d %H:%M')} to {format(p_end, '%Y-%m-%d %H:%M')} UTC.")
  sets <- if (rf$type == "aerosol") {
    list(window_convective = function(pt) pt$regime %in% "convective" & pt$in_window,
         convective = function(pt) pt$regime %in% "convective")
  } else {
    list(convective = function(pt) pt$regime %in% "convective", all_regimes = function(pt) rep(TRUE, nrow(pt)))
  }
  tbs <- list()
  for (k in cases) {
    tb <- measured_time_base(meas, std[[k]], era)
    tbs[[k]] <- tb
    pt <- pair_table(tb$mod, exact_ref(tb$meas), k)
    a <- apparent_solar_hour(pt$timestamp_utc + tb$mid_s, meta$lon)
    pt$in_window <- a >= opts$aerosol_window[1] & a < opts$aerosol_window[2]
    for (s in names(sets)) {
      x <- pt[sets[[s]](pt), ]
      wide[[length(wide) + 1]] <- metric_row(x, site_id, k, rf$src, s)
      for (se in levels(x$season)) {
        wide[[length(wide) + 1]] <- metric_row(x[x$season %in% se, ], site_id, k, rf$src, paste0(s, ":", se))
      }
    }
    if (rf$type == "turbulence") {
      i <- match(as.numeric(tb$mod$timestamp_utc), as.numeric(tb$meas$timestamp_utc))
      nd <- dplyr::tibble(timestamp_utc = tb$mod$timestamp_utc, mod = tb$mod$pblh_m, kind = tb$meas$kind[i],
                          ref = tb$meas$pblh_m[i], z_min = tb$meas$z_min[i], regime = tb$mod$x_regime,
                          season = season_of(tb$mod$timestamp_utc))
      nd <- nd[!is.na(nd$mod) & !is.na(nd$kind), ]
      for (s in c("all_regimes", "convective", "stable")) {
        x <- if (s == "all_regimes") nd else nd[nd$regime %in% s, ]
        nulls[[length(nulls) + 1]] <- dplyr::bind_cols(dplyr::tibble(site_id = site_id, case = k, regime_set = s),
                                                       censored_tail_null(x, opts))
      }
    }
  }
  cli::cli_inform("Measured comparison done ({format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}).")

  # --- clear-day growth curves, every case ---
  icos <- suppressMessages(suppressWarnings(read_icos_site(site_id, cfg, sites, sites_tv, exclusions)))
  d0 <- raw$default
  t_d0 <- as.numeric(as_utc_time(d0$timestamp_start_utc))
  day_of <- as.Date(d0$day_id)[match(as.numeric(icos$timestamp_start_utc), t_d0)]
  kt <- daily_clearness(day_of, icos$met_SW_IN, icos$fn_SW_IN_POT)
  clear <- kt$day_id[!is.na(kt$kt) & kt$coverage >= kt_cov & kt$kt >= kt_clear]
  tb <- tbs$default
  g <- dplyr::tibble(timestamp_utc = tb$meas$timestamp_utc,
                     measured = ifelse(tb$meas$kind %in% "exact", tb$meas$pblh_m, NA_real_))
  at <- function(x) x$pblh_m[match(as.numeric(g$timestamp_utc), as.numeric(x$timestamp_utc))]
  g$era5_blh_center <- at(tb$era$era5_blh_center)
  for (k in cases) g[[k]] <- at(tbs[[k]]$mod)
  g$day_id <- as.Date(d0$day_id)[match(as.numeric(g$timestamp_utc), t_d0)]
  vars <- c("measured", "era5_blh_center", cases)
  gc <- growth_curves(g[g$day_id %in% clear, ], vars, meta$lon, mid_s = tb$mid_s)
  readr::write_csv(dplyr::mutate(gc$summary, site_id = site_id, reference = rf$src, .before = 1),
                   file.path(tables_dir, paste0("11_growth_clearday_", site_id, ".csv")))
  cli::cli_inform("Clear days (kt >= {kt_clear}): {length(clear)}; hours with all values: {nrow(gc$rows)}.")
  gp <- gc$summary[gc$summary$n >= 10, ]
  gp$x <- gp$ast_hour + 0.5
  gp$comparable <- if (rf$type == "aerosol") gp$x >= opts$aerosol_window[1] & gp$x < opts$aerosol_window[2] else TRUE
  pal <- c(measured = "black", era5_blh_center = "#D55E00", case_pal)
  labs_v <- c(measured = rf$label, era5_blh_center = "ERA5 (tower cell)", stats::setNames(cases, cases))
  gp$variable <- factor(gp$variable, levels = names(pal))
  nd <- dplyr::distinct(gc$summary[c("season", "n_clear_days")])
  nd$label <- paste0(nd$n_clear_days, " clear days")
  fig <- ggplot2::ggplot(gp) +
    geom_line_comparable(gp, ggplot2::aes(x, median, colour = variable), c("season", "variable"), "x") +
    ggplot2::geom_text(data = nd, ggplot2::aes(x = 5.2, y = Inf, label = label), hjust = 0, vjust = 1.3, size = 2.6) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, labels = labs_v, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(4, 22, 2), limits = c(5, 21)) +
    ggplot2::labs(x = "Apparent solar time [h, interval midpoint]", y = "Median PBLH / MLH [m]",
                  subtitle = paste0(site_id, ": clear days (kt >= ", kt_clear, "); hours with every value valid",
                                    if (rf$type == "turbulence") "; measured: exact values only (biased low by day)" else "",
                                    "; PBLHmod any regime; n >= 10")) +
    theme_pub() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3))
  pth <- save_figure(fig, paste0("11_growth_clearday_", site_id), 19, 15)
  cli::cli_alert_success("{.path {pth[['png']]}}")
}

# --- tables ---
wide <- dplyr::bind_rows(wide)
long_cols <- c(mets, "n_days", "mlr_lo", "mlr_hi", "slope_sym", "slope_sym_lo", "slope_sym_hi", "ba_slope", "ba_slope_lo",
               "ba_slope_hi")
long <- dplyr::bind_rows(lapply(long_cols, function(m) {
  dplyr::tibble(site_id = wide$site_id, case = wide$case, reference = wide$reference, stratum_set = wide$stratum_set,
                metric = m, value = as.numeric(wide[[m]]))
}))
long$note <- paste("fb > 0: reference deeper (Sinclair2022); mlr_lo/hi, slope_sym_lo/hi, ba_slope_lo/hi: day-block",
                   "bootstrap 95 %; measured references: exact values (FI-Hyy) or quality-filtered aerosol MLH")
readr::write_csv(long, file.path(tables_dir, "11_gamma_cases_metrics.csv"))
reg_h <- dplyr::bind_rows(reg_h)
readr::write_csv(reg_h, file.path(tables_dir, "11_gamma_regime_heights.csv"))
if (length(nulls) > 0) {
  nulls <- dplyr::bind_rows(nulls)
  readr::write_csv(nulls, file.path(tables_dir, "11_gamma_censored_null.csv"))
}

# --- gamma diagnostics: daily vs climatology, by month ---
gd <- dplyr::bind_rows(lapply(site_ids, function(s) {
  d <- readr::read_csv(file.path(g_dir, paste0(s, "_daily.csv")), show_col_types = FALSE)
  cl <- readr::read_csv(file.path(g_dir, paste0(s, "_monthly.csv")), show_col_types = FALSE)
  d$month <- as.integer(format(d$date, "%m"))
  d$clim <- cl$median[match(d$month, cl$month)]
  dplyr::mutate(d, site_id = s)
}))
dvc <- dplyr::summarise(dplyr::group_by(gd[!is.na(gd$gamma), ], site_id), n_days = dplyr::n(),
                        median_abs_rel_dev = stats::median(abs(gamma / clim - 1)),
                        share_outside_25pct = mean(abs(gamma / clim - 1) > 0.25),
                        sd_daily = stats::sd(gamma), sd_clim = stats::sd(clim), .groups = "drop")
readr::write_csv(dvc, file.path(tables_dir, "11_gamma_daily_vs_clim.csv"))
gm <- gd[!is.na(gd$gamma), ]
gm$site_id <- factor(gm$site_id, levels = site_ids)
fig_g <- ggplot2::ggplot(gm, ggplot2::aes(factor(month), gamma)) +
  ggplot2::geom_hline(yintercept = 0.01, linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = 0.0025, linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3, fill = "#009E73", alpha = 0.4) +
  ggplot2::facet_wrap(~site_id, nrow = 2) +
  ggplot2::labs(x = "Month", y = "Daily morning gamma [K m-1] (ERA5, 500-3500 m above the site)",
                subtitle = "Lines: 0.01 (PBLHmod default, solid) and 0.0025 K m-1 (dashed)") +
  theme_pub()
pth <- save_figure(fig_g, paste0("11_gamma_by_month_", paste(site_ids, collapse = "_")), 19, 13)
cli::cli_alert_success("{.path {pth[['png']]}}")

# --- case figures: median log ratio with bootstrap interval ---
for (s in site_ids) {
  w <- wide[wide$site_id == s & !grepl(":", wide$stratum_set, fixed = TRUE), ]
  w$panel <- paste0(ifelse(w$reference == "era5_blh_center", "vs ERA5 tower cell", paste("vs", w$reference)), ", ",
                    w$stratum_set)
  w$case <- factor(w$case, levels = cases)
  fig_c <- ggplot2::ggplot(w, ggplot2::aes(case, median_log_ratio, colour = case)) +
    ratio_lines() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mlr_lo, ymax = mlr_hi), width = 0.2, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~panel) +
    ggplot2::scale_colour_manual(values = case_pal, guide = "none") +
    ggplot2::labs(x = NULL, y = "Median ln(PBLHmod / reference) [-] (day-block bootstrap 95 %)", subtitle = s) +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  pth <- save_figure(fig_c, paste0("11_gamma_cases_", s), 19, 13)
  cli::cli_alert_success("{.path {pth[['png']]}}")
}

# --- report ---
cli::cli_h1("Summary")
cli::cli_h2("Primary sets: median ln ratio [95 %], within factor 2, FB, nRMSD, symmetric slope")
prim <- wide[(wide$reference != "era5_blh_center" & wide$stratum_set %in% c("window_convective", "convective")) |
               (wide$reference == "era5_blh_center" & !grepl(":", wide$stratum_set, fixed = TRUE)), ]
print_tbl(prim[c("site_id", "reference", "stratum_set", "case", "n", "median_log_ratio", "mlr_lo", "mlr_hi",
                 "frac_within_factor2", "fb", "nrmsd", "slope_sym")])
cli::cli_h2("Regime shares and heights (default, ERA5 clim, ERA5 clim in B&G only)")
print_tbl(reg_h)
if (is.data.frame(nulls)) {
  cli::cli_h2("FI-Hyy censored-tail test by case")
  print_tbl(nulls[c("case", "regime_set", "tail", "n", "obs_share", "obs_lo", "obs_hi", "sim_share", "sim_lo", "sim_hi",
                    "sim_bound")])
}
cli::cli_h2("gamma: daily vs climatology")
print_tbl(dvc)
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
