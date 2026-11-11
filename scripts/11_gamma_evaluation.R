# Spec 11 D and amendment 1 (G 2026-11-11): evaluation of the gamma cases at the measured sites, with the generic
# spec-10 comparison (measured_time_base()), for the decision on the EC-only constant gamma (spec 09 decision 3b;
# evidence only, the default stays 0.01).
#   Reference of each site = its primary measured source (config/measured_sites.csv, measured_reference(), as
#   scripts/10_measured_compare.R): aerosol (SIRTA CHM15k hourly; MeteoSwiss STRATfinder / CABAM half-hourly) on the
#   window set (apparent solar time 10-16 h, PBLHmod convective) and all convective intervals; turbulence (FI-Hyy
#   Cloudnet epsilon MLH, exact values + censored-tail test) on the convective set and all regimes (where the
#   distance class allows). Sites without a primary source (BE-Bra) are compared with ERA5 only. The standard files
#   of the measured sources already exclude config/measured_exclusions.csv; tables of the decision carry the distance
#   class and the borderline flag.
#   Cases (--cases): default (PBLHmod, gamma = config pblh_gamma_default), clim (PBLHmod-clim, ERA5 climatology),
#   gamma_med (network median of the site medians of the daily ERA5 gamma), gamma_site (the site's median),
#   gamma0025 (legacy ETC value) by default; also possible: gamma001, gamma_era5daily, gamma_era5clim_zmin1000,
#   gamma_era5clim_bgonly (derived PBLHmod outputs of script 02, standardized in memory; a missing case aborts:
#   script 02 runs them, gamma_med only after scripts/09_collect.R wrote the network gamma).
#   Per site x case: metrics vs the measurement (pblh_metrics(), median log ratio with day-block bootstrap interval,
#   symmetric slope), all seasons and by season, the FI-Hyy censored-tail test, PBLHmod vs ERA5 tower cell
#   (convective, by season, and all regimes); ERA5 tower cell vs the measurement on the same sets (decision table
#   and season table only); regime shares and stable / neutral heights (gamma_era5clim_bgonly vs clim: the effect
#   of gamma in N of the neutral Eq. B3); clear-day growth curves by apparent solar hour (kt >= 0.65, SW_IN on
#   >= 90 % of the daylight half-hours; G 2026-10-30); gamma diagnostics (daily vs climatology).
#   Network context (PBLHmod, PBLHmod-clim, gamma_med, gamma_site vs ERA5 at every site): pairs 1, 6, 7, 8 of
#   tables/09_network_era5_metrics.csv (scripts/09_collect.R); not repeated here.
#
# Output (<output_dir>): tables/11_gamma_cases_metrics.csv (long), 11_gamma_cases_season.csv, 11_gamma_decision.csv
# (one row per site), 11_gamma_censored_null.csv, 11_gamma_regime_heights.csv, 11_gamma_daily_vs_clim.csv,
# 11_growth_clearday_<site>.csv; figures/11_gamma_decision.png, 11_gamma_by_month_<sites>.png,
# 11_gamma_cases_<site>.png, 11_growth_clearday_<site>.png (main: PBLHmod, PBLHmod-clim, ERA5, measured) and
# 11s_ref_gamma_cases_<site>.png, 11s_ref_growth_clearday_<site>.png (supplement: gamma_med, gamma_site and the
# reference gamma values; G 2026-10-30, 2026-11-11).
# Robustness (2026-11-11, FI-Lom): one site never stops the job. Each site runs inside an error handler; the
# comparison period is the overlap of the valid PBLHmod default and measured records (comparison_period()), and an
# empty or inverted overlap writes rows with n = 0; the status and reason of each site are columns of
# 11_gamma_decision.csv and printed at the end; the log lists the valid range of every input. Exit status 1 if any
# site has status error / error_growth (after writing every output), else 0.
# With the sites and cases of spec 11 D (FR-Gri FI-Hyy CH-Dav BE-Bra, --cases=default,clim,gamma001,gamma0025,
# gamma_era5daily,gamma_era5clim_zmin1000,gamma_era5clim_bgonly) the tables of spec 11 D are unchanged.
#
# Usage: Rscript scripts/11_gamma_evaluation.R [SITE ...] [--cases=a,b,...] [--n-boot=N]
#   Default sites: BE-Bra and the sites of config/measured_run_sites.txt; default cases: default, clim, gamma_med,
#   gamma_site, gamma0025. Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
nb_arg <- grep("^--n-boot=", args, value = TRUE)
case_arg <- grep("^--cases=", args, value = TRUE)
site_ids <- setdiff(args, c(nb_arg, case_arg))
bad_arg <- grep("^--", site_ids, value = TRUE)
if (length(bad_arg) > 0) cli::cli_abort("Unknown option{?s} {.val {bad_arg}}.")
opts <- measured_compare_options()
if (length(nb_arg) == 1) opts$n_boot <- as.integer(sub("^--n-boot=", "", nb_arg))
kt_clear <- 0.65 # clear day (spec 10 item 3b, confirmed G 2026-10-30)
kt_cov <- 0.9
# low_n: fewer days than this in the bootstrap of a set (spec 11 amendment 1; no minimum existed: flagged for G)
min_days <- 10L
cfg <- get_cfg()
if (length(site_ids) == 0) {
  ms_sites <- trimws(readLines(here::here("config", "measured_run_sites.txt"), warn = FALSE))
  site_ids <- c("BE-Bra", ms_sites[nzchar(ms_sites) & !startsWith(ms_sites, "#")])
}
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

all_cases <- c("default", "clim", "gamma001", "gamma0025", "gamma_era5daily", "gamma_era5clim_zmin1000",
               "gamma_era5clim_bgonly", "gamma_med", "gamma_site")
cases <- if (length(case_arg) == 1) strsplit(sub("^--cases=", "", case_arg), ",")[[1]] else
  c("default", "clim", "gamma_med", "gamma_site", "gamma0025")
bad <- setdiff(cases, all_cases)
if (length(bad) > 0) cli::cli_abort("Unknown case{?s} {.val {bad}} (known: {.val {all_cases}}).")
if (!"default" %in% cases) cli::cli_abort("Case {.val default} is required (regimes and the time base).")
# Figures (G 2026-10-30, 2026-11-11): main = PBLHmod (default) and PBLHmod-clim (clim) with ERA5 and the measurement;
# supplementary = gamma_med, gamma_site and the reference gamma values (gamma001, gamma0025), same layout; the other
# cases in tables only.
styled <- c("default", "clim", "gamma_med", "gamma_site", "gamma001", "gamma0025")
case_pal <- stats::setNames(style_of(paste0("pblhmod_", styled)), styled)
case_lab <- stats::setNames(style_of(paste0("pblhmod_", styled), "label"), styled)
ref_cases <- intersect(c("gamma_med", "gamma_site", "gamma001", "gamma0025"), cases)
fig_sets <- list(main = list(cases = intersect(c("default", "clim"), cases), prefix = "11_", title = NULL),
                 ref = list(cases = ref_cases, prefix = "11s_ref_",
                            title = if (any(c("gamma_med", "gamma_site") %in% ref_cases)) {
                              "γ network and site medians, reference γ values"
                            } else {
                              "reference γ values"
                            }))
fig_sets <- Filter(function(fs) length(fs$cases) > 0, fig_sets)
# Reference per site: the primary measured source (NULL = ERA5 only)
ms <- read_measured_sites(resolve_path(cfg$measured_sites_file))
refs <- stats::setNames(lapply(site_ids, function(s) measured_reference(ms, s)), site_ids)
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

na_t <- as.POSIXct(NA_real_, tz = "UTC")
empty_pt <- dplyr::tibble(timestamp_utc = na_t[0], mod = numeric(), ref = numeric())
fmt_t <- function(t) ifelse(is.na(t), "none", format(t, "%Y-%m-%d %H:%M", tz = "UTC"))

# One site: its rows of every table, in the environment o of the caller (wide, wide_em = ERA5 tower cell vs the
# measurement for the decision and season tables only, nulls, reg_h, period, status, reason, stage). Status "ok",
# "era5_only" (no primary measured source), "no_overlap" (empty or inverted overlap of the valid PBLHmod and measured
# records: rows with n = 0 and the reason). Errors are caught by the loop below, which still sees o: one site never
# stops the job (G 2026-11-11).
eval_site <- function(site_id, o) {
  o$wide <- list()
  o$wide_em <- list()
  o$nulls <- list()
  o$reg_h <- list()
  o$period <- NULL
  o$status <- "ok"
  o$reason <- NA_character_
  push <- function(what, x) o[[what]][[length(o[[what]]) + 1]] <- x
  meta <- sites[sites$site_id == site_id, ]
  if (nrow(meta) != 1) cli::cli_abort("{site_id}: {nrow(meta)} rows in the site table.")
  std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
  rd <- function(src) read_standard_pblh(file.path(std_dir, paste0(site_id, "__", src, ".csv.gz")))
  raw <- list()
  std <- list()
  for (k in cases) {
    p <- file.path(pbl_dir, paste0(site_id, "_pblhmod_", k, ".csv.gz"))
    if (!file.exists(p)) {
      cli::cli_abort(c("{.path {p}} not found (script 02).",
                       i = if (k == "gamma_med") "gamma_med needs the network gamma of scripts/09_collect.R (second pass)."))
    }
    raw[[k]] <- readr::read_csv(p, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
    std[[k]] <- suppressMessages(standardize_pblhmod(raw[[k]], k)$data)
  }
  era <- list(era5_blh_center = rd("era5_blh_center"), era5_blh_bilinear = rd("era5_blh_bilinear"))

  # --- regime shares and stable / neutral heights (bgonly vs clim vs default; gamma_med, gamma_site) ---
  for (k in intersect(c("default", "clim", "gamma0025", "gamma_era5clim_bgonly", "gamma_med", "gamma_site"), cases)) {
    d <- raw[[k]]
    push("reg_h", dplyr::bind_rows(lapply(c("stable", "unstable_diag", "convective"), function(r) {
      i <- d$regime %in% r & !is.na(d$pblh_m)
      dplyr::tibble(site_id = site_id, case = k, regime = r, share = mean(d$regime %in% r), n = sum(i),
                    share_B2 = mean(d$method[i] %in% "B2"), share_B3 = mean(d$method[i] %in% "B3"),
                    p25_m = stats::quantile(d$pblh_m[i], 0.25, names = FALSE), median_m = stats::median(d$pblh_m[i]),
                    p75_m = stats::quantile(d$pblh_m[i], 0.75, names = FALSE),
                    median_gamma_bg = stats::median(d$gamma_used[i]), median_gamma_n = stats::median(d$gamma_n_used[i]))
    })))
  }

  # --- PBLHmod vs ERA5 tower cell (ERA5 hours), every case ---
  era_c <- dplyr::rename(era$era5_blh_center, timestamp_utc = timestamp_start_utc)
  for (k in cases) {
    pt <- pair_table(match_halfhours_to_hours(std[[k]], era_c$timestamp_utc), era_c, k)
    x <- pt[pt$regime %in% "convective", ]
    push("wide", metric_row(x, site_id, k, "era5_blh_center", "convective"))
    for (se in levels(x$season)) {
      push("wide", metric_row(x[x$season %in% se, ], site_id, k, "era5_blh_center", paste0("convective:", se)))
    }
    push("wide", metric_row(pt, site_id, k, "era5_blh_center", "all_regimes"))
  }
  cli::cli_inform("PBLHmod vs ERA5 done ({format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}).")

  # --- measured reference ---
  rf <- refs[[site_id]]
  if (is.null(rf)) {
    cli::cli_inform("{site_id}: no primary measured source in {.file measured_sites.csv}: ERA5 only.")
    o$status <- "era5_only"
    return(o)
  }
  meas_all <- rd(rf$src)
  # time step on the full file: the period cut below may leave no row
  step_s <- measured_step_s(meas_all)
  # valid record of every input (log only; the period is the overlap of PBLHmod default and the measurement)
  rng_in <- c(stats::setNames(std, cases), list(era5_blh_center = era$era5_blh_center), stats::setNames(list(meas_all), rf$src))
  cp_all <- do.call(comparison_period, rng_in)
  cli::cli_inform("Valid records of the inputs:")
  print(as.data.frame(dplyr::mutate(cp_all$ranges, first_valid_utc = fmt_t(first_valid_utc),
                                    last_valid_utc = fmt_t(last_valid_utc))), row.names = FALSE)
  cp <- do.call(comparison_period, stats::setNames(list(std$default, meas_all), c("pblhmod_default", rf$src)))
  o$period <- dplyr::tibble(site_id = site_id, period_start_utc = cp$start, period_end_utc = cp$end)
  set_fun <- list(window_convective = function(pt) pt$regime %in% "convective" & pt$in_window,
                  convective = function(pt) pt$regime %in% "convective",
                  all_regimes = function(pt) rep(TRUE, nrow(pt)))
  sets <- set_fun[rf$sets]
  if (!cp$ok) {
    # empty or inverted overlap (FI-Lom: EC record from 2024-11, CABAM record ends 2024-05): one row per case x set
    # with n = 0 (no season rows), ERA5 vs measured likewise; the reason goes to the decision table
    cli::cli_alert_warning("{site_id}: {cp$reason}; rows with n = 0 written, no growth curve.")
    for (k in cases) for (s in names(sets)) push("wide", metric_row(empty_pt, site_id, k, rf$src, s))
    for (s in names(sets)) push("wide_em", metric_row(empty_pt, site_id, "era5_blh_center", rf$src, s))
    o$status <- "no_overlap"
    o$reason <- cp$reason
    return(o)
  }
  p_start <- cp$start
  p_end <- cp$end
  meas <- meas_all[meas_all$timestamp_start_utc >= p_start & meas_all$timestamp_start_utc <= p_end, ]
  cli::cli_inform("{style_of(rf$src, 'label')} ({rf$provider_site}, {rf$distance_km} km, {rf$distance_class}{if (rf$borderline) ', borderline' else ''}): comparison period {format(p_start, '%Y-%m-%d %H:%M')} to {format(p_end, '%Y-%m-%d %H:%M')} UTC.")
  tbs <- list()
  for (k in cases) {
    tb <- measured_time_base(meas, std[[k]], era, step_s = step_s)
    tbs[[k]] <- tb
    pt <- pair_table(tb$mod, exact_ref(tb$meas), k)
    a <- apparent_solar_hour(pt$timestamp_utc + tb$mid_s, meta$lon)
    pt$in_window <- a >= opts$aerosol_window[1] & a < opts$aerosol_window[2]
    for (s in names(sets)) {
      x <- pt[sets[[s]](pt), ]
      if (nrow(x) == 0) cli::cli_alert_warning("{site_id} {k} [{s}]: no pair; row with n = 0 written.")
      push("wide", metric_row(x, site_id, k, rf$src, s))
      for (se in levels(x$season)) {
        push("wide", metric_row(x[x$season %in% se, ], site_id, k, rf$src, paste0(s, ":", se)))
      }
    }
    if (k == "default") {
      # ERA5 tower cell vs the measurement on the same sets (PBLHmod regime of the ERA5 base; regimes do not depend
      # on gamma), as pair 2 of scripts/10_measured_compare.R
      pe <- pair_table(tb$era$era5_blh_center, exact_ref(tb$era_meas), "era5_blh_center")
      a <- apparent_solar_hour(pe$timestamp_utc + tb$era_mid_s, meta$lon)
      pe$in_window <- a >= opts$aerosol_window[1] & a < opts$aerosol_window[2]
      for (s in names(sets)) {
        x <- pe[sets[[s]](pe), ]
        push("wide_em", metric_row(x, site_id, "era5_blh_center", rf$src, s))
        for (se in levels(x$season)) {
          push("wide_em", metric_row(x[x$season %in% se, ], site_id, "era5_blh_center", rf$src, paste0(s, ":", se)))
        }
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
        push("nulls", dplyr::bind_cols(dplyr::tibble(site_id = site_id, case = k, regime_set = s),
                                       censored_tail_null(x, opts)))
      }
    }
  }
  cli::cli_inform("Measured comparison done ({format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}).")

  # --- clear-day growth curves, every case ---
  o$stage <- "growth"
  icos <-suppressMessages(suppressWarnings(read_icos_site(site_id, cfg, sites, sites_tv, exclusions)))
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
  g_clear <- g[g$day_id %in% clear, ]
  g_ok <- stats::complete.cases(g_clear[vars])
  if (!any(g_ok)) {
    # e.g. a short overlap in the polar winter: no growth curve, no table, no figure
    cli::cli_alert_warning("{site_id}: {length(clear)} clear day{?s}, no hour with every value valid: growth curves skipped.")
    return(o)
  }
  gc <- growth_curves(g_clear, vars, meta$lon, mid_s = tb$mid_s)
  readr::write_csv(dplyr::mutate(gc$summary, site_id = site_id, reference = rf$src, .before = 1),
                   file.path(tables_dir, paste0("11_growth_clearday_", site_id, ".csv")))
  cli::cli_inform("Clear days (kt >= {kt_clear}): {length(clear)}; hours with all values: {nrow(gc$rows)}.")
  gp <- gc$summary[gc$summary$n >= 10, ]
  if (nrow(gp) == 0) {
    cli::cli_alert_warning("{site_id}: no apparent solar hour with n >= 10: growth figures skipped.")
    return(o)
  }
  gp$x <- gp$ast_hour + 0.5
  gp$comparable <- if (rf$type == "aerosol") gp$x >= opts$aerosol_window[1] & gp$x < opts$aerosol_window[2] else TRUE
  nd <- dplyr::distinct(gc$summary[c("season", "n_clear_days")])
  nd$label <- paste0(nd$n_clear_days, " clear days")
  for (fs in fig_sets) {
    pal <- c(measured = unname(style_of(rf$src)), era5_blh_center = unname(style_of("era5_blh_center")), case_pal[fs$cases])
    labs_v <- c(measured = unname(style_of(rf$src, "label")), era5_blh_center = unname(style_of("era5_blh_center", "label")),
                case_lab[fs$cases])
    gs <- gp[gp$variable %in% names(pal), ]
    gs$variable <- factor(gs$variable, levels = names(pal))
    fig <- ggplot2::ggplot(gs) +
      geom_line_comparable(gs, ggplot2::aes(x, median, colour = variable), c("season", "variable"), "x") +
      ggplot2::geom_text(data = nd, ggplot2::aes(x = 5.2, y = Inf, label = label), hjust = 0, vjust = 1.3, size = 2.6) +
      ggplot2::facet_wrap(~season, nrow = 2) +
      ggplot2::scale_colour_manual(values = pal, labels = labs_v, name = NULL) +
      ggplot2::scale_x_continuous(breaks = seq(4, 22, 2), limits = c(5, 21)) +
      ggplot2::labs(x = "Apparent solar time [h, interval midpoint]", y = "Median PBLH / MLH [m]",
                    title = fs$title,
                    subtitle = paste0(site_id, ": clear days (kt >= ", kt_clear, "); hours with every value valid",
                                      " (all cases)",
                                      if (rf$type == "turbulence") "; measured: exact values only (biased low by day)" else "",
                                      "; PBLHmod any regime; n >= 10")) +
      theme_pblh() +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
    pth <- save_figure(fig, paste0(fs$prefix, "growth_clearday_", site_id), 19, 15)
    cli::cli_alert_success("{.path {pth[['png']]}}")
  }
  o
}

wide <- list()
wide_em <- list()
nulls <- list()
reg_h <- list()
periods <- list()
site_status <- list()
for (site_id in site_ids) {
  cli::cli_h1("Gamma cases: {site_id}")
  o <- new.env()
  o$stage <- "metrics"
  res <- tryCatch(eval_site(site_id, o), error = function(e) e)
  if (inherits(res, "error")) {
    msg <- gsub("\\s+", " ", paste(conditionMessage(res), collapse = " "))
    if (identical(o$stage, "growth")) {
      # the metrics of the site are complete: keep them, report the growth-curve error
      cli::cli_alert_danger("{site_id}: error in the clear-day growth curves (metrics kept): {msg}")
      o$status <- "error_growth"
      o$reason <- paste("growth curves:", msg)
    } else {
      # partial metrics would be inconsistent between cases: the site gets no metric rows, only its status
      cli::cli_alert_danger("{site_id}: error, site skipped (no metric rows): {msg}")
      o$wide <- o$wide_em <- o$nulls <- o$reg_h <- list()
      o$period <- NULL
      o$status <- "error"
      o$reason <- msg
    }
  }
  wide <- c(wide, o$wide)
  wide_em <- c(wide_em, o$wide_em)
  nulls <- c(nulls, o$nulls)
  reg_h <- c(reg_h, o$reg_h)
  if (!is.null(o$period)) periods[[site_id]] <- o$period
  site_status[[site_id]] <- dplyr::tibble(site_id = site_id, status = o$status, no_data_reason = o$reason)
}
site_status <- dplyr::bind_rows(site_status)
if (all(site_status$status == "error")) cli::cli_abort("Every site failed: nothing to write.")

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
wide_em <- dplyr::bind_rows(wide_em)
wide_all <- dplyr::bind_rows(wide, wide_em)
# Season table (amendment 1): every case and ERA5 against the measurement, and every case against ERA5
season_tab <- gamma_season_table(wide_all, min_days)
pairing <- function(s, what) {
  vapply(s, function(x) if (is.null(refs[[x]])) NA_character_ else as.character(refs[[x]][[what]]), "")
}
season_tab$distance_class <- unname(pairing(season_tab$site_id, "distance_class"))
season_tab$borderline <- as.logical(unname(pairing(season_tab$site_id, "borderline")))
season_tab$note <- paste("ln(source / reference); case era5_blh_center = ERA5 tower cell vs the measurement;",
                         "low_n: fewer than", min_days, "days in the day-block bootstrap")
readr::write_csv(season_tab, file.path(tables_dir, "11_gamma_cases_season.csv"))
# Decision table (amendment 1): one row per site
dec_cases <- c(cases, "era5_blh_center")
decision <- gamma_decision_table(wide_all, site_ids, refs, dec_cases, min_days)
per <- dplyr::bind_rows(c(list(dplyr::tibble(site_id = character(), period_start_utc = as.POSIXct(character(), tz = "UTC"),
                                             period_end_utc = as.POSIXct(character(), tz = "UTC"))), periods))
decision <- dplyr::left_join(decision, per, by = "site_id")
decision <- dplyr::relocate(decision, period_start_utc, period_end_utc, .after = borderline)
decision <- dplyr::left_join(decision, site_status, by = "site_id")
decision <- dplyr::relocate(decision, status, no_data_reason, .after = period_end_utc)
decision$note <- paste("mlr = median ln(PBLHmod case / measured) on the primary set (aerosol: window",
                       paste(opts$aerosol_window, collapse = "-"), "h apparent solar time and PBLHmod convective;",
                       "turbulence: convective, exact values), day-block bootstrap 95 %; era5_blh_center = ERA5 tower",
                       "cell vs measured; season_max = season with the largest |mlr| among seasons with >=", min_days,
                       "days; era5_mlr = median ln(PBLHmod case / ERA5 tower cell), convective; low_n: fewer than",
                       min_days, "days; status: ok, era5_only, no_overlap (n = 0, reason in no_data_reason), error",
                       "(no metric rows), error_growth (metrics kept)")
readr::write_csv(decision, file.path(tables_dir, "11_gamma_decision.csv"))
reg_h <- dplyr::bind_rows(reg_h)
readr::write_csv(reg_h, file.path(tables_dir, "11_gamma_regime_heights.csv"))
if (length(nulls) > 0) {
  nulls <- dplyr::bind_rows(nulls)
  readr::write_csv(nulls, file.path(tables_dir, "11_gamma_censored_null.csv"))
}

# --- gamma diagnostics: daily vs climatology, by month ---
gd <- dplyr::bind_rows(lapply(site_ids, function(s) {
  f <- file.path(g_dir, paste0(s, c("_daily.csv", "_monthly.csv")))
  if (!all(file.exists(f))) {
    cli::cli_alert_warning("{s}: gamma table{?s} {.path {basename(f[!file.exists(f)])}} missing (script 11a): no gamma diagnostics.")
    return(NULL)
  }
  d <- readr::read_csv(f[1], show_col_types = FALSE)
  cl <- readr::read_csv(f[2], show_col_types = FALSE)
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
  ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3, fill = style_of("gamma_era5", "colour", "element"), alpha = 0.4) +
  ggplot2::facet_wrap(~site_id, nrow = 2) +
  ggplot2::labs(x = "Month", y = "Daily morning gamma [K m-1] (ERA5, 500-3500 m above the site)",
                subtitle = "Lines: 0.01 (PBLHmod default, solid) and 0.0025 K m-1 (dashed)") +
  theme_pblh()
pth <- save_figure(fig_g, paste0("11_gamma_by_month_", paste(site_ids, collapse = "_")), 19, 13)
cli::cli_alert_success("{.path {pth[['png']]}}")

# --- case figures: median log ratio with bootstrap interval ---
for (s in site_ids) {
  for (fs in fig_sets) {
    w <- wide[wide$site_id == s & wide$case %in% fs$cases & !grepl(":", wide$stratum_set, fixed = TRUE), ]
    if (nrow(w) == 0) next # site with an error: no rows
    w$panel <- paste0(ifelse(w$reference == "era5_blh_center", "vs ERA5 tower cell", paste("vs", w$reference)), ", ",
                      w$stratum_set)
    w$case <- factor(w$case, levels = fs$cases)
    fig_c <- ggplot2::ggplot(w, ggplot2::aes(case, median_log_ratio, colour = case)) +
      ratio_lines() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mlr_lo, ymax = mlr_hi), width = 0.2, linewidth = 0.4) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::facet_wrap(~panel) +
      ggplot2::scale_x_discrete(labels = case_lab) +
      ggplot2::scale_colour_manual(values = case_pal, guide = "none") +
      ggplot2::labs(x = NULL, y = "Median ln(PBLHmod / reference) [-] (day-block bootstrap 95 %)", title = fs$title,
                    subtitle = s) +
      theme_pblh()
    # a figure of one site never stops the job (the tables are already written)
    pth <- tryCatch(save_figure(fig_c, paste0(fs$prefix, "gamma_cases_", s), 19, 13), error = function(e) {
      cli::cli_alert_danger("{s}: case figure {fs$prefix} not written: {conditionMessage(e)}")
      NULL
    })
    if (!is.null(pth)) cli::cli_alert_success("{.path {pth[['png']]}}")
  }
}

# --- decision overview (amendment 1): sites on the y axis, median ln ratio with interval per case ---
dm <- wide_all[wide_all$reference != "era5_blh_center" & wide_all$case %in% c(intersect(styled, cases), "era5_blh_center"), ]
dm <- dm[dm$stratum_set == vapply(dm$site_id, function(s) refs[[s]]$primary_set, ""), ]
if (nrow(dm) > 0) {
  ov_ids <- c(paste0("pblhmod_", intersect(styled, cases)), "era5_blh_center")
  ov_ids <- style_ids("source", ov_ids)
  dm$source <- ifelse(dm$case == "era5_blh_center", "era5_blh_center", paste0("pblhmod_", dm$case))
  dm$source <- factor(dm$source, levels = rev(ov_ids)) # first source on top after the flip
  ms_sites <- intersect(site_ids, unique(dm$site_id))
  lab_site <- vapply(ms_sites, function(s) {
    rf <- refs[[s]]
    nd <- decision$n_days_default[decision$site_id == s]
    st <- site_status$status[site_status$site_id == s]
    paste0(s, " (", unname(style_of(rf$src, "label")), ")", if (rf$borderline) " [borderline]" else "",
           if (identical(st, "no_overlap")) " [no overlap]" else if (isTRUE(nd < min_days)) paste0(" [low n: ", nd, " d]") else "")
  }, "")
  dm$site <- factor(lab_site[dm$site_id], levels = rev(lab_site))
  dm$low_n <- ifelse(is.na(dm$n_days) | dm$n_days < min_days, "low n", "ok")
  pos <- ggplot2::position_dodge(width = 0.75)
  fig_d <- ggplot2::ggplot(dm, ggplot2::aes(site, median_log_ratio, colour = source, group = source)) +
    ratio_lines() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mlr_lo, ymax = mlr_hi), width = 0, linewidth = 0.4, position = pos,
                           na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(shape = low_n), size = 1.8, position = pos, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = style_of(ov_ids), labels = style_of(ov_ids, "label"), breaks = ov_ids,
                                 name = NULL) +
    ggplot2::scale_shape_manual(values = c(ok = 16, `low n` = 1), name = NULL,
                                labels = c(ok = paste0(">= ", min_days, " days"), `low n` = paste0("< ", min_days, " days"))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Median ln(source / measured) [-] (day-block bootstrap 95 %)",
                  subtitle = sprintf(paste("Aerosol: %g-%g h apparent solar time, PBLHmod convective;\n",
                                           "Cloudnet (FI-Hyy): PBLHmod convective, exact values"), opts$aerosol_window[1],
                                     opts$aerosol_window[2])) +
    theme_pblh() +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
  pth <- save_figure(fig_d, "11_gamma_decision", 19, 4 + 1.6 * length(ms_sites))
  cli::cli_alert_success("{.path {pth[['png']]}}")
}

# --- report ---
cli::cli_h1("Summary")
cli::cli_h2("Decision table: median ln ratio on the primary set [95 %], n days")
dec_print <- decision[c("site_id", "stratum_set", "borderline")]
for (k in dec_cases) {
  dec_print[[k]] <- ifelse(is.na(decision[[paste0("mlr_", k)]]), NA_character_,
                           sprintf("%+.2f [%+.2f, %+.2f] %s d%s", decision[[paste0("mlr_", k)]],
                                   decision[[paste0("mlr_lo_", k)]], decision[[paste0("mlr_hi_", k)]],
                                   decision[[paste0("n_days_", k)]],
                                   ifelse(decision[[paste0("low_n_", k)]] %in% TRUE, " low_n", "")))
}
print(as.data.frame(dec_print), row.names = FALSE)
cli::cli_h2("Primary sets: median ln ratio [95 %], within factor 2, FB, nRMSD, symmetric slope")
prim <- wide[(wide$reference != "era5_blh_center" & wide$stratum_set %in% c("window_convective", "convective")) |
               (wide$reference == "era5_blh_center" & !grepl(":", wide$stratum_set, fixed = TRUE)), ]
print_tbl(prim[c("site_id", "reference", "stratum_set", "case", "n", "median_log_ratio", "mlr_lo", "mlr_hi",
                 "frac_within_factor2", "fb", "nrmsd", "slope_sym")])
cli::cli_h2("Regime shares and heights (default, clim, gamma0025, ERA5 clim in B&G only)")
print_tbl(reg_h)
if (is.data.frame(nulls)) {
  cli::cli_h2("FI-Hyy censored-tail test by case")
  print_tbl(nulls[c("case", "regime_set", "tail", "n", "obs_share", "obs_lo", "obs_hi", "sim_share", "sim_lo", "sim_hi",
                    "sim_bound")])
}
cli::cli_h2("gamma: daily vs climatology")
print_tbl(dvc)
cli::cli_h2("Site status (also in 11_gamma_decision.csv)")
print(as.data.frame(site_status), row.names = FALSE)
n_err <- sum(startsWith(site_status$status, "error"))
if (n_err > 0) cli::cli_alert_danger("{n_err} site{?s} with an error: see the status table above and the log.")
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
# Exit status (G 2026-11-11): 1 if any site has status error / error_growth (all tables and figures are written
# first), so that sacct shows FAILED; 0 otherwise (no_overlap and era5_only are not errors).
if (n_err > 0) quit(save = "no", status = 1)
