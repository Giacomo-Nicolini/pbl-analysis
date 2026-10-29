# Spec 10, PBLHmod deficit against the aerosol MLH (G 2026-10-29, item 3): diagnostics from existing outputs only,
# no model change (gamma from ERA5 is the next spec).
#   (a) Window metrics (pblh_metrics(), PBLHmod convective, apparent solar time 10-16 h at the interval midpoint)
#       of every PBLHmod case already run (<derived_dir>/pblhmod/<site>_pblhmod_<label>.csv.gz) against CHM15k,
#       on each case's own window hours and on the window hours of the default case.
#   (b) Batchvarova-Gryning start height h0 at sunrise (bg_start_heights()): distribution by season, ERA5 BLH at
#       the same time (context), and the window log ratio of the day against h0.
#   (c) Window hours on days where B&G stopped at a forcing gap (flag_bg_stopped_gap).
#   (d) Growth curves on clear days: median by apparent solar hour of PBLHmod (any regime), CHM15k and ERA5 tower
#       cell on the hours where all three are valid. Clear day: daily clearness index kt = sum SW_IN / sum SW_IN_POT
#       >= 0.65 over the daylight half-hours with SW_IN present on >= 90 % of them (choice of this script, flagged).
#   Matching as script 10 (measured_time_base(): hourly intervals [t, t + 1 h)); comparison period = overlap of
#   valid PBLHmod default, CHM15k and CL31, as script 10.
#
# Output (<output_dir>): tables/10b_cases_window_<site>.csv, 10b_h0_<site>.csv (per day), 10b_h0_summary_<site>.csv,
# 10b_stopped_days_<site>.csv, 10b_growth_clear_days_<site>.csv; figures/10b_growth_clear_days_<site>.png.
#
# Usage: Rscript scripts/10b_pblhmod_deficit.R SITE   (SITE with a SIRTA-ReOBS pairing, i.e. FR-Gri)
#   Profile: env var R_CONFIG_ACTIVE (default | cluster).

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

site_id <- commandArgs(trailingOnly = TRUE)
if (length(site_id) != 1) cli::cli_abort("Pass exactly one ICOS station code, e.g. {.val FR-Gri}.")
opts <- measured_compare_options()
kt_clear <- 0.65 # clearness threshold (this script)
kt_cov <- 0.9    # minimum SW_IN coverage of the daylight half-hours
cfg <- get_cfg()
sites <- read_sites(resolve_path(cfg$sites_file))
meta_site <- sites[sites$site_id == site_id, ]
std_dir <- file.path(resolve_path(cfg$standard_dir), site_id)
pbl_dir <- file.path(resolve_path(cfg$derived_dir), "pblhmod")
out_root <- resolve_path(cfg$output_dir, must_exist = FALSE)
tables_dir <- file.path(out_root, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
print_tbl <- function(x) print(as.data.frame(dplyr::mutate(x, dplyr::across(dplyr::where(is.double), ~ signif(.x, 4)))),
                               row.names = FALSE)
rd <- function(src) read_standard_pblh(file.path(std_dir, paste0(site_id, "__", src, ".csv.gz")))
wr <- function(x, name) readr::write_csv(x, file.path(tables_dir, paste0(name, "_", site_id, ".csv")))
t0 <- Sys.time()
cli::cli_h1("PBLHmod deficit diagnostics (spec 10): {site_id}")

chm <- rd("sirta_mlh_chm15k")
cl <- rd("sirta_mlh_cl31")
mod_def <- rd("pblhmod_default")
era <- list(era5_blh_center = rd("era5_blh_center"))
v_rng <- function(d) range(d$timestamp_start_utc[!is.na(d$pblh_m)])
p_start <- max(v_rng(mod_def)[1], v_rng(chm)[1], v_rng(cl)[1])
p_end <- min(v_rng(mod_def)[2], v_rng(chm)[2], v_rng(cl)[2])
chm <- chm[chm$timestamp_start_utc >= p_start & chm$timestamp_start_utc <= p_end, ]
cli::cli_inform("Comparison period: {format(p_start, '%Y-%m-%d %H:%M')} to {format(p_end, '%Y-%m-%d %H:%M')} UTC (as script 10).")
win_of <- function(t_start) {
  a <- apparent_solar_hour(t_start + 1800, meta_site$lon)
  a >= opts$aerosol_window[1] & a < opts$aerosol_window[2]
}
exact_ref <- function(v) v[v$kind %in% "exact", c("timestamp_utc", "pblh_m")]
mets <- c("n", "bias", "median_log_ratio", "mad_log_ratio", "frac_within_factor2", "r_spearman", "fb", "nrmsd")

# --- (a) all PBLHmod cases already run ---
files <- list.files(pbl_dir, pattern = paste0("^", site_id, "_pblhmod_.*\\.csv\\.gz$"), full.names = TRUE)
labels <- sub(paste0("^", site_id, "_pblhmod_(.*)\\.csv\\.gz$"), "\\1", basename(files))
cli::cli_h2("(a) Window metrics of {length(files)} PBLHmod case{?s} vs CHM15k")
raw <- list()
pairs <- list()
for (k in seq_along(files)) {
  d <- readr::read_csv(files[k], show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
  if (labels[k] == "default") raw$default <- d
  std <- suppressMessages(standardize_pblhmod(d, labels[k])$data)
  tb <- measured_time_base(chm, std, era)
  pt <- pair_table(tb$mod, exact_ref(tb$meas), labels[k])
  pt$in_window <- win_of(pt$timestamp_utc)
  pairs[[labels[k]]] <- pt
}
def_hours <- pairs$default$timestamp_utc[pairs$default$regime %in% "convective" & pairs$default$in_window]
cases <- dplyr::bind_rows(lapply(names(pairs), function(k) {
  pt <- pairs[[k]]
  own <- pt[pt$regime %in% "convective" & pt$in_window, ]
  dh <- pt[pt$timestamp_utc %in% def_hours, ]
  dplyr::bind_rows(dplyr::bind_cols(dplyr::tibble(case = k, hours = "own window convective"), pblh_metrics(own$mod, own$ref)),
                   dplyr::bind_cols(dplyr::tibble(case = k, hours = "default window convective"), pblh_metrics(dh$mod, dh$ref)))
}))
cases$reference <- "sirta_mlh_chm15k"
wr(cases, "10b_cases_window")
print_tbl(cases[c("case", "hours", mets)])

# --- (b) B&G start height ---
cli::cli_h2("(b) B&G start height h0 at sunrise (default case)")
d <- raw$default
h0 <- bg_start_heights(d)
h0 <- h0[h0$sunrise_utc >= p_start & h0$sunrise_utc <= p_end, ]
h0$season <- season_of(h0$sunrise_utc)
h0$bg_run <- h0$n_convective > 0
# reconstruction check: where wtheta <= 0 at the first daylight half-hour, the first window value is h0 (held)
w1 <- d$wtheta[match(as.numeric(h0$sunrise_utc), as.numeric(as_utc_time(d$timestamp_start_utc)))]
chk <- h0$bg_run & !is.na(w1) & w1 <= 0 & !is.na(h0$h0_m)
cli::cli_inform("Reconstruction check: first window value = h0 on {sum(abs(h0$first_h_bg_m[chk] - h0$h0_m[chk]) < 1e-6, na.rm = TRUE)} of {sum(chk)} days with wtheta <= 0 at sunrise.")
# ERA5 BLH at the time of h0 (end of the h0 half-hour = sunrise minus (lag - 1) half-hours, nearest hour)
t_h0 <- h0$sunrise_utc - (h0$h0_lag_halfhours - 1) * 1800
hr <- as.POSIXct(round(as.numeric(t_h0) / 3600) * 3600, origin = "1970-01-01", tz = "UTC")
h0$era5_blh_m <- era$era5_blh_center$pblh_m[match(as.numeric(hr), as.numeric(era$era5_blh_center$timestamp_start_utc))]
# window log ratio of the day (median over window convective hours of the default case)
pd <- pairs$default[pairs$default$regime %in% "convective" & pairs$default$in_window & pairs$default$mod > 0, ]
pd$day_id <- as.Date(d$day_id)[match(as.numeric(pd$timestamp_utc), as.numeric(as_utc_time(d$timestamp_start_utc)))]
day_lr <- dplyr::summarise(dplyr::group_by(pd, day_id), n_win = dplyr::n(), mlr_day = stats::median(log(mod / ref)),
                           .groups = "drop")
h0 <- dplyr::left_join(h0, day_lr, by = "day_id")
wr(h0, "10b_h0")
qs <- function(x) stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE, names = FALSE)
h0s <- dplyr::bind_rows(lapply(c("all", levels(h0$season)), function(s) {
  x <- h0[h0$bg_run & (s == "all" | h0$season == s), ]
  q <- qs(x$h0_m)
  e <- qs(x$era5_blh_m)
  dplyr::tibble(season = s, n_days = nrow(x), h0_p05 = q[1], h0_p25 = q[2], h0_p50 = q[3], h0_p75 = q[4],
                h0_p95 = q[5], frac_h0_B1 = mean(x$h0_method %in% "B1"), frac_h0_lt_50m = mean(x$h0_m < 50, na.rm = TRUE),
                era5_p50 = e[3], era5_p25 = e[2], era5_p75 = e[4],
                rho_h0_vs_day_logratio = suppressWarnings(stats::cor(x$h0_m, x$mlr_day, method = "spearman",
                                                                     use = "complete.obs")))
}))
wr(h0s, "10b_h0_summary")
cli::cli_inform("Days in the period: {nrow(h0)}; with a B&G window: {sum(h0$bg_run)}; no h0: {sum(h0$flag_no_h0)}.")
print_tbl(h0s)
x <- h0[h0$bg_run & !is.na(h0$mlr_day), ]
if (nrow(x) >= 30) {
  x$h0_tercile <- cut(x$h0_m, stats::quantile(x$h0_m, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE), include.lowest = TRUE)
  cli::cli_inform("Daily median window ln(PBLHmod / CHM15k) by tercile of h0:")
  print_tbl(dplyr::summarise(dplyr::group_by(x, h0_tercile), n_days = dplyr::n(), mlr_day_median = stats::median(mlr_day),
                             .groups = "drop"))
}

# --- (c) days where B&G stopped at a gap ---
cli::cli_h2("(c) Window hours on days where B&G stopped at a forcing gap")
pa <- pairs$default
pa$day_id <- as.Date(d$day_id)[match(as.numeric(pa$timestamp_utc), as.numeric(as_utc_time(d$timestamp_start_utc)))]
pa$stopped <- h0$flag_bg_stopped_gap[match(pa$day_id, h0$day_id)] %in% TRUE
# all window hours with a CHM15k value (PBLHmod valid, any regime)
ww <- pa[pa$in_window, ]
stp <- dplyr::bind_rows(
  dplyr::tibble(item = "days with a B&G window", value = sum(h0$bg_run)),
  dplyr::tibble(item = "of which stopped at a gap", value = sum(h0$bg_run & h0$flag_bg_stopped_gap)),
  dplyr::tibble(item = "days stopped at a gap (all days in the period)", value = sum(h0$flag_bg_stopped_gap)),
  dplyr::tibble(item = "window hours with CHM15k and PBLHmod", value = nrow(ww)),
  dplyr::tibble(item = "  on stopped days", value = sum(ww$stopped)),
  dplyr::tibble(item = "  on stopped days, convective (before the stop)", value = sum(ww$stopped & ww$regime %in% "convective")),
  dplyr::tibble(item = "  on stopped days, not convective (after the stop: diagnostic)", value = sum(ww$stopped & !ww$regime %in% "convective")),
  dplyr::tibble(item = "  on other days, not convective", value = sum(!ww$stopped & !ww$regime %in% "convective")),
  dplyr::tibble(item = "window convective hours (primary set)", value = sum(ww$regime %in% "convective")),
  dplyr::tibble(item = "  share on stopped days", value = mean(ww$stopped[ww$regime %in% "convective"]))
)
wc <- ww[ww$regime %in% "convective", ]
stm <- dplyr::bind_rows(
  dplyr::bind_cols(dplyr::tibble(days = "stopped"), pblh_metrics(wc$mod[wc$stopped], wc$ref[wc$stopped])),
  dplyr::bind_cols(dplyr::tibble(days = "not stopped"), pblh_metrics(wc$mod[!wc$stopped], wc$ref[!wc$stopped]))
)
stm_long <- dplyr::bind_rows(lapply(mets, function(m) {
  dplyr::tibble(item = paste0(m, " (", stm$days, " days)"), value = as.numeric(stm[[m]]))
}))
wr(dplyr::bind_rows(dplyr::mutate(stp, table = "counts"), dplyr::mutate(stm_long, table = "metrics, window convective")),
   "10b_stopped_days")
print_tbl(stp)
print_tbl(stm[c("days", mets)])

# --- (d) growth curves on clear days ---
cli::cli_h2("(d) Growth curves on clear days")
sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
exclusions <- read_exclusions(resolve_path(cfg$exclusions_file))
icos <- suppressWarnings(read_icos_site(site_id, cfg, sites, sites_tv, exclusions))
day_of <- as.Date(d$day_id)[match(as.numeric(icos$timestamp_start_utc), as.numeric(as_utc_time(d$timestamp_start_utc)))]
kt <- daily_clearness(day_of, icos$met_SW_IN, icos$fn_SW_IN_POT)
clear <- kt$day_id[!is.na(kt$kt) & kt$coverage >= kt_cov & kt$kt >= kt_clear]
cli::cli_inform("Clear days (kt >= {kt_clear}, SW_IN coverage >= {kt_cov}): {length(clear)} of {sum(kt$coverage >= kt_cov)} days with radiation.")
tb <- measured_time_base(chm, mod_def, era)
g <- dplyr::tibble(timestamp_utc = tb$meas$timestamp_utc, sirta_mlh_chm15k = ifelse(tb$meas$kind %in% "exact", tb$meas$pblh_m, NA))
g$pblhmod_default <- tb$mod$pblh_m[match(as.numeric(g$timestamp_utc), as.numeric(tb$mod$timestamp_utc))]
g$era5_blh_center <- tb$era$era5_blh_center$pblh_m[match(as.numeric(g$timestamp_utc), as.numeric(tb$era$era5_blh_center$timestamp_utc))]
g$day_id <- as.Date(d$day_id)[match(as.numeric(g$timestamp_utc), as.numeric(as_utc_time(d$timestamp_start_utc)))]
g <- g[g$day_id %in% clear & stats::complete.cases(g[c("sirta_mlh_chm15k", "pblhmod_default", "era5_blh_center")]), ]
g$season <- season_of(g$timestamp_utc)
g$ast_hour <- floor(apparent_solar_hour(g$timestamp_utc + 1800, meta_site$lon))
vars <- c("sirta_mlh_chm15k", "pblhmod_default", "era5_blh_center")
gl <- dplyr::bind_rows(lapply(vars, function(v) {
  dplyr::tibble(season = g$season, ast_hour = g$ast_hour, variable = v, pblh_m = g[[v]])
}))
gs <- dplyr::summarise(dplyr::group_by(gl, season, ast_hour, variable), n = dplyr::n(),
                       median = stats::median(pblh_m), q25 = stats::quantile(pblh_m, 0.25, names = FALSE),
                       q75 = stats::quantile(pblh_m, 0.75, names = FALSE), .groups = "drop")
nd <- dplyr::summarise(dplyr::group_by(g, season), n_clear_days = dplyr::n_distinct(day_id), .groups = "drop")
gs <- dplyr::left_join(gs, nd, by = "season")
wr(gs, "10b_growth_clear_days")
pal <- aerosol_palette()
gp <- gs[gs$n >= 10, ]
gp$x <- gp$ast_hour + 0.5
gp$variable <- factor(gp$variable, levels = intersect(names(pal), unique(gp$variable)))
nd$label <- paste0(nd$n_clear_days, " clear days")
fig <- ggplot2::ggplot(gp, ggplot2::aes(x, median, colour = variable, fill = variable)) +
  ggplot2::annotate("rect", xmin = opts$aerosol_window[1], xmax = opts$aerosol_window[2], ymin = -Inf, ymax = Inf,
                    fill = "#F0E442", alpha = 0.25) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_text(data = nd, ggplot2::aes(x = 5.2, y = Inf, label = label), inherit.aes = FALSE, hjust = 0,
                     vjust = 1.3, size = 2.6) +
  ggplot2::facet_wrap(~season, nrow = 2) +
  ggplot2::scale_colour_manual(values = pal, labels = aerosol_labels(), name = NULL) +
  ggplot2::scale_fill_manual(values = pal, labels = aerosol_labels(), name = NULL) +
  ggplot2::scale_x_continuous(breaks = seq(4, 22, 2), limits = c(5, 21)) +
  ggplot2::labs(x = "Apparent solar time [h, interval midpoint]", y = "Median PBLH / MLH [m] (band: IQR)",
                subtitle = sprintf("Clear days (kt >= %.2f); hours with all three valid; PBLHmod any regime; n >= 10 per hour",
                                   kt_clear)) +
  theme_pub()
pth <- save_figure(fig, paste0("10b_growth_clear_days_", site_id), 17, 13)
cli::cli_alert_success("{.path {pth[['png']]}}")
for (s in c("MAM", "JJA")) {
  cat("--", s, "\n")
  x <- gs[gs$season == s & gs$ast_hour >= 5 & gs$ast_hour <= 20, ]
  w <- dplyr::tibble(ast_hour = sort(unique(x$ast_hour)))
  for (v in vars) {
    xv <- x[x$variable == v, ]
    w[[v]] <- round(xv$median[match(w$ast_hour, xv$ast_hour)])
  }
  w$n <- x$n[x$variable == "sirta_mlh_chm15k"][match(w$ast_hour, x$ast_hour[x$variable == "sirta_mlh_chm15k"])]
  print_tbl(w)
}
cli::cli_inform("Run time: {format(round(difftime(Sys.time(), t0, units = 'mins'), 1))}.")
