# Spec 17 sections 3 C, D, E: network collection of the EC-derived gamma (gamma_EC). Reads only per-site outputs:
# <derived_dir>/gamma_ec/<site>_{daily,monthly,ebr}.csv (17a), <derived_dir>/gamma/<site>_{daily,monthly}.csv (11a),
# <output_dir>/tables/04_metrics_<site>.csv (pairs 1, 6, 9), and for part E the footprint rows of the reference year
# (fp_site_rows(): 08a input table and 08b chunks) with the standard file pblhmod_gamma_ec (step 04).
#   C: monthly gamma_EC vs the ERA5 climatology, share of months with a value, used days per site; PBLHmod-gammaEC vs
#      ERA5 tower cell (pair 9) next to pairs 1 and 6.
#   D: diagnostic only, no correction: gamma_EC / gamma_ERA5 (median of the monthly ratios) against the daytime JJA
#      energy-balance ratio (H + LE) / (NETRAD - G) of the reference year (17a); availability per site.
#   E: bracketing check for decision D1: share of the footprint rows (common_pe, etc, valid_qc_used) where h of
#      PBLHmod-gammaEC is outside [min, max] of (ERA5, PBLHmod, PBLHmod-clim), median ln distance outside, per site x
#      zm/L class x day/night. No FFP run.
#
# Outputs (<output_dir>/tables, `#` header; <output_dir>/figures): 17_gamma_ec_vs_era5.csv (site x month),
#   17_gamma_ec_sites.csv (one row per site: C and D), 17_network_logratio.csv (pairs 1, 6, 9), 17_gamma_ec_bracketing.csv
#   (E); figures 17_gamma_ec_vs_era5, 17_gamma_ec_vs_era5_seasonal, 17s_gamma_ec_ebr.
# Supplement (G 2026-11-26, D1 rejected): --offset-h=H reads the 17a tables of fit-start offset H (gamma_ec_off<H>h/)
#   and writes the C/D tables and figures with the suffix _off<H>h (no 17_network_logratio). With offset 0 and part D,
#   also 17_gamma_ratio_vs_height.csv (site x offset: site gamma_EC / gamma_ERA5 with zm = EC height - 2/3 hc, EC
#   height, canopy height on the used days, IGBP; every offset folder found) and figure 17s_gamma_ratio_vs_height.
#
# Usage: Rscript scripts/17b_gamma_ec_collect.R [--sites-file=config/run_sites.txt | --sites=A,B] [--parts=C,D,E]
#          [--label=test] [--offset-h=H]

for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, def = NULL) {
  v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
  if (length(v)) v else def
}
bad <- setdiff(sub("=.*$", "", grep("^--", args, value = TRUE)), c("--sites-file", "--sites", "--parts", "--label", "--offset-h"))
if (length(bad)) cli::cli_abort("Unknown argument{?s}: {.val {bad}}.")
split <- function(x) trimws(strsplit(x, ",")[[1]])
parts <- split(opt("parts", "C,D,E"))
label <- opt("label", "")
off_h <- as.numeric(opt("offset-h", "0"))
ec_name <- gamma_ec_dir_name(off_h)
sfx <- if (ec_name == "gamma_ec") "" else paste0("_", sub("^gamma_ec_", "", ec_name))
if (nzchar(sfx)) label <- paste0(c(if (nzchar(label)) label, paste0("fit start +", format(off_h), " h")), collapse = ", ")
cfg <- get_cfg()
ids <- if (!is.null(opt("sites"))) split(opt("sites")) else {
  x <- trimws(readLines(opt("sites-file", here::here("config", "run_sites.txt")), warn = FALSE))
  x[nzchar(x) & !startsWith(x, "#")]
}
derived_dir <- resolve_path(cfg$derived_dir, must_exist = FALSE)
out_dir <- resolve_path(cfg$output_dir, must_exist = FALSE)
tables_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(p) if (file.exists(p)) readr::read_csv(p, comment = "#", show_col_types = FALSE, progress = FALSE) else NULL
print_tbl <- function(x) print(as.data.frame(dplyr::mutate(x, dplyr::across(dplyr::where(is.double), ~ signif(.x, 3)))),
                               row.names = FALSE)
tag <- if (nzchar(label)) paste0(" (", label, ")") else ""

if (any(c("C", "D") %in% parts)) {
  cli::cli_h1("Parts C, D: gamma_EC vs ERA5, energy-balance ratio{tag}")
  mon <- list(); summ <- list(); lr <- list(); miss <- character()
  for (s in ids) {
    ecd <- rd(file.path(derived_dir, ec_name, paste0(s, "_daily.csv")))
    ecm <- rd(file.path(derived_dir, ec_name, paste0(s, "_monthly.csv")))
    if (is.null(ecd) || is.null(ecm)) {
      miss <- c(miss, s)
      next
    }
    m <- gamma_ec_vs_era5_table(ecm, rd(file.path(derived_dir, "gamma", paste0(s, "_monthly.csv"))), s)
    mon[[s]] <- m
    summ[[s]] <- gamma_ec_site_summary(m, ecd, rd(file.path(derived_dir, "gamma", paste0(s, "_daily.csv"))),
                                       rd(file.path(derived_dir, ec_name, paste0(s, "_ebr.csv"))))
    met <- rd(file.path(tables_dir, paste0("04_metrics_", s, ".csv")))
    if (!is.null(met)) {
      x <- met[met$pair %in% c("1_pblhmod_default_vs_era5_center", "6_pblhmod_clim_vs_era5_center",
                               "9_pblhmod_gamma_ec_vs_era5_center") &
                 met$stratum %in% c("all", "regime=convective", "regime=stable") &
                 met$metric %in% c("n", "median_log_ratio", "mad_log_ratio"), c("pair", "stratum", "metric", "value")]
      lr[[s]] <- dplyr::mutate(x, site_id = s, .before = 1)
    }
  }
  if (length(miss)) cli::cli_alert_warning("No gamma_EC tables (17a) for {length(miss)} site{?s}: {paste(miss, collapse = ', ')}.")
  mon <- dplyr::bind_rows(mon)
  summ <- dplyr::bind_rows(summ)
  lr <- dplyr::bind_rows(lr)
  hd <- c(paste0("Spec 17 C/D (scripts/17b_gamma_ec_collect.R)", tag, ": gamma_EC = Method A of scripts/17a_gamma_ec.R on the used days; monthly median, p25, p75 from >= gamma_ec_min_days used days [K m-1]; ERA5 = monthly median of the daily morning gamma of 11a (500-3500 m above the site) [K m-1]."),
          "ratio_ec_era5 = gamma_EC / gamma_ERA5 [-]; ebr = daytime (SW_IN_POT > 0) JJA sum(H + LE) / sum(NETRAD - G) of the reference year with H, LE measured (QC 0); ebr_no_g = (H + LE) / NETRAD (diagnostic only, no correction).")
  write_csv_header(mon, file.path(tables_dir, paste0("17_gamma_ec_vs_era5", sfx, ".csv")), hd)
  write_csv_header(summ, file.path(tables_dir, paste0("17_gamma_ec_sites", sfx, ".csv")), hd)
  if (nrow(lr) && !nzchar(sfx)) write_csv_header(lr, file.path(tables_dir, "17_network_logratio.csv"),
                                 c(paste0("Spec 17 C", tag, ": pairs 1 (PBLHmod), 6 (PBLHmod-clim), 9 (PBLHmod-gammaEC) vs ERA5 tower cell from 04_metrics_<site>.csv; median_log_ratio = median ln(source / ERA5) [-]."),
                                   "Pair 9 is absent at a site where step 04 ran before the gamma_ec case of step 02."))
  cap <- fp_wrap(paste0("Monthly medians of the used days (γEC, Method A) and of the daily ERA5 morning γ; months with fewer than ",
                        gamma_ec_options_cfg(cfg)$min_days, " used γEC days are not shown. Dashed line: 0.01 K m-1 (PBLHmod).", tag))
  if (nrow(mon) && any(is.finite(mon$gamma_ec))) {
    save_figure(plot_gamma_ec_vs_era5(mon, cap), paste0("17_gamma_ec_vs_era5", sfx), 12, 11, dir = fig_dir)
    nsite <- length(unique(mon$site_id))
    save_figure(plot_gamma_ec_seasonal(mon, cap), paste0("17_gamma_ec_vs_era5_seasonal", sfx), 19, 6 + 3.2 * ceiling(nsite / 6), dir = fig_dir)
  }
  if ("D" %in% parts && nrow(summ) && any(is.finite(summ$ebr) | is.finite(summ$ebr_no_g))) {
    d <- summ[is.finite(summ$ratio_ec_era5), ]
    d$ebr_used <- ifelse(is.finite(d$ebr), d$ebr, d$ebr_no_g)
    d$g <- ifelse(is.finite(d$ebr), "with G", "without G (no G at the site)")
    fig <- ggplot2::ggplot(d, ggplot2::aes(ebr_used, ratio_ec_era5, shape = g)) +
      ggplot2::geom_hline(yintercept = 1, colour = style_of("reference_line", "colour", "element")) +
      ggplot2::geom_point(colour = style_of("pblhmod_gamma_ec"), size = 2) +
      ggplot2::geom_text(ggplot2::aes(label = site_id), size = 2.4, vjust = -0.9, show.legend = FALSE) +
      ggplot2::scale_shape_manual(values = c(`with G` = 16, `without G (no G at the site)` = 1), name = NULL) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Daytime JJA energy-balance ratio of the reference year [-]",
                    y = "γEC / γERA5 (median of the monthly ratios) [-]",
                    caption = fp_wrap(paste0("Diagnostic only (no H correction, G 2026-11-12): a γ fitted with the H that drives the model absorbs a systematic error of H.", tag))) +
      theme_pblh()
    save_figure(fig, paste0("17s_gamma_ec_ebr", sfx), 14, 12, dir = fig_dir)
    ct <- suppressWarnings(stats::cor.test(d$ebr_used, log(d$ratio_ec_era5), method = "spearman"))
    cli::cli_inform("Part D: Spearman rho(EBR, ln ratio) = {signif(ct$estimate, 3)} (p = {signif(ct$p.value, 3)}, {nrow(d)} sites).")
  }
  cli::cli_h2("Per site")
  print_tbl(summ)
  if (nrow(lr)) {
    cli::cli_h2("Median ln(source / ERA5): pairs 1, 6, 9")
    w <- lr[lr$metric == "median_log_ratio", ]
    print_tbl(dplyr::arrange(w, site_id, stratum, pair))
  }
}

if ("D" %in% parts && !nzchar(sfx)) {
  # Supplement (G 2026-11-26): site gamma_EC / gamma_ERA5 against the measurement height, every offset found
  cli::cli_h1("Supplement: gamma_EC / gamma_ERA5 vs measurement height{tag}")
  sites_md <- read_sites(resolve_path(cfg$sites_file))
  sites_tv <- read_sites_timevarying(resolve_path(cfg$sites_timevarying_file))
  dirs <- list.files(derived_dir, pattern = "^gamma_ec(_off[0-9.]+h)?$")
  offs <- ifelse(dirs == "gamma_ec", 0, as.numeric(sub("^gamma_ec_off([0-9.]+)h$", "\\1", dirs)))
  hv <- list()
  for (k in seq_along(dirs)) for (s in ids) {
    ecd <- rd(file.path(derived_dir, dirs[k], paste0(s, "_daily.csv")))
    ecm <- rd(file.path(derived_dir, dirs[k], paste0(s, "_monthly.csv")))
    if (is.null(ecd) || is.null(ecm)) next
    m <- gamma_ec_vs_era5_table(ecm, rd(file.path(derived_dir, "gamma", paste0(s, "_monthly.csv"))), s)
    sm <- gamma_ec_site_summary(m, ecd, rd(file.path(derived_dir, "gamma", paste0(s, "_daily.csv"))),
                                rd(file.path(derived_dir, dirs[k], paste0(s, "_ebr.csv"))))
    hv[[paste(s, k)]] <- dplyr::bind_cols(
      dplyr::tibble(site_id = s, offset_h = offs[k]),
      sm[c("n_used_days", "months_with_value", "gamma_ec_site", "gamma_era5_site", "ratio_ec_era5", "ebr", "ebr_no_g")],
      gamma_ec_heights(ecd, s, sites_md, sites_tv))
  }
  hv <- dplyr::arrange(dplyr::bind_rows(hv), offset_h, site_id)
  if (nrow(hv)) {
    write_csv_header(hv, file.path(tables_dir, "17_gamma_ratio_vs_height.csv"),
                     c(paste0("Spec 17 supplement (G 2026-11-26; D1 rejected)", tag, ": one row per site x fit-start offset_h [h] (0 = standard fit; other offsets from gamma_ec_off<H>h/)."),
                       "ratio_ec_era5 = median of the monthly gamma_EC / gamma_ERA5 [-]; gamma_ec_site = median gamma_A of the used days, gamma_era5_site = median daily ERA5 gamma [K m-1]; ebr, ebr_no_g as 17_gamma_ec_sites.",
                       "ec_height_m, canopy_height_m = medians on the used days (get_site_heights(), config/sites_timevarying.csv); zm_m = EC height - 2/3 canopy height [m]; igbp = land cover (config/sites.csv)."))
    for (o in sort(unique(hv$offset_h))) {
      x <- hv[hv$offset_h == o & is.finite(hv$ratio_ec_era5), ]
      rho <- function(v) {
        ok <- is.finite(v) & v > 0
        if (sum(ok) < 4) return(NA_real_)
        unname(suppressWarnings(stats::cor(v[ok], x$ratio_ec_era5[ok], method = "spearman")))
      }
      cli::cli_inform("Offset +{o} h: {nrow(x)} sites; Spearman rho(ratio, zm) = {signif(rho(x$zm_m), 2)}, rho(ratio, EC height) = {signif(rho(x$ec_height_m), 2)}.")
      if (nrow(x)) print_tbl(dplyr::summarise(dplyr::group_by(x, igbp), n_sites = dplyr::n(),
                                              median_ratio = stats::median(ratio_ec_era5), .groups = "drop"))
    }
    if (any(is.finite(hv$ratio_ec_era5))) {
      cap <- fp_wrap(paste0("Site median of the monthly γEC / γERA5 (spec 17 supplement; decision D1 rejected). zm = EC height - 2/3 canopy height on the used days; labels: site (IGBP land cover). Filled: standard fit; open: fit start later after the first positive heat flux.", tag), 110)
      save_figure(plot_gamma_ratio_vs_height(hv, cap), "17s_gamma_ratio_vs_height", 19, 11, dir = fig_dir)
    }
  }
}

if ("E" %in% parts) {
  cli::cli_h1("Part E: bracketing of PBLHmod-gammaEC by ERA5, PBLHmod, PBLHmod-clim on the footprint rows{tag}")
  ref <- read_reference_years(resolve_path(cfg$reference_years_file))
  fp_dir <- file.path(derived_dir, "footprint")
  std_dir <- resolve_path(cfg$standard_dir, must_exist = FALSE)
  br <- list()
  for (s in ids) {
    res <- tryCatch({
      yr <- ref$year[ref$site_id == s]
      if (length(yr) != 1) cli::cli_abort("no reference year")
      p <- file.path(std_dir, s, paste0(s, "__pblhmod_gamma_ec.csv.gz"))
      if (!file.exists(p)) cli::cli_abort("standard file {.path {basename(p)}} not found (step 04 after the gamma_ec case)")
      ec <- suppressMessages(read_standard_pblh(p))
      r <- fp_site_rows(s, yr, fp_dir, sets = "etc", subsets = "valid_qc_used")
      h <- ec$pblh_m[match(as.numeric(r$timestamp_start_utc), as.numeric(ec$timestamp_start_utc))]
      gamma_ec_bracketing(r, h)
    }, error = function(e) {
      cli::cli_alert_warning("{s}: skipped ({conditionMessage(e)}).")
      NULL
    })
    br[[s]] <- res
  }
  br <- dplyr::bind_rows(br)
  if (nrow(br)) {
    write_csv_header(br, file.path(tables_dir, "17_gamma_ec_bracketing.csv"),
                     c(paste0("Spec 17 E", tag, ": rows = spec-15 common_pe of the reference year (etc, valid_qc_used: valid with an FFP result in pblhmod, pblhmod_clim, era5) with h of PBLHmod-gammaEC (standard file pblhmod_gamma_ec)."),
                       "share_outside = share of rows with h_gammaEC < min or > max of (ERA5, PBLHmod, PBLHmod-clim) h [0-1]; median_dist_ln, p90_dist_ln = ln(min / h) below, ln(h / max) above, outside rows only [-].",
                       "regime = zm/L class (unstable < -0.1 <= neutral <= 0.1 < stable), all = margin; daynight from the input table."))
    cli::cli_h2("All rows per site")
    print_tbl(br[br$regime == "all" & br$daynight == "all", ])
    cli::cli_inform("Sites with share_outside >= 5 % (all rows): {paste(br$site_id[br$regime == 'all' & br$daynight == 'all' & br$share_outside >= 0.05], collapse = ', ')}.")
  } else cli::cli_alert_warning("Part E: no site with footprint rows and the pblhmod_gamma_ec standard file.")
}
cli::cli_alert_success("Spec 17 collection done ({paste(parts, collapse = ', ')}).")
