# Figures of the first PBLHmod vs ERA5 comparison (spec 04 section 8).
# Functions return ggplot objects (or lists of them); saving: save_figure().

#' Diurnal cycle by season (figure 1)
#'
#' Median and IQR ribbon per hour UTC, one panel per season.
#'
#' @param ds Output of [diurnal_summary()] with `variable` in the names of
#'   [source_palette()] (`icos_oper_all` drawn dashed, without ribbon).
#' @return ggplot object.
plot_diurnal_by_season <- function(ds) {
  pal <- source_palette()
  lab <- source_labels()
  ds <- ds[ds$variable %in% names(pal), ]
  ds$variable <- factor(ds$variable, levels = intersect(names(pal), unique(ds$variable)))
  rib <- ds[ds$variable != "icos_oper_all", ]
  ggplot2::ggplot(ds, ggplot2::aes(hour_utc, median, colour = variable, linetype = variable)) +
    ggplot2::geom_ribbon(data = rib, ggplot2::aes(ymin = p25, ymax = p75, fill = variable),
                         alpha = 0.15, colour = NA, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, labels = lab, name = NULL) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_linetype_manual(values = c(pblhmod_default = "solid", pblhmod_clim = "solid",
                                              pblhmod_gamma001 = "solid", pblhmod_gamma0025 = "solid",
                                              era5_blh_center = "solid", era5_blh_bilinear = "solid",
                                              icos_oper_qc0 = "solid", icos_oper_all = "dashed"),
                                   labels = lab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::labs(x = "Hour [UTC]", y = "PBLH [m]") +
    theme_pub() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), linetype = ggplot2::guide_legend(nrow = 2))
}

#' 2-D density scatter of a pair on log axes (figure 2)
#'
#' 1:1 line (solid) and factor-2 lines (dashed); one panel per regime.
#'
#' @param pt Output of [pair_table()] (one pair).
#' @param xlab,ylab Axis labels.
#' @return ggplot object.
plot_scatter_pair <- function(pt, xlab = "ERA5 BLH [m]", ylab = "PBLHmod [m]") {
  pt <- pt[pt$mod > 0 & pt$ref > 0 & !is.na(pt$regime), ]
  ggplot2::ggplot(pt, ggplot2::aes(ref, mod)) +
    ggplot2::geom_bin2d(bins = 60) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = c(-1, 1) * log10(2), linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Hours") +
    ggplot2::facet_wrap(~regime, nrow = 1, labeller = ggplot2::as_labeller(regime_labels())) +
    ggplot2::labs(x = xlab, y = ylab) +
    theme_pub() +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1)
}

#' Log ratio by hour UTC and season (figure 3)
#'
#' @param pt Output of [pair_table()] (one pair).
#' @param ylab Axis label.
#' @return ggplot object.
plot_logratio_by_hour <- function(pt, ylab = "log(PBLHmod / ERA5) [-]") {
  pt <- pt[pt$mod > 0 & pt$ref > 0, ]
  pt$lr <- log(pt$mod / pt$ref)
  ggplot2::ggplot(pt, ggplot2::aes(factor(hour_utc), lr)) +
    ratio_lines() +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey85", linewidth = 0.3) +
    ggplot2::coord_cartesian(ylim = c(-3, 3)) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_x_discrete(breaks = as.character(seq(0, 21, 3))) +
    ggplot2::labs(x = "Hour [UTC]", y = ylab) +
    theme_pub()
}

#' ERA5 3x3 spread and cell sensitivity (figure 4)
#'
#' @param era_c Standard table `era5_blh_center` (`timestamp_start_utc`,
#'   `pblh_m`, `x_blh_cv9`).
#' @param era_b Standard table `era5_blh_bilinear`.
#' @return List of two ggplot objects: cv9 per hour, centre vs bilinear.
plot_era5_spread <- function(era_c, era_b) {
  a <- dplyr::tibble(hour_utc = as.integer(format(era_c$timestamp_start_utc, "%H", tz = "UTC")),
                     cv9 = era_c$x_blh_cv9)
  p1 <- ggplot2::ggplot(a[!is.na(a$cv9), ], ggplot2::aes(factor(hour_utc), cv9)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey85", linewidth = 0.3) +
    ggplot2::coord_cartesian(ylim = c(0, 0.6)) +
    ggplot2::scale_x_discrete(breaks = as.character(seq(0, 21, 3))) +
    ggplot2::labs(x = "Hour [UTC]", y = "CV of ERA5 BLH over 3x3 cells [-]", title = "a") +
    theme_pub()
  b <- dplyr::inner_join(
    dplyr::tibble(t = era_c$timestamp_start_utc, center = era_c$pblh_m),
    dplyr::tibble(t = era_b$timestamp_start_utc, bilin = era_b$pblh_m),
    by = "t"
  )
  b <- b[b$center > 0 & b$bilin > 0 & !is.na(b$center) & !is.na(b$bilin), ]
  p2 <- ggplot2::ggplot(b, ggplot2::aes(center, bilin)) +
    ggplot2::geom_bin2d(bins = 60) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Hours") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "ERA5 BLH, tower cell [m]", y = "ERA5 BLH, bilinear [m]", title = "b") +
    theme_pub() +
    ggplot2::theme(legend.position = "right")
  list(p1, p2)
}

#' EC H vs ERA5 H (figure 5)
#'
#' @param hd Tibble with `timestamp_utc` (end of the hour), `H_ec_Wm2`,
#'   `H_era5_Wm2` [W m-2], `season`.
#' @return List of two ggplot objects: scatter, mean diurnal cycles by season.
plot_h_ec_vs_era5 <- function(hd) {
  hd <- hd[!is.na(hd$H_ec_Wm2) & !is.na(hd$H_era5_Wm2), ]
  lim <- range(c(hd$H_ec_Wm2, hd$H_era5_Wm2))
  p1 <- ggplot2::ggplot(hd, ggplot2::aes(H_era5_Wm2, H_ec_Wm2)) +
    ggplot2::geom_bin2d(bins = 70) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Hours") +
    ggplot2::coord_equal(xlim = lim, ylim = lim) +
    ggplot2::labs(x = "ERA5 H, tower cell [W m-2]", y = "EC H (PBLHmod forcing) [W m-2]", title = "a") +
    theme_pub() +
    ggplot2::theme(legend.position = "right")
  hd$hour_mid <- (as.integer(format(hd$timestamp_utc, "%H", tz = "UTC")) - 0.5) %% 24
  cyc <- dplyr::summarise(hd, .by = c(season, hour_mid), EC = mean(H_ec_Wm2), ERA5 = mean(H_era5_Wm2))
  cyc <- dplyr::bind_rows(
    dplyr::tibble(season = cyc$season, hour_mid = cyc$hour_mid, source = "EC", H = cyc$EC),
    dplyr::tibble(season = cyc$season, hour_mid = cyc$hour_mid, source = "ERA5", H = cyc$ERA5)
  )
  p2 <- ggplot2::ggplot(cyc[order(cyc$hour_mid), ], ggplot2::aes(hour_mid, H, colour = source)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = c(EC = unname(source_palette()["pblhmod_default"]),
                                            ERA5 = unname(source_palette()["era5_blh_center"])), name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::labs(x = "Hour [UTC, interval midpoint]", y = "Mean H [W m-2]", title = "b") +
    theme_pub()
  list(p1, p2)
}

#' Log ratio by cloud-base class and regime (figure 6)
#'
#' @param pc Pair table with `mod`, `ref`, `regime` and `cloud` (output of
#'   [cloud_class()]).
#' @return ggplot object.
plot_logratio_by_cloud <- function(pc) {
  pc <- pc[pc$mod > 0 & pc$ref > 0 & !is.na(pc$regime), ]
  pc$lr <- log(pc$mod / pc$ref)
  ggplot2::ggplot(pc, ggplot2::aes(cloud, lr)) +
    ratio_lines() +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey85", linewidth = 0.3, varwidth = TRUE) +
    ggplot2::coord_cartesian(ylim = c(-3, 3)) +
    ggplot2::facet_wrap(~regime, labeller = ggplot2::as_labeller(regime_labels())) +
    ggplot2::scale_x_discrete(labels = c(no_cloud = "no cloud", cbh_le_2000 = "base <= 2 km",
                                         cbh_gt_2000 = "base > 2 km")) +
    ggplot2::labs(x = "ERA5 cloud-base class", y = "log(PBLHmod / ERA5) [-]") +
    theme_pub()
}

#' Monthly climatology of gamma_850_700 (figure 7)
#'
#' Median and IQR per calendar month at 00 and 12 UTC; reference lines at the
#' constant gamma values used in PBLHmod (0.01 and 0.0025 K m-1).
#'
#' @param gm Output of [gamma_monthly()].
#' @return ggplot object.
plot_gamma_monthly <- function(gm) {
  gm$hour <- factor(sprintf("%02d UTC", gm$hour_utc))
  ggplot2::ggplot(gm, ggplot2::aes(month, median, colour = hour, fill = hour)) +
    ggplot2::geom_hline(yintercept = c(0.01, 0.0025), linetype = "dashed", linewidth = 0.3) +
    ggplot2::annotate("text", x = 12.4, y = c(0.01, 0.0025), label = c("0.01", "0.0025"),
                      hjust = 0, vjust = -0.3, size = 3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = p25, ymax = p75), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_x_continuous(breaks = 1:12, labels = month.abb, limits = c(1, 12.9)) +
    ggplot2::scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
    ggplot2::scale_fill_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
    ggplot2::labs(x = NULL, y = expression(gamma[850-700] ~ "[K m"^-1 * "]")) +
    theme_pub()
}

# Reference lines for log ratios: 0 (solid) and +-log(2) (factor 2, dashed).
ratio_lines <- function() {
  list(
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3),
    ggplot2::geom_hline(yintercept = c(-log(2), log(2)), linetype = "dashed", linewidth = 0.3)
  )
}

# Readable PBLHmod regime names for facet strips.
regime_labels <- function() {
  c(stable = "stable", unstable_diag = "unstable (diagnostic)", convective = "convective (B&G)",
    mixed = "mixed (transition hour)")
}
