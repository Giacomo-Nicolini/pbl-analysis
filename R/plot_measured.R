# Figures of the comparison with the measured MLH (spec 07 section 5).
# Functions return ggplot objects; saving: save_figure().

#' Diurnal cycle by season: measured MLH (censored and plain medians) and other sources
#'
#' @param dm Output of [measured_diurnal()] (`season`, `hour_utc`,
#'   `median_censored`, `median_exact`).
#' @param ds Output of [diurnal_summary()] for the other sources (`variable` in
#'   the names of [source_palette()]).
#' @return ggplot object.
plot_measured_diurnal <- function(dm, ds) {
  pal <- c(source_palette(), cloudnet_censored = "black", cloudnet_exact = "grey50")
  lab <- c(source_labels(), cloudnet_censored = "Measured MLH (censored median)",
           cloudnet_exact = "Measured MLH (qc 0/2 median)")
  meas <- dplyr::bind_rows(
    dplyr::transmute(dm, season, hour_utc, variable = "cloudnet_censored", median = median_censored),
    dplyr::transmute(dm, season, hour_utc, variable = "cloudnet_exact", median = median_exact)
  )
  d <- dplyr::bind_rows(meas, ds[ds$variable %in% names(pal), c("season", "hour_utc", "variable", "median")])
  d$variable <- factor(d$variable, levels = intersect(c("cloudnet_censored", "cloudnet_exact", names(pal)),
                                                      unique(d$variable)))
  lt <- stats::setNames(rep("solid", length(pal)), names(pal))
  lt["cloudnet_exact"] <- "dashed"
  ggplot2::ggplot(d, ggplot2::aes(hour_utc, median, colour = variable, linetype = variable)) +
    ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, labels = lab, name = NULL) +
    ggplot2::scale_linetype_manual(values = lt, labels = lab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::labs(x = "Hour [UTC]", y = "Median PBLH / MLH [m]") +
    theme_pub() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), linetype = ggplot2::guide_legend(nrow = 2))
}

#' Each source vs the measured MLH (qc 0/2), by regime, log axes
#'
#' @param pt Pair tables ([pair_table()]) bound together, with a `source` column.
#' @return ggplot object.
plot_measured_scatter <- function(pt) {
  pt <- pt[pt$mod > 0 & pt$ref > 0 & !is.na(pt$regime), ]
  ggplot2::ggplot(pt, ggplot2::aes(ref, mod)) +
    ggplot2::geom_bin2d(bins = 50) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = c(-1, 1) * log10(2), linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Count") +
    ggplot2::facet_grid(source ~ regime, labeller = ggplot2::labeller(regime = regime_labels())) +
    ggplot2::labs(x = "Measured MLH, qc 0/2 [m]", y = "Source [m]") +
    theme_pub() +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1)
}

#' Share of measured half-hours by qc flag, hour x season
#'
#' @param hq Standard table `cloudnet_eps_mlh` (`timestamp_start_utc`, `qc_flag`).
#' @return ggplot object.
plot_measured_qc <- function(hq) {
  qlab <- c(`0` = "0 measured", `2` = "2 partly censored", `3` = "3 signal-limited (lower bound)",
            `1` = "1 censored below z_min", `NA` = "NA coverage < min")
  d <- dplyr::tibble(season = season_of(hq$timestamp_start_utc),
                     hour_utc = as.integer(format(hq$timestamp_start_utc, "%H", tz = "UTC")),
                     qc = factor(ifelse(is.na(hq$qc_flag), "NA", hq$qc_flag), levels = names(qlab)))
  ggplot2::ggplot(d, ggplot2::aes(hour_utc, fill = qc)) +
    ggplot2::geom_bar(position = "fill", width = 0.9) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_fill_manual(values = c(`0` = "#0072B2", `2` = "#56B4E9", `3` = "#E69F00", `1` = "#999999",
                                          `NA` = "#DDDDDD"), labels = qlab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(100 * x, " %")) +
    ggplot2::labs(x = "Hour [UTC]", y = "Share of half-hours") +
    theme_pub() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))
}

#' log(source / measured MLH) by hour and season (qc 0/2)
#'
#' @param pt Pair tables ([pair_table()]) bound together, with a `source` column.
#' @return ggplot object.
plot_measured_logratio <- function(pt) {
  pt <- pt[pt$mod > 0 & pt$ref > 0, ]
  pt$lr <- log(pt$mod / pt$ref)
  ggplot2::ggplot(pt, ggplot2::aes(factor(hour_utc), lr, fill = source)) +
    ratio_lines() +
    ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.25, position = ggplot2::position_dodge(0.8),
                          width = 0.7) +
    ggplot2::coord_cartesian(ylim = c(-3, 3)) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_x_discrete(breaks = as.character(seq(0, 21, 3))) +
    ggplot2::labs(x = "Hour [UTC]", y = "log(source / measured MLH) [-]", fill = NULL) +
    theme_pub()
}
