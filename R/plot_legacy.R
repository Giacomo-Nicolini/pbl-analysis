# Figures of the legacy regression (spec 03 section 6). Functions return ggplot objects
# (or lists of them); saving: save_figure() (R/plot_utils.R).

#' Reproduction of the operational PBLH (figure 1)
#'
#' @param t Times [POSIXct UTC].
#' @param oper Operational PBLH [m].
#' @param leg Legacy reproduction [m].
#' @return List of two ggplot objects: 2-D density (log axes, 1:1 line) and the
#'   difference time series.
plot_legacy_reproduction <- function(t, oper, leg) {
  d <- dplyr::tibble(t = t, oper = oper, leg = leg, diff = leg - oper)
  ok <- d[!is.na(d$oper) & !is.na(d$leg) & d$oper > 0 & d$leg > 0, ]
  p1 <- ggplot2::ggplot(ok, ggplot2::aes(oper, leg)) +
    ggplot2::geom_bin2d(bins = 80) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.3) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Half-hours") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Operational L2 PBLH [m]", y = "Legacy code on L2 inputs [m]", title = "a") +
    theme_pub() + ggplot2::theme(legend.position = "right")
  p2 <- ggplot2::ggplot(d[!is.na(d$diff), ], ggplot2::aes(t, diff)) +
    ggplot2::geom_point(size = 0.2, alpha = 0.4) +
    ggplot2::labs(x = NULL, y = "Legacy - operational [m]", title = "b") +
    theme_pub()
  list(p1, p2)
}

#' Median PBLH by regime along the cumulative path (figure 2)
#'
#' @param med Tibble `step` (ordered factor), `regime`, `median` [m].
#' @param ref Tibble `regime`, `median` [m] of PBLHmod v2 (reference lines).
#' @return ggplot object.
plot_legacy_waterfall <- function(med, ref) {
  ggplot2::ggplot(med, ggplot2::aes(step, median, colour = regime, group = regime)) +
    ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = median, colour = regime), linetype = "dashed",
                        linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = c(all = "black", stable = "#0072B2", unstable_diag = "#E69F00",
                                            convective = "#D55E00"), name = NULL) +
    ggplot2::labs(x = "Cumulative change from the operational legacy code (dashed: PBLHmod v2)", y = "Median PBLH [m]") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
}

#' Median log ratio to PBLHmod v2 per switch, one at a time (figure 3)
#'
#' @param lr Tibble `variant`, `regime`, `median_log_ratio` [-].
#' @return ggplot object.
plot_legacy_one_at_a_time <- function(lr) {
  ggplot2::ggplot(lr, ggplot2::aes(median_log_ratio, variant, colour = regime)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_colour_manual(values = c(all = "black", stable = "#0072B2", unstable_diag = "#E69F00",
                                            convective = "#D55E00"), name = NULL) +
    ggplot2::labs(x = "Median log(variant / PBLHmod v2) [-]", y = NULL) +
    theme_pub()
}

#' Median diurnal cycle by season (figure 4)
#'
#' @param ds Output of [diurnal_summary()] with `variable` in `operational`,
#'   `legacy_baseline`, `legacy_no_postproc`, `pblhmod_v2`.
#' @return ggplot object.
plot_legacy_diurnal <- function(ds) {
  pal <- c(operational = "#009E73", legacy_baseline = "#CC79A7", legacy_no_postproc = "#E69F00", pblhmod_v2 = "#0072B2")
  lab <- c(operational = "Operational L2", legacy_baseline = "Legacy code (L2 inputs)",
           legacy_no_postproc = "Legacy without fill / floor / smoothing", pblhmod_v2 = "PBLHmod v2")
  ds$variable <- factor(ds$variable, levels = names(pal))
  ggplot2::ggplot(ds, ggplot2::aes(hour_utc, median, colour = variable, linetype = variable)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, labels = lab, name = NULL) +
    ggplot2::scale_linetype_manual(values = c(operational = "solid", legacy_baseline = "22", legacy_no_postproc = "solid",
                                              pblhmod_v2 = "solid"), labels = lab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::labs(x = "Hour [UTC]", y = "Median PBLH [m]") +
    theme_pub() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), linetype = ggplot2::guide_legend(nrow = 2))
}
