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
    ggplot2::labs(x = "Hour [UTC, bin centred on the hour]", y = "Median PBLH / MLH [m]") +
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

#' Share of measured half-hours by qc flag, hour x season (hour bins centred on the
#' full hour, from the half-hour midpoint)
#'
#' @param hq Standard table `cloudnet_eps_mlh` (`timestamp_start_utc`, `qc_flag`).
#' @return ggplot object.
plot_measured_qc <- function(hq) {
  qlab <- c(`0` = "0 measured", `2` = "2 partly censored", `3` = "3 signal-limited (lower bound)",
            `1` = "1 censored below z_min", `NA` = "NA coverage < min")
  d <- dplyr::tibble(season = season_of(hq$timestamp_start_utc),
                     hour_utc = hour_bin_centred(hq$timestamp_start_utc + 900), # midpoint, centred bins
                     qc = factor(ifelse(is.na(hq$qc_flag), "NA", hq$qc_flag), levels = names(qlab)))
  ggplot2::ggplot(d, ggplot2::aes(hour_utc, fill = qc)) +
    ggplot2::geom_bar(position = "fill", width = 0.9) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_fill_manual(values = c(`0` = "#0072B2", `2` = "#56B4E9", `3` = "#E69F00", `1` = "#999999",
                                          `NA` = "#DDDDDD"), labels = qlab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(100 * x, " %")) +
    ggplot2::labs(x = "Hour [UTC, bin centred on the hour]", y = "Share of half-hours") +
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
    ggplot2::labs(x = "Hour [UTC, bin centred on the hour]", y = "log(source / measured MLH) [-]", fill = NULL) +
    theme_pub()
}

# --- Spec 10: aerosol-based references (hourly interval means) ---

# Colours and labels of the spec-10 figures (Okabe & Ito 2008, as source_palette(); kept separate so that the
# spec-07 figures are unchanged).
aerosol_palette <- function() {
  c(sirta_mlh_chm15k = "black", sirta_mlh_cl31 = "#CC79A7", pblhmod_default = "#0072B2",
    era5_blh_center = "#D55E00")
}
aerosol_labels <- function() {
  c(sirta_mlh_chm15k = "SIRTA CHM15k (STRATfinder)", sirta_mlh_cl31 = "SIRTA CL31 (CABAM)",
    pblhmod_default = "PBLHmod (default)", era5_blh_center = "ERA5 (tower cell)")
}

#' Diurnal cycle by season with an aerosol-based reference, window and night marked
#'
#' @param ds Long tibble: `season`, `x_utc` [h, interval midpoint or instant], `variable` (names of
#'   `aerosol_palette()`), `median` [m].
#' @param bands Tibble per season: `season`, `win_start`, `win_end` [h UTC, aerosol window], `sunrise`, `sunset`
#'   [h UTC, median over the season].
#' @return ggplot object.
plot_aerosol_diurnal <- function(ds, bands) {
  pal <- aerosol_palette()
  ds$variable <- factor(ds$variable, levels = intersect(names(pal), unique(ds$variable)))
  ymax <- max(ds$median, na.rm = TRUE)
  bands$night_label <- "night:\nnot comparable"
  ggplot2::ggplot(ds) +
    ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = win_start, xmax = win_end, ymin = -Inf, ymax = Inf),
                       fill = "#F0E442", alpha = 0.25, inherit.aes = FALSE) +
    ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = 0, xmax = sunrise, ymin = -Inf, ymax = Inf),
                       fill = "grey60", alpha = 0.2, inherit.aes = FALSE) +
    ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = sunset, xmax = 24, ymin = -Inf, ymax = Inf),
                       fill = "grey60", alpha = 0.2, inherit.aes = FALSE) +
    ggplot2::geom_text(data = bands, ggplot2::aes(x = 0.3, y = ymax * 0.97, label = night_label),
                       size = 2.3, colour = "grey30", hjust = 0, vjust = 1, inherit.aes = FALSE) +
    ggplot2::geom_line(ggplot2::aes(x_utc, median, colour = variable), linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, labels = aerosol_labels(), name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24), expand = c(0, 0)) +
    ggplot2::labs(x = "Hour [UTC, hourly interval means at the midpoint; ERA5 = mean of t and t + 1 h]",
                  y = "Median PBLH / MLH [m]") +
    theme_pub() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
}

#' Source vs aerosol-based MLH in the window, log axes, with the Theil-Sen line
#'
#' @param pt Pair rows with `panel` (facet label), `mod`, `ref` [m].
#' @param cs Tibble per `panel`: `slope`, `intercept` [ln-ln], `slope_lo`, `slope_hi`.
#' @return ggplot object.
plot_compression_scatter <- function(pt, cs) {
  pt <- pt[pt$mod > 0 & pt$ref > 0, ]
  rng <- range(c(pt$mod, pt$ref))
  xs <- exp(seq(log(rng[1]), log(rng[2]), length.out = 50))
  ln <- dplyr::bind_rows(lapply(seq_len(nrow(cs)), function(i) {
    dplyr::tibble(panel = cs$panel[i], x = xs, y = exp(cs$intercept[i] + cs$slope[i] * log(xs)))
  }))
  cs$label <- sprintf("Theil-Sen slope %.2f [%.2f, %.2f]", cs$slope, cs$slope_lo, cs$slope_hi)
  ggplot2::ggplot(pt, ggplot2::aes(ref, mod)) +
    ggplot2::geom_bin2d(bins = 40) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = c(-1, 1) * log10(2), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(data = ln, ggplot2::aes(x, y), colour = "#D55E00", linewidth = 0.7) +
    ggplot2::geom_text(data = cs, ggplot2::aes(x = rng[1], y = rng[2], label = label), hjust = 0, vjust = 1,
                       size = 2.6) +
    ggplot2::scale_x_log10(limits = rng) +
    ggplot2::scale_y_log10(limits = rng) +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Count") +
    ggplot2::facet_wrap(~panel, nrow = 2) +
    ggplot2::labs(x = "Reference MLH [m]", y = "Source [m]") +
    theme_pub() +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1)
}

#' Median log ratio by tercile of the measured MLH
#'
#' @param lt Tibble: `pair_label`, `stratum_set`, `tercile`, `ref_lo`, `ref_hi` [m], `median_log_ratio` [-], `n`.
#' @return ggplot object.
plot_logratio_tercile <- function(lt) {
  ggplot2::ggplot(lt, ggplot2::aes(tercile, median_log_ratio, colour = stratum_set, group = stratum_set)) +
    ratio_lines() +
    ggplot2::geom_point(size = 1.8, position = ggplot2::position_dodge(0.4)) +
    ggplot2::geom_line(linewidth = 0.4, position = ggplot2::position_dodge(0.4)) +
    ggplot2::facet_wrap(~pair_label, nrow = 2) +
    ggplot2::scale_colour_manual(values = c(window_convective = "black", convective = "grey55"),
                                 labels = c(window_convective = "Window, convective (primary)",
                                            convective = "All convective hours"), name = NULL) +
    ggplot2::labs(x = "Tercile of the measured MLH (T1 shallowest)", y = "Median ln(source / reference) [-]") +
    theme_pub()
}
