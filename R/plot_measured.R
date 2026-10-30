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
  pal <- source_palette()
  lab <- source_labels()
  meas <- dplyr::bind_rows(
    dplyr::transmute(dm, season, hour_utc, variable = "cloudnet_censored", median = median_censored),
    dplyr::transmute(dm, season, hour_utc, variable = "cloudnet_exact", median = median_exact)
  )
  d <- dplyr::bind_rows(meas, ds[ds$variable %in% names(pal), c("season", "hour_utc", "variable", "median")])
  d$variable <- factor(d$variable, levels = intersect(c("cloudnet_censored", "cloudnet_exact", names(pal)),
                                                      unique(d$variable)))
  lt <- style_of(names(pal), "linetype")
  ggplot2::ggplot(d, ggplot2::aes(hour_utc, median, colour = variable, linetype = variable)) +
    ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, labels = lab, name = NULL) +
    ggplot2::scale_linetype_manual(values = lt, labels = lab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::labs(x = "Hour [UTC, bin centred on the hour]", y = "Median PBLH / MLH [m]") +
    theme_pblh() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), linetype = ggplot2::guide_legend(nrow = 2))
}

#' Each source vs the measured MLH (qc 0/2), by regime, log axes
#'
#' @param pt Pair tables ([pair_table()]) bound together, with a `source` column; optional logical `comparable`
#'   (spec 11 E: rows outside the comparable set drawn as faded points, counts from the comparable rows only).
#' @return ggplot object.
plot_measured_scatter <- function(pt) {
  pt <- pt[pt$mod > 0 & pt$ref > 0 & !is.na(pt$regime), ]
  if (!"comparable" %in% names(pt)) pt$comparable <- TRUE
  ctx <- geom_context_points(pt, ggplot2::aes(ref, mod))
  pt <- pt[pt$comparable %in% TRUE, ]
  ggplot2::ggplot(pt, ggplot2::aes(ref, mod)) +
    ctx +
    ggplot2::geom_bin2d(bins = 50) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = c(-1, 1) * log10(2), linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Count") +
    ggplot2::facet_grid(source ~ regime, labeller = ggplot2::labeller(regime = regime_labels())) +
    ggplot2::labs(x = "Measured MLH, qc 0/2 [m]", y = "Source [m]") +
    theme_pblh() +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1)
}

#' Share of measured half-hours by qc flag, hour x season (hour bins centred on the
#' full hour, from the half-hour midpoint)
#'
#' @param hq Standard table `cloudnet_eps_mlh` (`timestamp_start_utc`, `qc_flag`).
#' @return ggplot object.
plot_measured_qc <- function(hq) {
  qlab <- style_of(style_ids("qc"), "label", "qc")
  d <- dplyr::tibble(season = season_of(hq$timestamp_start_utc),
                     hour_utc = hour_bin_centred(hq$timestamp_start_utc + 900), # midpoint, centred bins
                     qc = factor(ifelse(is.na(hq$qc_flag), "NA", hq$qc_flag), levels = names(qlab)))
  ggplot2::ggplot(d, ggplot2::aes(hour_utc, fill = qc)) +
    ggplot2::geom_bar(position = "fill", width = 0.9) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_fill_manual(values = style_of(names(qlab), "colour", "qc"), labels = qlab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6)) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(100 * x, " %")) +
    ggplot2::labs(x = "Hour [UTC, bin centred on the hour]", y = "Share of half-hours") +
    theme_pblh() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))
}

#' log(source / measured MLH) by hour and season (qc 0/2)
#'
#' @param pt Pair tables ([pair_table()]) bound together, with a `source` column; optional logical `comparable`
#'   (spec 11 E: boxes of the rows outside the comparable set drawn faded behind).
#' @return ggplot object.
plot_measured_logratio <- function(pt) {
  pt <- pt[pt$mod > 0 & pt$ref > 0, ]
  pt$lr <- log(pt$mod / pt$ref)
  if (!"comparable" %in% names(pt)) pt$comparable <- TRUE
  ctx <- pt[!(pt$comparable %in% TRUE), ]
  ctx_layers <- if (nrow(ctx) == 0) list() else {
    ctx$.set <- comparable_set(ctx$comparable)
    list(ggplot2::geom_boxplot(data = ctx, ggplot2::aes(alpha = .data$.set), outlier.shape = NA, linewidth = 0.25,
                               position = ggplot2::position_dodge(0.8), width = 0.7),
         scale_alpha_comparable())
  }
  pt <- pt[pt$comparable %in% TRUE, ]
  ggplot2::ggplot(pt, ggplot2::aes(factor(hour_utc), lr, fill = source)) +
    ratio_lines() +
    ctx_layers +
    ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.25, position = ggplot2::position_dodge(0.8),
                          width = 0.7) +
    ggplot2::coord_cartesian(ylim = c(-3, 3)) +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_x_discrete(breaks = as.character(seq(0, 21, 3))) +
    ggplot2::labs(x = "Hour [UTC, bin centred on the hour]", y = "log(source / measured MLH) [-]", fill = NULL) +
    theme_pblh()
}

# --- Spec 10: aerosol-based references (hourly interval means) ---

# Colours and labels of the spec-10 figures: the sources plotted there (pblh_sources_style).
aerosol_palette <- function() style_of(c("sirta_mlh_chm15k", "sirta_mlh_cl31", "pblhmod_default", "era5_blh_center"))
aerosol_labels <- function() style_of(names(aerosol_palette()), "label")

#' Diurnal cycle by season with an aerosol-based reference, window marked
#'
#' Spec 11 E: the hours outside the aerosol window (the comparable set) are drawn faded.
#'
#' @param ds Long tibble: `season`, `x_utc` [h, interval midpoint or instant], `variable` (names of
#'   `aerosol_palette()`), `median` [m].
#' @param bands Tibble per season: `season`, `win_start`, `win_end` [h UTC, aerosol window].
#' @return ggplot object.
plot_aerosol_diurnal <- function(ds, bands) {
  pal <- aerosol_palette()
  ds$variable <- factor(ds$variable, levels = intersect(names(pal), unique(ds$variable)))
  b <- bands[match(as.character(ds$season), as.character(bands$season)), ]
  ds$comparable <- ds$x_utc >= b$win_start & ds$x_utc < b$win_end
  ggplot2::ggplot(ds) +
    ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = win_start, xmax = win_end, ymin = -Inf, ymax = Inf),
                       fill = style_of("aerosol_window", "colour", "element"), alpha = 0.25, inherit.aes = FALSE) +
    geom_line_comparable(ds, ggplot2::aes(x_utc, median, colour = variable), c("season", "variable"), "x_utc") +
    ggplot2::facet_wrap(~season, nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, labels = aerosol_labels(), name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24), expand = c(0, 0)) +
    ggplot2::labs(x = "Hour [UTC, hourly interval means at the midpoint; ERA5 = mean of t and t + 1 h]",
                  y = "Median PBLH / MLH [m]") +
    theme_pblh() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
}

#' Source vs aerosol-based MLH in the window, log axes, with the symmetric (RMA) line
#'
#' @param pt Pair rows with `panel` (facet label), `mod`, `ref` [m]; optional logical `comparable` (spec 11 E: the
#'   other rows are drawn as faded points; counts and lines from the comparable rows only).
#' @param cs Tibble per `panel`: `slope_sym`, `intercept_sym` [ln-ln], `slope_sym_lo`, `slope_sym_hi`, `ts_slope`,
#'   `ts_slope_reverse` [-] ([compression_slope()]).
#' @return ggplot object.
plot_compression_scatter <- function(pt, cs) {
  pt <- pt[pt$mod > 0 & pt$ref > 0, ]
  if (!"comparable" %in% names(pt)) pt$comparable <- TRUE
  rng <- range(c(pt$mod, pt$ref))
  ctx <- geom_context_points(pt, ggplot2::aes(ref, mod))
  pt <- pt[pt$comparable %in% TRUE, ]
  xs <- exp(seq(log(rng[1]), log(rng[2]), length.out = 50))
  ln <- dplyr::bind_rows(lapply(seq_len(nrow(cs)), function(i) {
    dplyr::tibble(panel = cs$panel[i], x = xs, y = exp(cs$intercept_sym[i] + cs$slope_sym[i] * log(xs)))
  }))
  cs$label <- sprintf("symmetric slope %.2f [%.2f, %.2f]\nTheil-Sen (diagnostic): %.2f, reverse %.2f",
                      cs$slope_sym, cs$slope_sym_lo, cs$slope_sym_hi, cs$ts_slope, cs$ts_slope_reverse)
  ggplot2::ggplot(pt, ggplot2::aes(ref, mod)) +
    ctx +
    ggplot2::geom_bin2d(bins = 40) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
    ggplot2::geom_abline(slope = 1, intercept = c(-1, 1) * log10(2), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(data = ln, ggplot2::aes(x, y), colour = style_of("trend", "colour", "element"), linewidth = 0.7) +
    ggplot2::geom_text(data = cs, ggplot2::aes(x = rng[1], y = rng[2], label = label), hjust = 0, vjust = 1,
                       size = 2.6) +
    ggplot2::scale_x_log10(limits = rng) +
    ggplot2::scale_y_log10(limits = rng) +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Count") +
    ggplot2::facet_wrap(~panel, nrow = 2) +
    ggplot2::labs(x = "Reference MLH [m]", y = "Source [m]") +
    theme_pblh() +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1)
}

#' Bland-Altman plot of the log ratio against the geometric mean
#'
#' `ln(source / reference)` against `exp((ln source + ln reference) / 2)` (log axis; BlandAltman1999): counts,
#' binned medians with interquartile range ([bland_altman_bins()]), and the least-squares trend with its day-block
#' bootstrap band ([compression_slope()]). A trend falling to the right is range compression; a flat trend away from
#' 0 is a constant ratio.
#'
#' @param pt Pair rows with `panel`, `mod`, `ref` [m]; optional logical `comparable` (spec 11 E: the other rows are
#'   drawn as faded points; counts from the comparable rows only).
#' @param bins Tibble with `panel` and the columns of [bland_altman_bins()].
#' @param band Tibble with `panel`, `a` [ln m], `fit`, `lo`, `hi` [-] (attribute `ba_band` of [compression_slope()]).
#' @param cs Tibble per `panel`: `ba_slope`, `ba_slope_lo`, `ba_slope_hi`, `slope_sym`, `slope_sym_lo`,
#'   `slope_sym_hi` [-].
#' @param subtitle Plot subtitle [character] (e.g. the stratum and the values used).
#' @return ggplot object.
plot_bland_altman <- function(pt, bins, band, cs, subtitle = NULL) {
  pt <- pt[pt$mod > 0 & pt$ref > 0, ]
  pt$gm <- sqrt(pt$mod * pt$ref)
  pt$d <- log(pt$mod / pt$ref)
  if (!"comparable" %in% names(pt)) pt$comparable <- TRUE
  band$gm <- exp(band$a)
  cs$label <- sprintf("trend %.3f [%.3f, %.3f] per ln m\nsymmetric slope %.2f [%.2f, %.2f]", cs$ba_slope,
                      cs$ba_slope_lo, cs$ba_slope_hi, cs$slope_sym, cs$slope_sym_lo, cs$slope_sym_hi)
  xr <- range(pt$gm)
  ctx <- geom_context_points(pt, ggplot2::aes(gm, d))
  pt <- pt[pt$comparable %in% TRUE, ]
  ggplot2::ggplot(pt, ggplot2::aes(gm, d)) +
    ctx +
    ggplot2::geom_bin2d(bins = 40) +
    ratio_lines() +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(x = gm, ymin = lo, ymax = hi), inherit.aes = FALSE,
                         fill = style_of("trend", "colour", "element"), alpha = 0.3) +
    ggplot2::geom_line(data = band, ggplot2::aes(gm, fit), colour = style_of("trend", "colour", "element"), linewidth = 0.7) +
    ggplot2::geom_errorbar(data = bins, ggplot2::aes(x = gm_m, ymin = q25_d, ymax = q75_d), inherit.aes = FALSE,
                           width = 0, linewidth = 0.4) +
    ggplot2::geom_point(data = bins, ggplot2::aes(gm_m, median_d), size = 1.6, shape = 21, fill = style_of("point_fill", "colour", "element")) +
    ggplot2::geom_text(data = cs, ggplot2::aes(x = xr[1], y = Inf, label = label), hjust = 0, vjust = 1.2,
                       size = 2.4) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "Count") +
    ggplot2::facet_wrap(~panel, nrow = 2) +
    ggplot2::labs(x = "Geometric mean of source and reference [m]", y = "ln(source / reference) [-]",
                  subtitle = subtitle) +
    theme_pblh() +
    ggplot2::theme(legend.position = "right")
}
