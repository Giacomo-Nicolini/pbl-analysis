# Figures of the network footprint summary (spec 15; paper candidates). Colours, labels and line types from
# pblh_sources_style (groups fp_pair, fp_share, footprint_case, element), theme from theme_pblh().

# Site labels: exceptions marked with " *"; network rows at the bottom.
fp_site_axis <- function(site_order, exceptions = character()) {
  lab <- ifelse(site_order %in% exceptions, paste0(site_order, " *"), site_order)
  stats::setNames(lab, site_order)
}

#' Per-site difference figure (target-area difference or ln x_80 ratio)
#'
#' p10-p90 bar and median point per site and case pair, facets regime x day/night (crossed strata), the network
#' (pooled rows) as the bottom row, exceptions marked with " *".
#'
#' @param x Output of [fp_target_area_table()] or [fp_fetch_table()] (one set and subset).
#' @param lo,mid,hi Column names of p10, median, p90 [character].
#' @param site_order Sites from top to bottom [character].
#' @param exceptions Sites to mark [character].
#' @param xlab Axis label [character].
#' @param subtitle Subtitle [character].
#' @return A ggplot.
plot_fp_site_diff <- function(x, lo, mid, hi, site_order, exceptions = character(), xlab, subtitle = NULL) {
  d <- x[x$scope %in% c("site", "network_pooled") & x$regime != "all" & x$daynight != "all", ]
  lev <- c("network", rev(site_order)) # first level at the bottom: network last, first site on top
  labs <- c(rev(fp_site_axis(site_order, exceptions)), network = "Network (pooled)")
  d$site <- factor(d$site_id, levels = lev)
  d <- d[!is.na(d$site), ]
  d$regime <- factor(d$regime, levels = c("unstable", "neutral", "stable"))
  d$daynight <- factor(d$daynight, levels = c("day", "night"))
  pr <- style_ids("fp_pair", unique(d$pair))
  d$pair <- factor(d$pair, levels = pr)
  pos <- ggplot2::position_dodge(width = 0.7)
  ggplot2::ggplot(d, ggplot2::aes(.data[[mid]], site, colour = pair)) +
    ggplot2::geom_vline(xintercept = 0, colour = style_of("reference_line", "colour", "element"), linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data[[lo]], xmax = .data[[hi]]), width = 0, linewidth = 0.35, orientation = "y",
                            position = pos) +
    ggplot2::geom_point(size = 1.1, position = pos) +
    ggplot2::facet_grid(daynight ~ regime) +
    ggplot2::scale_y_discrete(labels = labs, drop = TRUE) +
    ggplot2::scale_colour_manual(values = style_of(pr, "colour", "fp_pair"), labels = style_of(pr, "label", "fp_pair"),
                                 name = NULL) +
    ggplot2::labs(x = xlab, y = NULL, subtitle = subtitle) +
    theme_pblh()
}

#' Top-affected rows: stacked shares per site
#'
#' @param top Output of [fp_topaffected_table()] (one set and pair).
#' @param site_order Sites from top to bottom [character].
#' @param subtitle Subtitle [character].
#' @return A ggplot: three panels (stability class, day / night, removal by QC / used), one stacked bar per site and
#'   one for the network with the pooled shares, so that each bar sums to 100 % (fix-14c, G 2026-11-12; the median of
#'   the site shares stays in the table).
plot_fp_topaffected <- function(top, site_order, subtitle = NULL) {
  d <- top[top$scope %in% c("site", "network_pooled"), ]
  d$row <- ifelse(d$scope == "site", d$site_id, "net_pooled")
  lev <- c("net_pooled", rev(site_order)) # network row at the bottom
  labs <- c(stats::setNames(rev(site_order), rev(site_order)), net_pooled = "Network (pooled)")
  long <- dplyr::bind_rows(
    dplyr::tibble(row = d$row, panel = "stability (zm/L)", cat = "zl_unstable", share = d$share_unstable),
    dplyr::tibble(row = d$row, panel = "stability (zm/L)", cat = "zl_neutral", share = d$share_neutral),
    dplyr::tibble(row = d$row, panel = "stability (zm/L)", cat = "zl_stable", share = d$share_stable),
    dplyr::tibble(row = d$row, panel = "day / night", cat = "day", share = 1 - d$share_night),
    dplyr::tibble(row = d$row, panel = "day / night", cat = "night", share = d$share_night),
    dplyr::tibble(row = d$row, panel = "removal", cat = "kept", share = 1 - d$share_removed_by_qc_or_used),
    dplyr::tibble(row = d$row, panel = "removal", cat = "removed_used",
                  share = d$share_removed_by_qc_or_used - d$share_removed_by_qc),
    dplyr::tibble(row = d$row, panel = "removal", cat = "removed_qc", share = d$share_removed_by_qc))
  long$row <- factor(long$row, levels = lev)
  long <- long[!is.na(long$row), ]
  long$panel <- factor(long$panel, levels = c("stability (zm/L)", "day / night", "removal"))
  cats <- style_ids("fp_share")
  long$cat <- factor(long$cat, levels = rev(cats))
  ggplot2::ggplot(long, ggplot2::aes(100 * share, row, fill = cat)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_y_discrete(labels = labs) +
    ggplot2::scale_fill_manual(values = style_of(cats, "colour", "fp_share"), labels = style_of(cats, "label", "fp_share"),
                               breaks = cats, name = NULL) +
    ggplot2::labs(x = "Share of the top 10 % half-hours [%]", y = NULL, subtitle = subtitle) +
    theme_pblh() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3))
}

#' ICOS criteria per site and case
#'
#' @param ic Output of [fp_icos_criteria_table()].
#' @param site_order Sites from top to bottom [character].
#' @param subtitle Subtitle [character].
#' @return A ggplot: share of the half-hours with >= 90 % from the target area (all rows; stable rows; 50 % line,
#'   Rebmann et al. 2018) and with >= 70 % among the good data (ICOS labelling rule, 70 % line), points per case.
plot_fp_icos_criteria <- function(ic, site_order, subtitle = NULL) {
  pn <- c("TA >= 90 %, all rows (Rebmann 2018)", "TA >= 90 %, stable rows", "TA >= 70 %, good data (ICOS labelling)")
  long <- dplyr::bind_rows(
    dplyr::tibble(site_id = ic$site_id, case = ic$case, panel = pn[1], share = ic$share_ta90),
    dplyr::tibble(site_id = ic$site_id, case = ic$case, panel = pn[2], share = ic$share_ta90_stable),
    dplyr::tibble(site_id = ic$site_id, case = ic$case, panel = pn[3], share = ic$icos_lab_share))
  long$site <- factor(long$site_id, levels = rev(site_order))
  long <- long[!is.na(long$site), ]
  long$panel <- factor(long$panel, levels = pn)
  cs <- style_ids("footprint_case", unique(long$case))
  ln <- dplyr::tibble(panel = factor(pn, levels = pn), x = c(50, 50, 70), id = c("ta90_line", "ta90_line", "lab70_line"))
  ggplot2::ggplot(long, ggplot2::aes(100 * share, site, colour = case)) +
    ggplot2::geom_vline(data = ln, ggplot2::aes(xintercept = x), colour = style_of(ln$id, "colour", "element"),
                        linetype = style_of(ln$id, "linetype", "element"), linewidth = 0.4) +
    ggplot2::geom_point(size = 1.3, alpha = 0.85, na.rm = TRUE) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_colour_manual(values = style_of(cs, "colour", "footprint_case"),
                                 labels = style_of(cs, "label", "footprint_case"), breaks = cs, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 100)) +
    ggplot2::labs(x = "Share of the half-hours [%]", y = NULL, subtitle = subtitle) +
    theme_pblh()
}

#' Site gamma with IQR, the network gamma and the literature constant
#'
#' @param st Output of [gamma_site_table()].
#' @param net Output of [gamma_network_table()].
#' @return A ggplot.
plot_gamma_site <- function(st, net) {
  st <- st[!is.na(st$gamma_site), ]
  st$site <- factor(st$site_id, levels = st$site_id[order(st$gamma_site)])
  ln <- dplyr::tibble(x = c(net$gamma_med, 0.01), id = c("pblhmod_gamma_med", "pblhmod_gamma001"))
  ggplot2::ggplot(st, ggplot2::aes(gamma_site, site)) +
    ggplot2::annotate("rect", xmin = net$gamma_med_lo95, xmax = net$gamma_med_hi95, ymin = -Inf, ymax = Inf,
                      fill = style_of("pblhmod_gamma_med"), alpha = 0.15) +
    ggplot2::geom_vline(data = ln, ggplot2::aes(xintercept = x, colour = id), linewidth = 0.5) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = p25, xmax = p75), width = 0, linewidth = 0.4, orientation = "y",
                            colour = style_of("context", "colour", "element")) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = style_of(ln$id), labels = c(pblhmod_gamma_med = "γ network median (95 % site bootstrap)",
                                                                       pblhmod_gamma001 = "γ = 0.01 K m-1 (literature)"), name = NULL) +
    ggplot2::labs(x = "Daily morning γ [K m-1] (ERA5, 500-3500 m above the site): site median and IQR", y = NULL) +
    theme_pblh()
}
