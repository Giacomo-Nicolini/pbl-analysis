# Figures of the footprint uncertainty (spec 16). Colours, labels, order and line types from pblh_sources_style
# (groups mc_input, element); theme from theme_pblh().

# Facet factors of the strata (regime x day/night); `all` dropped.
mc_strata_facets <- function(d) {
  d <- d[d$stratum %in% mc_strata_levels(), ]
  d$regime <- factor(sub("_.*$", "", d$stratum), levels = c("unstable", "neutral", "stable"))
  d$daynight <- factor(sub("^.*_", "", d$stratum), levels = c("day", "night"))
  d
}

#' Input ranking: 1-sigma effect of each input per stratum, with the PBLH-source effect as reference
#'
#' @param s Bound site and network summaries ([mc_site_summary()], [mc_network_summary()]).
#' @param quantity `"ta_pp"` (TA fraction [pp]) or `"lnx80"` (ln x80 [-]).
#' @param subtitle,caption Texts (wrapped).
#' @param rsl_off Inputs whose effect is taken from the runs without the e_rsl test (items `<input>|rsl_off`, spec 16
#'   amendment 1), labelled "rsl test off". Their share of runs failing e_rsl under the etc set is shown per site in
#'   [plot_mc_validity_flips()] (close-out 2026-11-19: no network median here).
#' @return A ggplot: inputs on y (ordered by the network median over all strata, largest on top), the 1-sigma effect
#'   on x; sites as small points, network as the large point with the IQR over sites; vertical reference lines at the
#'   network median |ERA5 - PBLHmod| (and |PBLHmod-clim - PBLHmod|) of the stratum. One legend below the panels:
#'   inputs by colour (points) and the PBLH-source references (lines).
plot_mc_input_ranking <- function(s, quantity = "ta_pp", subtitle = NULL, caption = NULL, rsl_off = c("hc", "z0f")) {
  d <- s[s$kind == "oat_1sigma" & s$quantity == quantity, ]
  ro <- paste0(rsl_off, "|rsl_off")
  have <- rsl_off[ro %in% d$item]
  d <- d[!d$item %in% have, ] # the etc rows of these inputs replaced by their rsl_off rows
  d$item[d$item %in% paste0(have, "|rsl_off")] <- sub("\\|rsl_off$", "", d$item[d$item %in% paste0(have, "|rsl_off")])
  net_all <- d[d$scope == "network" & d$stratum == "all", ]
  d <- mc_strata_facets(d)
  src_all <- c("mc_source_era5", "mc_source_clim")
  ref <- mc_strata_facets(s[s$kind == "source_effect" & s$quantity == quantity & s$scope == "network", ])
  ref$id <- ifelse(ref$item == "era5_vs_pblhmod", "mc_source_era5", "mc_source_clim")
  ref <- ref[!is.na(ref$median), ]
  src <- src_all[src_all %in% ref$id]
  ids <- style_ids("mc_input", unique(d$item))
  ord <- ids[order(net_all$median[match(ids, net_all$item)], na.last = FALSE)]
  labs <- stats::setNames(style_of(ord, "label", "mc_input"), ord)
  labs[ord %in% have] <- paste0(labs[ord %in% have], ", ", style_of("mc_rsl_off", "label", "element"))
  d$item <- factor(d$item, levels = ord)
  site <- d[d$scope == "site", ]
  net <- d[d$scope == "network", ]
  xl <- if (quantity == "ta_pp") "1σ effect on the target-area fraction [pp]" else "1σ effect on ln x80 [-]"
  # one legend: inputs (top to bottom as on the y axis) as points, then the source references as lines
  brk <- c(rev(ord), src)
  n_in <- length(ord)
  cols <- stats::setNames(c(style_of(ids, "colour", "mc_input"), style_of(src_all, "colour", "element")), c(ids, src_all))
  lts <- stats::setNames(style_of(src_all, "linetype", "element"), src_all)
  g <- ggplot2::ggplot(mapping = ggplot2::aes(y = item)) +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = median, linetype = id, colour = id), linewidth = 0.5,
                        key_glyph = "path") +
    ggplot2::geom_point(data = site, ggplot2::aes(x = median), colour = style_of("mc_site", "colour", "element"),
                        size = 0.9, na.rm = TRUE) +
    ggplot2::geom_errorbar(data = net, ggplot2::aes(xmin = p25, xmax = p75, colour = item), width = 0, linewidth = 0.5,
                           orientation = "y", na.rm = TRUE, show.legend = FALSE) +
    ggplot2::geom_point(data = net, ggplot2::aes(x = median, colour = item), size = 2.2, na.rm = TRUE) +
    ggplot2::facet_grid(daynight ~ regime) +
    ggplot2::scale_colour_manual(values = cols, breaks = brk, labels = c(labs[rev(ord)], style_of(src, "label", "element")),
                                 name = NULL, guide = ggplot2::guide_legend(
                                   ncol = 3, byrow = FALSE,
                                   override.aes = list(linetype = c(rep("blank", n_in), lts[src]),
                                                       shape = c(rep(16, n_in), rep(NA, length(src))),
                                                       size = c(rep(2.2, n_in), rep(NA, length(src))),
                                                       linewidth = c(rep(0, n_in), rep(0.6, length(src)))))) +
    ggplot2::scale_linetype_manual(values = lts, guide = "none") +
    ggplot2::scale_y_discrete(labels = labs) +
    ggplot2::labs(x = xl, y = NULL, subtitle = fp_wrap(subtitle), caption = fp_wrap(caption, 130)) +
    theme_pblh()
  ex <- ggplot2::expansion(mult = c(0.03, 0.05))
  g <- g + if (quantity == "ta_pp") ggplot2::scale_x_continuous(trans = "sqrt", expand = ex) else
    ggplot2::scale_x_continuous(expand = ex)
  g + ggplot2::theme(panel.spacing.x = ggplot2::unit(1.2, "lines"), legend.position = "bottom",
                     legend.key.width = ggplot2::unit(1.6, "lines"), legend.key.height = ggplot2::unit(0.8, "lines"),
                     legend.key.spacing.y = ggplot2::unit(0, "lines"))
}

#' Joint 1-sigma width against the PBLH-source effect, per row
#'
#' @param mc Output of [mc_mc_table()] with `site_id`; @param src Output of [mc_source_effects()] with `site_id`.
#' @param h_option h option [character]; @param floor_pp Smallest value shown on the log axes [pp] (zeros there).
#' @param subtitle,caption Texts.
#' @return A ggplot (log-log, facets by stratum, 1:1 line).
plot_mc_width <- function(mc, src, h_option = "h_pct30", floor_pp = 0.001, subtitle = NULL, caption = NULL) {
  m <- mc[mc$h_option == h_option, ]
  k <- match(paste(m$site_id, m$row_id), paste(src$site_id, src$row_id))
  d <- dplyr::tibble(stratum = m$stratum, width = pmax(m$width_ta_pp, floor_pp),
                     src = pmax(abs(src$src_ta_pp__era5_vs_pblhmod[k]), floor_pp))
  d <- mc_strata_facets(d[!is.na(d$width) & !is.na(d$src), ])
  ggplot2::ggplot(d, ggplot2::aes(width, src)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = style_of("mc_one_to_one", "colour", "element"),
                         linetype = style_of("mc_one_to_one", "linetype", "element")) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6, colour = style_of("mc_network", "colour", "element")) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_grid(daynight ~ regime) +
    ggplot2::labs(x = paste0("Joint 1σ width of the TA fraction [pp] (", h_option, ")"),
                  y = "|ERA5 − PBLHmod| TA fraction [pp]", subtitle = fp_wrap(subtitle), caption = fp_wrap(caption, 130)) +
    theme_pblh()
}

#' Share of perturbed runs that become invalid under the etc set (design A), per site and input
#'
#' @param fs Output of [mc_flip_shares()].
#' @param subtitle,caption Texts.
#' @return A ggplot: one panel per site, bars per input (share of the runs at -1 and +1 sigma invalid, stratum-weighted),
#'   the value written beside each non-zero bar with the share failing the e_rsl test in brackets.
plot_mc_validity_flips <- function(fs, subtitle = NULL, caption = NULL) {
  ids <- style_ids("mc_input", unique(fs$input))
  fs$input <- factor(fs$input, levels = rev(ids))
  pc <- function(p) ifelse(p < 5e-4, "< 0.1", ifelse(p < 0.01, sprintf("%.1f", 100 * p), sprintf("%.0f", 100 * p)))
  fs$txt <- ifelse(fs$share_invalid > 0,
                   paste0(pc(fs$share_invalid), " %",
                          ifelse(fs$share_rsl > 0, paste0(" (", style_of("mc_rsl_flip", "label", "element"), " ",
                                                          pc(fs$share_rsl), " %)"), "")), "")
  ggplot2::ggplot(fs, ggplot2::aes(100 * share_invalid, input)) +
    ggplot2::geom_col(ggplot2::aes(fill = input), width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = txt), hjust = -0.05, size = 2.4,
                       colour = style_of("mc_rsl_flip", "colour", "element")) +
    ggplot2::facet_wrap(~site_id, ncol = 3) +
    ggplot2::scale_fill_manual(values = style_of(ids, "colour", "mc_input"), guide = "none") +
    ggplot2::scale_y_discrete(labels = style_of(ids, "label", "mc_input")) +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "Perturbed runs (−1σ and +1σ) invalid under the etc set [%]", y = NULL,
                  subtitle = fp_wrap(subtitle), caption = fp_wrap(caption, 130)) +
    theme_pblh() + ggplot2::theme(panel.spacing.x = ggplot2::unit(1.2, "lines"))
}
