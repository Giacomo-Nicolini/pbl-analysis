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
#'   amendment 1), labelled "rsl test off"; for these (and any input with a flip) the network share of runs that fail
#'   e_rsl is written beside the points, so that a flip never reads as a zero effect.
#' @return A ggplot: inputs on y (ordered by the network median over all strata, largest on top), the 1-sigma effect
#'   on x; sites as small points, network as the large point with the IQR over sites; vertical reference lines at the
#'   network median |ERA5 - PBLHmod| (and |PBLHmod-clim - PBLHmod|) of the stratum.
plot_mc_input_ranking <- function(s, quantity = "ta_pp", subtitle = NULL, caption = NULL, rsl_off = c("hc", "z0f")) {
  d <- s[s$kind == "oat_1sigma" & s$quantity == quantity, ]
  ro <- paste0(rsl_off, "|rsl_off")
  have <- rsl_off[ro %in% d$item]
  d <- d[!d$item %in% have, ] # the etc rows of these inputs replaced by their rsl_off rows
  d$item[d$item %in% paste0(have, "|rsl_off")] <- sub("\\|rsl_off$", "", d$item[d$item %in% paste0(have, "|rsl_off")])
  net_all <- d[d$scope == "network" & d$stratum == "all", ]
  d <- mc_strata_facets(d)
  ref <- mc_strata_facets(s[s$kind == "source_effect" & s$quantity == quantity & s$scope == "network", ])
  ref$id <- ifelse(ref$item == "era5_vs_pblhmod", "mc_source_era5", "mc_source_clim")
  ref <- ref[!is.na(ref$median), ]
  ids <- style_ids("mc_input", unique(d$item))
  ord <- ids[order(net_all$median[match(ids, net_all$item)], na.last = FALSE)]
  labs <- style_of(ord, "label", "mc_input")
  labs[ord %in% have] <- paste0(labs[ord %in% have], ", ", style_of("mc_rsl_off", "label", "element"))
  d$item <- factor(d$item, levels = ord)
  site <- d[d$scope == "site", ]
  net <- d[d$scope == "network", ]
  fl <- net[!is.na(net$share_rsl_flip) & (net$item %in% have | net$share_rsl_flip > 0), ]
  fl$txt <- paste0(style_of("mc_rsl_flip", "label", "element"), " ", round(100 * fl$share_rsl_flip), " %")
  xl <- if (quantity == "ta_pp") "1σ effect on the target-area fraction [pp]" else "1σ effect on ln x80 [-]"
  g <- ggplot2::ggplot(mapping = ggplot2::aes(y = item)) +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = median, linetype = id, colour = id), linewidth = 0.5) +
    ggplot2::geom_point(data = site, ggplot2::aes(x = median), colour = style_of("mc_site", "colour", "element"),
                        size = 0.9, na.rm = TRUE) +
    ggplot2::geom_errorbar(data = net, ggplot2::aes(xmin = p25, xmax = p75, colour = item), width = 0, linewidth = 0.5,
                           orientation = "y", na.rm = TRUE) +
    ggplot2::geom_point(data = net, ggplot2::aes(x = median, colour = item), size = 2.2, na.rm = TRUE) +
    ggplot2::geom_text(data = fl, ggplot2::aes(x = Inf, label = txt), hjust = 1.03, vjust = 0.5, size = 2.5,
                       colour = style_of("mc_rsl_flip", "colour", "element")) +
    ggplot2::facet_grid(daynight ~ regime) +
    ggplot2::scale_colour_manual(values = c(style_of(ids, "colour", "mc_input"),
                                            style_of(c("mc_source_era5", "mc_source_clim"), "colour", "element")),
                                 breaks = c("mc_source_era5", "mc_source_clim"),
                                 labels = style_of(c("mc_source_era5", "mc_source_clim"), "label", "element"), name = NULL) +
    ggplot2::scale_linetype_manual(values = style_of(c("mc_source_era5", "mc_source_clim"), "linetype", "element"),
                                   labels = style_of(c("mc_source_era5", "mc_source_clim"), "label", "element"), name = NULL) +
    ggplot2::scale_y_discrete(labels = labs) +
    ggplot2::labs(x = xl, y = NULL, subtitle = fp_wrap(subtitle), caption = fp_wrap(caption, 130)) +
    theme_pblh()
  # right margin of each panel kept free for the flip shares (written on the input's row, never over a point)
  ex <- ggplot2::expansion(mult = c(0.03, if (nrow(fl)) 0.6 else 0.05))
  g <- g + if (quantity == "ta_pp") ggplot2::scale_x_continuous(trans = "sqrt", expand = ex) else
    ggplot2::scale_x_continuous(expand = ex)
  g <- g + ggplot2::theme(panel.spacing.x = ggplot2::unit(1.2, "lines"))
  g
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

#' Share of perturbed runs that become invalid (design A), by input and stratum
#'
#' @param oat Bound [mc_oat_table()] outputs (all sites).
#' @param subtitle,caption Texts.
#' @return A ggplot (bars per input, facets by stratum).
plot_mc_validity_flips <- function(oat, subtitle = NULL, caption = NULL) {
  o <- oat[oat$input != "base", ]
  if ("validity" %in% names(o)) o <- o[o$validity == "etc", ] # flips under the etc set (rsl_off runs excluded)
  d <- dplyr::summarise(dplyr::group_by(o, stratum, input), share = mean(!valid_etc), n = dplyr::n(), .groups = "drop")
  d <- mc_strata_facets(d)
  ids <- style_ids("mc_input", unique(d$input))
  d$input <- factor(d$input, levels = rev(ids))
  ggplot2::ggplot(d, ggplot2::aes(100 * share, input, fill = input)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_grid(daynight ~ regime) +
    ggplot2::scale_fill_manual(values = style_of(ids, "colour", "mc_input"), guide = "none") +
    ggplot2::scale_y_discrete(labels = style_of(ids, "label", "mc_input")) +
    ggplot2::labs(x = "Perturbed runs (−1σ and +1σ) invalid under the etc set [%]", y = NULL,
                  subtitle = fp_wrap(subtitle), caption = fp_wrap(caption, 130)) +
    theme_pblh()
}
