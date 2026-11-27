# Spec 18 collection: network tables and figures from the per-site tables of scripts/18a_pblh_ec_site.R.

pec_parts <- function() {
  c("variables", "sl_share", "ustar_cross", "nee_residual", "storage_h", "storage_nights", "panofsky", "panofsky_bins",
    "scf_h")
}

#' Bind the per-site spec-18 tables of one part
#'
#' @param dir Folder with `18_<part>_<site>.csv`; @param part Part name ([pec_parts()]).
#' @return Tibble (all sites; empty tibble when no file).
pec_bind_part <- function(dir, part) {
  f <- list.files(dir, pattern = paste0("^18_", part, "_[A-Z]{2}-[A-Za-z0-9]{3}\\.csv$"), full.names = TRUE)
  chr <- c("site_id", "bin", "stratum", "ustar_bin", "ta_tercile", "part_1", "part_2", "part_3", "part_4", "part_5",
           "part_6", "measured_source", "measured_type", "igbp", "ustar_thr_source", "first_utc", "last_utc",
           "labelling_date")
  x <- lapply(f, function(p) {
    d <- readr::read_csv(p, comment = "#", show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
    dplyr::mutate(d, dplyr::across(dplyr::any_of(chr), as.character))
  })
  x <- x[vapply(x, nrow, integer(1)) > 0]
  if (!length(x)) return(dplyr::tibble())
  # a column with only NA in one site's file is read as logical: make it numeric (or character for `chr`)
  allna <- function(v) is.logical(v) && all(is.na(v))
  dplyr::bind_rows(lapply(x, function(d) dplyr::mutate(d, dplyr::across(dplyr::where(allna), as.numeric),
                                                       dplyr::across(dplyr::any_of(chr), as.character))))
}

#' Network rows: median over sites of numeric columns, by grouping columns
#'
#' @param x Site table; @param by Grouping columns; @param cols Columns summarised.
#' @return Tibble: `site_id = "network"`, `by`, `n_sites`, the medians of `cols` (sites with a finite value).
pec_network_median <- function(x, by, cols) {
  if (!nrow(x)) return(dplyr::tibble())
  g <- dplyr::group_by(x, dplyr::across(dplyr::all_of(by)))
  out <- dplyr::summarise(g, n_sites = dplyr::n_distinct(site_id[is.finite(.data[[cols[1]]])]),
                          dplyr::across(dplyr::all_of(cols), ~ med_fin(.x)), .groups = "drop")
  dplyr::bind_cols(dplyr::tibble(site_id = rep("network", nrow(out))), out)
}

#' Sites whose stable-night share above the surface layer exceeds a limit for any source
#'
#' @param sl `18_sl_share` (site rows); @param limit Share [-]; @param threshold zm/h threshold column.
#' @return Tibble: `site_id`, `sources` (`;`-separated source ids above the limit, common rows), `max_share`.
pec_stable_night_flags <- function(sl, limit = 0.25, threshold = "share_gt_0.1") {
  x <- sl[sl$row_basis == "common" & sl$daynight == "night" & sl$zl == "stable" & sl$season == "all" &
            sl$site_id != "network", ]
  x <- x[is.finite(x[[threshold]]) & x[[threshold]] > limit, ]
  if (!nrow(x)) return(dplyr::tibble(site_id = character(), sources = character(), max_share = numeric()))
  dplyr::summarise(dplyr::group_by(x, site_id), sources = paste(source, collapse = ";"),
                   max_share = max(.data[[threshold]]), .groups = "drop")
}

#' How often each source ranks first in the Panofsky test
#'
#' @param pan `18_panofsky` (site rows).
#' @return Tibble: `row_basis`, `source`, `n_sites` (sites ranked), `n_first`, `share_first` [-].
pec_panofsky_rank <- function(pan) {
  x <- pan[pan$row_basis != "own" & !is.na(pan$rank) & pan$site_id != "network", ]
  if (!nrow(x)) return(dplyr::tibble())
  dplyr::summarise(dplyr::group_by(x, row_basis, source), n_sites = dplyr::n(), n_first = sum(rank == 1),
                   share_first = n_first / n_sites, .groups = "drop")
}

pec_circ_note <- function() {
  "PBLHmod's night-time h is computed from u* and L (spec 02), so with PBLHmod zm/h > 0.1 is close to low u* (circular); ERA5 is independent of the station's turbulence."
}

#' Figure 18_sl_share: stable-night share above the surface layer per site, three sources side by side
#'
#' @param sl `18_sl_share`; @param vars `18_variables` (`site_id`, `zm_median_m`); @param caption Caption.
#' @return ggplot.
plot_pec_sl_share <- function(sl, vars, caption = NULL) {
  srcs <- pec_sources()
  x <- sl[sl$row_basis == "common" & sl$daynight == "night" & sl$zl == "stable" & sl$season == "all" &
            sl$site_id != "network" & sl$source %in% srcs, ]
  zm <- stats::setNames(vars$zm_median_m, vars$site_id)
  ord <- names(sort(zm[unique(x$site_id)]))
  x$site <- factor(paste0(x$site_id, " (", formatC(zm[x$site_id], format = "f", digits = 1), " m)"),
                   levels = paste0(ord, " (", formatC(zm[ord], format = "f", digits = 1), " m)"))
  x$source <- factor(x$source, levels = style_ids("source", srcs))
  ggplot2::ggplot(x, ggplot2::aes(share_gt_0.1, site, colour = source, shape = source)) +
    ggplot2::geom_vline(xintercept = 0.25, colour = style_of("reference_line", "colour", "element"),
                        linetype = "dashed") +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 1.8) +
    ggplot2::scale_colour_manual(values = stats::setNames(style_of(levels(x$source)), levels(x$source)),
                                 labels = stats::setNames(style_of(levels(x$source), "label"), levels(x$source)), name = NULL) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 1, 17)[seq_along(levels(x$source))], levels(x$source)),
                                labels = stats::setNames(style_of(levels(x$source), "label"), levels(x$source)), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "share of stable-night half-hours with zm / h > 0.1 [-]", y = "site (zm, ordered by zm)",
                  caption = caption) +
    theme_pblh()
}

#' Figure 18_nee_residual: median NEE residual difference (zm/h > 0.1 minus <= 0.1) by u* bin
#'
#' @param nr `18_nee_residual`; @param order_sites Site order (top to bottom); @param flux,measure Rows shown.
#' @param caption Caption.
#' @return ggplot.
plot_pec_nee_residual <- function(nr, order_sites, flux = "NEE", measure = "abs", caption = NULL) {
  x <- nr[nr$flux == flux & nr$measure == measure & nr$stratum != "cell" & nr$site_id != "network", ]
  x$ustar_bin[x$stratum == "all"] <- "all"
  x$ustar_bin <- factor(x$ustar_bin, levels = c("1-1.5", "1.5-2", ">2", "all"))
  x$source <- factor(x$source, levels = style_ids("source", unique(x$source)))
  x <- x[is.finite(x$diff), ]
  x$site_id <- factor(x$site_id, levels = intersect(order_sites, unique(x$site_id)))
  pd <- ggplot2::position_dodge(width = 0.5)
  ggplot2::ggplot(x, ggplot2::aes(ustar_bin, diff, colour = source)) +
    ggplot2::geom_hline(yintercept = 0, colour = style_of("reference_line", "colour", "element")) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lo90, ymax = hi90), position = pd, size = 0.25) +
    ggplot2::scale_colour_manual(values = stats::setNames(style_of(levels(x$source)), levels(x$source)),
                                 labels = stats::setNames(style_of(levels(x$source), "label"), levels(x$source)),
                                 name = "h source", drop = FALSE) +
    ggplot2::facet_wrap(~site_id) +
    ggplot2::labs(x = "u* / u*_thr bin",
                  y = if (measure == "abs") paste0("median ", flux, " - RECO_NT,\nabove - within [µmol m-2 s-1]")
                  else paste0("median (", flux, " - RECO_NT) / RECO_NT,\nabove - within [-]"), caption = caption) +
    theme_pblh()
}

#' Figure 18_storage_h: night storage accumulation against the night median h
#'
#' @param nights `18_storage_nights`; @param st `18_storage_h` (binned medians); @param caption Caption.
#' @return ggplot.
plot_pec_storage_h <- function(nights, st, caption = NULL) {
  srcs <- c("era5_blh_center", "pblhmod_default")
  long <- dplyr::bind_rows(lapply(srcs, function(s) {
    dplyr::tibble(site_id = nights$site_id, source = s, h = nights[[paste0("h_night_", s)]], accum = nights$accum)
  }))
  long <- long[is.finite(long$h) & long$h > 0 & is.finite(long$accum), ]
  long$source <- factor(long$source, levels = style_ids("source", srcs))
  b <- st[st$part == "accumulation_bin" & is.finite(st$value), ]
  if (nrow(b)) {
    lo <- as.numeric(sub("^\\[([^,]+),.*$", "\\1", b$stratum))
    hi <- as.numeric(sub("^.*,([^)]+)\\)$", "\\1", b$stratum))
    b$h_mid <- ifelse(is.finite(hi), sqrt(pmax(lo, 1) * hi), lo * 1.5)
    b$source <- factor(b$source, levels = levels(long$source))
  }
  ggplot2::ggplot(long, ggplot2::aes(h, accum / 1e3, colour = source)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.25) +
    { if (nrow(b)) ggplot2::geom_line(data = b, ggplot2::aes(h_mid, value / 1e3, colour = source), linewidth = 0.8) } +
    ggplot2::scale_colour_manual(values = stats::setNames(style_of(levels(long$source)), levels(long$source)),
                                 labels = stats::setNames(style_of(levels(long$source), "label"), levels(long$source)),
                                 name = "h source") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~site_id, scales = "free_y") +
    ggplot2::labs(x = "night median h [m]", y = "night CO2 storage accumulation [mmol m-2]", caption = caption) +
    theme_pblh()
}

#' Figure 18_panofsky: binned sigma_v / u* against -h/L per source, with the Panofsky et al. (1977) curve
#'
#' @param bins `18_panofsky_bins`; @param sites Sites shown; @param caption Caption; @param min_n Bins with fewer rows
#'   not drawn.
#' @return ggplot.
plot_pec_panofsky <- function(bins, sites, caption = NULL, min_n = 20) {
  x <- bins[bins$site_id %in% sites & bins$n >= min_n & is.finite(bins$sv_ustar_median), ]
  x$site_id <- factor(x$site_id, levels = intersect(sites, unique(x$site_id)))
  x$source <- factor(x$source, levels = style_ids("source", unique(x$source)))
  xr <- range(x$x_mid, na.rm = TRUE)
  cv <- dplyr::tibble(x_mid = exp(seq(log(max(xr[1], 0.1)), log(xr[2]), length.out = 60)))
  cv$y <- (12 + 0.5 * cv$x_mid)^(1 / 3)
  ggplot2::ggplot(x, ggplot2::aes(x_mid, sv_ustar_median, colour = source)) +
    ggplot2::geom_line(data = cv, ggplot2::aes(x_mid, y), inherit.aes = FALSE,
                       colour = style_of("reference_line", "colour", "element"), linetype = "dashed") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = stats::setNames(style_of(levels(x$source)), levels(x$source)),
                                 labels = stats::setNames(style_of(levels(x$source), "label"), levels(x$source)),
                                 name = "h source") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~site_id) +
    ggplot2::labs(x = "-h / L [-]", y = "median σv / u* [-]", caption = caption) +
    theme_pblh()
}

#' Figure 18s_scf_h: spectral correction factors against zm/h (unstable rows), per site and network median
#'
#' @param scf `18_scf_h` (site and network rows); @param source h source shown; @param caption Caption.
#' @return ggplot.
plot_pec_scf <- function(scf, source = "era5_blh_center", caption = NULL, min_n = 20) {
  x <- scf[scf$source == source & is.finite(scf$median) & is.finite(scf$x_mid) & scf$x_mid > 0, ]
  x <- x[x$site_id == "network" | x$n >= min_n, ]
  x$x_var <- factor(ifelse(x$x_var == "zm_h", "zm / h [-]", "-h / L [-]"), levels = c("zm / h [-]", "-h / L [-]"))
  site <- x[x$site_id != "network", ]
  net <- x[x$site_id == "network", ]
  ggplot2::ggplot(site, ggplot2::aes(x_mid, median, group = site_id)) +
    ggplot2::geom_line(colour = style_of("context", "colour", "element"), linewidth = 0.3, alpha = 0.7) +
    ggplot2::geom_line(data = net, ggplot2::aes(group = 1), colour = style_of(source), linewidth = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_grid(scf_var ~ x_var, scales = "free") +
    ggplot2::labs(x = NULL, y = "median total spectral correction factor [-]", caption = caption) +
    theme_pblh()
}
