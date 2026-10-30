# Plot helpers: figure style (one table and one theme), saving (png + pdf), faded-context helpers.
# Rules: CLAUDE.md (Style), docs/specs/04_era5_standard_compare.md section 8, spec 09 Docs (figure style, G 2026-10-30).
#
# Figure style (G 2026-10-30): every figure takes its colours, labels, order and line types only from
# `pblh_sources_style` and its theme only from theme_pblh(), so that all figures can be restyled by editing these
# two. The colours below are placeholders (colour-blind-safe Okabe & Ito 2008 set where possible); G will replace
# them with the final theme.

# Placeholder colours, each defined once (edit here to restyle every figure).
.pblh_col <- c(
  blue = "#0072B2", green = "#009E73", vermillion = "#D55E00", orange = "#E69F00", sky = "#56B4E9",
  purple = "#CC79A7", yellow = "#F0E442", black = "#000000", grey = "#999999", grey_light = "#DDDDDD",
  grey_dark = "grey45", blue_light = "#6BAED6", blue_dark = "#08519C"
)

#' Style table of every plotted category
#'
#' One row per category: PBLH sources (standard source ids), footprint cases, regimes, footprint subsets, qc classes
#' of the measured MLH, legacy-regression variants and figure elements (bands, trend lines, highlighted classes).
#'
#' @format Tibble: `id`, `group` (`source`, `footprint_case`, `regime`, `subset`, `qc`, `legacy`, `element`),
#'   `colour` [hex or R colour name], `label`, `order` [integer, legend and plotting order within the group],
#'   `linetype` (ggplot2 line type).
pblh_sources_style <- local({
  k <- .pblh_col
  row <- function(id, group, colour, label, order, linetype = "solid") {
    dplyr::tibble(id = id, group = group, colour = unname(colour), label = label, order = as.integer(order),
                  linetype = linetype)
  }
  dplyr::bind_rows(
    # PBLH sources (main figures: PBLHmod, PBLHmod-clim, ERA5, measured; supplement: reference gamma values, oper)
    row("pblhmod_default", "source", k["blue"], "PBLHmod", 1),
    row("pblhmod_clim", "source", k["green"], "PBLHmod-clim", 2),
    row("era5_blh_center", "source", k["vermillion"], "ERA5 (tower cell)", 3),
    row("era5_blh_bilinear", "source", k["vermillion"], "ERA5 (bilinear)", 4, "dashed"),
    row("cloudnet_eps_mlh", "source", k["black"], "Cloudnet ε MLH", 5),
    row("sirta_mlh_chm15k", "source", k["black"], "SIRTA CHM15k (STRATfinder)", 5),
    row("meteoswiss_stratfinder_mlh", "source", k["black"], "MeteoSwiss STRATfinder", 5),
    row("sirta_mlh_cl31", "source", k["grey"], "SIRTA CL31 (CABAM)", 6),
    row("cloudnet_censored", "source", k["black"], "Measured MLH (censored median)", 5),
    row("cloudnet_exact", "source", k["grey"], "Measured MLH (qc 0/2 median)", 6, "dashed"),
    row("pblhmod_gamma001", "source", k["sky"], "PBLHmod (γ 0.01)", 7),
    row("pblhmod_gamma0025", "source", k["orange"], "PBLHmod (γ 0.0025)", 8),
    row("icos_oper_qc0", "source", k["purple"], "ICOS operational (qc 0)", 9),
    row("icos_oper_all", "source", k["purple"], "ICOS operational (all)", 10, "dashed"),
    # footprint cases (spec 08 / 09 A5)
    row("pblhmod", "footprint_case", k["blue"], "PBLHmod", 1),
    row("pblhmod_clim", "footprint_case", k["green"], "PBLHmod-clim", 2),
    row("era5", "footprint_case", k["vermillion"], "ERA5 (tower cell)", 3),
    row("era5_bilinear", "footprint_case", k["vermillion"], "ERA5 (bilinear)", 4, "dashed"),
    row("measured", "footprint_case", k["black"], "Measured MLH", 5),
    row("bound", "footprint_case", k["grey"], "Measured MLH (lower bound)", 6, "dashed"),
    row("pblhmod_x0.8", "footprint_case", k["blue_light"], "PBLHmod × 0.8", 7, "dotted"),
    row("pblhmod_x1.2", "footprint_case", k["blue_dark"], "PBLHmod × 1.2", 8, "dashed"),
    row("pblhmod_gamma001", "footprint_case", k["sky"], "PBLHmod (γ 0.01)", 9),
    row("pblhmod_gamma0025", "footprint_case", k["orange"], "PBLHmod (γ 0.0025)", 10),
    row("oper", "footprint_case", k["purple"], "ICOS operational", 11),
    # regimes (labels also used for facets)
    row("all", "regime", k["black"], "all", 0),
    row("stable", "regime", k["blue"], "stable", 1),
    row("unstable_diag", "regime", k["orange"], "unstable (diagnostic)", 2),
    row("convective", "regime", k["vermillion"], "convective (B&G)", 3),
    row("mixed", "regime", k["grey"], "mixed (transition hour)", 4),
    # footprint subsets
    row("valid", "subset", k["black"], "valid", 1),
    row("valid_qc", "subset", k["blue"], "valid, NEE QC good", 2),
    row("valid_qc_used", "subset", k["green"], "valid, NEE used (u* filter)", 3),
    row("valid_qc_h", "subset", k["orange"], "valid, H measured", 4),
    # qc classes of the Cloudnet MLH (spec 07)
    row("0", "qc", k["blue"], "0 measured", 1),
    row("2", "qc", k["sky"], "2 partly censored", 2),
    row("3", "qc", k["orange"], "3 signal-limited (lower bound)", 3),
    row("1", "qc", k["grey"], "1 censored below z_min", 4),
    row("NA", "qc", k["grey_light"], "NA coverage < min", 5),
    # legacy regression (spec 03)
    row("operational", "legacy", k["purple"], "Operational L2", 1),
    row("legacy_baseline", "legacy", k["orange"], "Legacy code (L2 inputs)", 2, "22"),
    row("legacy_no_postproc", "legacy", k["orange"], "Legacy without fill / floor / smoothing", 3),
    row("pblhmod_v2", "legacy", k["blue"], "PBLHmod v2", 4),
    # figure elements
    row("aerosol_window", "element", k["yellow"], "aerosol window", 1),
    row("trend", "element", k["vermillion"], "trend", 2),
    row("context", "element", k["grey_dark"], "not comparable (shown for context)", 3),
    row("reference_line", "element", "grey50", "reference line", 4),
    row("box_fill", "element", "grey85", "box", 5),
    row("target_area", "element", "grey40", "ICOS target area", 6, "dashed"),
    row("top_all", "element", "grey60", "all", 7),
    row("top25", "element", k["orange"], "top 25 %", 8),
    row("top10", "element", k["vermillion"], "top 10 %", 9),
    row("gamma_era5", "element", k["green"], "ERA5 γ", 10),
    row("hour_00", "element", k["blue"], "00 UTC", 11),
    row("hour_12", "element", k["vermillion"], "12 UTC", 12),
    row("point_fill", "element", "white", "open point", 13)
  )
})

#' Style values of plotted categories
#'
#' @param ids Category ids [character].
#' @param what Column of [pblh_sources_style]: `"colour"`, `"label"`, `"linetype"` or `"order"`.
#' @param group Group of the ids (`"source"`, `"footprint_case"`, ...).
#' @return Named vector (names = `ids`); aborts on an id without a style row.
style_of <- function(ids, what = "colour", group = "source") {
  st <- pblh_sources_style[pblh_sources_style$group == group, ]
  miss <- setdiff(ids, st$id)
  if (length(miss)) cli::cli_abort("No style for {.val {miss}} in group {.val {group}} ({.code pblh_sources_style}).")
  stats::setNames(st[[what]][match(ids, st$id)], ids)
}

#' Ids of a style group in plotting order
#'
#' @param group Style group.
#' @param ids Optional subset of ids to keep [character].
#' @return Character vector of ids, ordered by `order`.
style_ids <- function(group = "source", ids = NULL) {
  st <- pblh_sources_style[pblh_sources_style$group == group, ]
  st <- st[order(st$order, seq_len(nrow(st))), ]
  if (is.null(ids)) st$id else st$id[st$id %in% ids]
}

#' Theme of all figures
#'
#' @param base_size Base font size [pt].
#' @return A ggplot2 theme.
theme_pblh <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1.2, "lines"),
      panel.spacing = ggplot2::unit(0.6, "lines")
    )
}

#' Colours of the PBLH sources
#'
#' @return Named character vector of colours, names = source ids ([pblh_sources_style], group `source`).
source_palette <- function() style_of(style_ids("source"), "colour")

#' Readable source labels
#'
#' @return Named character vector, names = source ids ([pblh_sources_style], group `source`).
source_labels <- function() style_of(style_ids("source"), "label")

#' Save a figure (png by default; pdf only for paper figures)
#'
#' A list of ggplot objects is laid out side by side (base `grid` viewports;
#' no extra package). PNG at 300 dpi by default; add `"pdf"` to `formats`
#' only for the final paper figures (CLAUDE.md, Style).
#'
#' @param p A ggplot object or a list of them.
#' @param name File name without extension [character].
#' @param width_cm,height_cm Size [cm].
#' @param dir Output directory; default `<output_dir>/figures` from `config.yml`.
#' @param formats Output formats, any of `"png"`, `"pdf"` [character];
#'   default `"png"`.
#' @return Named character vector of the written paths (names = formats),
#'   invisibly.
save_figure <- function(p, name, width_cm = 16, height_cm = 10, dir = NULL, formats = "png") {
  bad <- setdiff(formats, c("png", "pdf"))
  if (length(formats) == 0 || length(bad) > 0) cli::cli_abort("{.arg formats} must be any of {.val png}, {.val pdf}.")
  if (is.null(dir)) dir <- file.path(resolve_path(get_cfg()$output_dir, must_exist = FALSE), "figures")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  plots <- if (inherits(p, "ggplot")) list(p) else p
  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(plots))))
    for (i in seq_along(plots)) {
      print(plots[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i), newpage = FALSE)
    }
    invisible(NULL)
  }
  w <- width_cm / 2.54
  h <- height_cm / 2.54
  paths <- c()
  if ("png" %in% formats) {
    paths["png"] <- file.path(dir, paste0(name, ".png"))
    grDevices::png(paths[["png"]], width = w, height = h, units = "in", res = 300)
    tryCatch(draw(), finally = grDevices::dev.off())
  }
  if ("pdf" %in% formats) {
    paths["pdf"] <- file.path(dir, paste0(name, ".pdf"))
    pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    pdf_dev(paths[["pdf"]], width = w, height = h)
    tryCatch(draw(), finally = grDevices::dev.off())
  }
  invisible(paths)
}

# --- Spec 11 E: plotting convention (G 2026-10-29). In every comparison figure, values outside the comparable set
# (distance class, aerosol window, regime) are drawn but faded (alpha 0.3) behind the comparable ones, with the
# legend entry "not comparable (shown for context)". Metrics are computed on the comparable set only. When every
# value is comparable the helpers add nothing, so the figure is unchanged.

#' Legend label of the faded context values
#'
#' @return Character.
context_label <- function() "not comparable (shown for context)"

#' Comparable / context factor
#'
#' @param comparable Logical per row (`NA` counts as not comparable).
#' @return Factor with levels `comparable`, `context`.
comparable_set <- function(comparable) {
  factor(ifelse(comparable %in% TRUE, "comparable", "context"), levels = c("comparable", "context"))
}

#' Alpha scale of the comparable set (legend shows the context entry only)
#'
#' @param context_alpha Opacity of the context values [0-1].
#' @return ggplot2 scale.
scale_alpha_comparable <- function(context_alpha = 0.3) {
  ggplot2::scale_alpha_manual(values = c(comparable = 1, context = context_alpha), breaks = "context",
                              labels = c(context = context_label()), name = NULL, drop = TRUE)
}

#' Faded points for the context rows
#'
#' @param data Data frame with the plotted columns and a logical `comparable` column.
#' @param mapping `ggplot2::aes()` of x and y.
#' @param size Point size [mm].
#' @return List of layers (empty when every row is comparable).
geom_context_points <- function(data, mapping, size = 0.5) {
  d <- data[!(data$comparable %in% TRUE), ]
  if (nrow(d) == 0) return(list())
  d$.set <- comparable_set(d$comparable)
  list(ggplot2::geom_point(data = d, mapping = utils::modifyList(mapping, ggplot2::aes(alpha = .data$.set)),
                           colour = style_of("context", "colour", "element"), size = size, stroke = 0, inherit.aes = FALSE),
       scale_alpha_comparable())
}

#' Lines faded outside the comparable set
#'
#' The whole line is drawn faded and its comparable stretches on top at full opacity (each stretch a separate group,
#' so no segment joins two stretches across a context gap).
#'
#' @param data Data frame with the x, y, colour columns, the grouping columns and a logical `comparable` column.
#' @param mapping `ggplot2::aes()` of x, y and colour (and linetype).
#' @param group_cols Columns defining one line (e.g. `c("season", "variable")`) [character].
#' @param x Name of the x column (order along the line) [character].
#' @param linewidth Line width [mm].
#' @return List of layers: a plain `geom_line()` when every row is comparable.
geom_line_comparable <- function(data, mapping, group_cols, x, linewidth = 0.6) {
  if (all(data$comparable %in% TRUE)) return(list(ggplot2::geom_line(data = data, mapping = mapping,
                                                                     linewidth = linewidth, na.rm = TRUE)))
  data <- data[do.call(order, c(unname(as.list(data[group_cols])), list(data[[x]]))), ]
  key <- interaction(data[group_cols], drop = TRUE, lex.order = TRUE)
  run <- stats::ave(as.integer(data$comparable %in% TRUE), key, FUN = function(v) cumsum(c(1L, diff(v) != 0)))
  data$.line <- interaction(key, drop = TRUE)
  data$.seg <- interaction(key, run, drop = TRUE)
  data$.set_ctx <- comparable_set(FALSE)
  com <- data[data$comparable %in% TRUE, ]
  com$.set <- comparable_set(TRUE)
  list(
    ggplot2::geom_line(data = data, mapping = utils::modifyList(mapping, ggplot2::aes(group = .data$.line,
                                                                                    alpha = .data$.set_ctx)),
                       linewidth = linewidth, na.rm = TRUE),
    ggplot2::geom_line(data = com, mapping = utils::modifyList(mapping, ggplot2::aes(group = .data$.seg,
                                                                                   alpha = .data$.set)),
                       linewidth = linewidth, na.rm = TRUE),
    scale_alpha_comparable()
  )
}
