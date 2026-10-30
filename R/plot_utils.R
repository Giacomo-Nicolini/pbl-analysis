# Plot helpers: publication theme, source palette, saving (png + pdf).
# Rules: CLAUDE.md (Style) and docs/specs/04_era5_standard_compare.md section 8.

#' Publication theme
#'
#' @param base_size Base font size [pt].
#' @return A ggplot2 theme.
theme_pub <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1.2, "lines")
    )
}

#' Colours of the PBLH sources
#'
#' One colour-blind-safe palette (Okabe & Ito 2008, "Color Universal Design")
#' for all figures. Main figures show PBLHmod, PBLHmod-clim, ERA5 and the
#' measurement; the reference gamma values (`pblhmod_gamma001`,
#' `pblhmod_gamma0025`) only in supplementary figures (G 2026-10-30), so
#' `pblhmod_gamma001` shares the colour of `pblhmod_default` (equal for now).
#'
#' @return Named character vector of hex colours, names = source ids.
source_palette <- function() {
  c(
    pblhmod_default = "#0072B2",   # blue
    pblhmod_clim = "#009E73",      # bluish green
    pblhmod_gamma001 = "#0072B2",  # blue (supplementary figures only)
    pblhmod_gamma0025 = "#56B4E9", # sky blue
    era5_blh_center = "#D55E00",   # vermillion
    era5_blh_bilinear = "#E69F00", # orange
    icos_oper_qc0 = "#CC79A7",     # reddish purple (internal, --with-oper)
    icos_oper_all = "#CC79A7"      # same source, dashed line
  )
}

#' Readable source labels
#'
#' @return Named character vector, names = source ids.
source_labels <- function() {
  c(
    pblhmod_default = "PBLHmod",
    pblhmod_clim = "PBLHmod-clim",
    pblhmod_gamma001 = "PBLHmod (γ 0.01)",
    pblhmod_gamma0025 = "PBLHmod (γ 0.0025)",
    era5_blh_center = "ERA5 (tower cell)",
    era5_blh_bilinear = "ERA5 (bilinear)",
    icos_oper_qc0 = "ICOS operational (qc 0)",
    icos_oper_all = "ICOS operational (all)"
  )
}

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
                           colour = "grey45", size = size, stroke = 0, inherit.aes = FALSE),
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
