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
#' for all figures.
#'
#' @return Named character vector of hex colours, names = source ids.
source_palette <- function() {
  c(
    pblhmod_default = "#0072B2",   # blue
    pblhmod_gamma0025 = "#56B4E9", # sky blue
    era5_blh_center = "#D55E00",   # vermillion
    era5_blh_bilinear = "#E69F00", # orange
    icos_oper_qc0 = "#009E73",     # bluish green
    icos_oper_all = "#009E73"      # same source, dashed line
  )
}

#' Readable source labels
#'
#' @return Named character vector, names = source ids.
source_labels <- function() {
  c(
    pblhmod_default = "PBLHmod (default)",
    pblhmod_gamma0025 = "PBLHmod (gamma 0.0025)",
    era5_blh_center = "ERA5 (tower cell)",
    era5_blh_bilinear = "ERA5 (bilinear)",
    icos_oper_qc0 = "ICOS operational (qc 0)",
    icos_oper_all = "ICOS operational (all)"
  )
}

#' Save a figure as png (300 dpi) and pdf
#'
#' A list of ggplot objects is laid out side by side (base `grid` viewports;
#' no extra package).
#'
#' @param p A ggplot object or a list of them.
#' @param name File name without extension [character].
#' @param width_cm,height_cm Size [cm].
#' @param dir Output directory; default `<output_dir>/figures` from `config.yml`.
#' @return Paths of the two files, invisibly.
save_figure <- function(p, name, width_cm = 16, height_cm = 10, dir = NULL) {
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
  png_path <- file.path(dir, paste0(name, ".png"))
  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  grDevices::png(png_path, width = w, height = h, units = "in", res = 300)
  tryCatch(draw(), finally = grDevices::dev.off())
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  pdf_dev(pdf_path, width = w, height = h)
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(c(png = png_path, pdf = pdf_path))
}
