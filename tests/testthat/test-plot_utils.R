test_that("comparable_set: NA and FALSE are context", {
  expect_equal(as.character(comparable_set(c(TRUE, FALSE, NA))), c("comparable", "context", "context"))
  expect_equal(levels(comparable_set(TRUE)), c("comparable", "context"))
})

test_that("geom_context_points: only the context rows, faded; nothing when all comparable", {
  d <- data.frame(x = 1:6, y = 1:6, comparable = c(TRUE, TRUE, FALSE, NA, TRUE, FALSE))
  expect_length(geom_context_points(transform(d, comparable = TRUE), ggplot2::aes(x, y)), 0)
  l <- geom_context_points(d, ggplot2::aes(x, y))
  expect_equal(l[[1]]$data$x, c(3, 4, 6))
  p <- ggplot2::ggplot(d[d$comparable %in% TRUE, ], ggplot2::aes(x, y)) + l + ggplot2::geom_point()
  b <- ggplot2::ggplot_build(p)
  expect_equal(unique(b$data[[1]]$alpha), 0.3)
  expect_equal(nrow(b$data[[2]]), 3) # comparable points at full opacity
})

test_that("geom_line_comparable: faded full line, comparable stretches on top as separate groups", {
  d <- data.frame(h = 0:7, v = c(1, 2, 3, 4, 5, 6, 7, 8), var = "a", s = "JJA",
                  comparable = c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))
  plain <- geom_line_comparable(transform(d, comparable = TRUE), ggplot2::aes(h, v, colour = var), c("s", "var"), "h")
  expect_length(plain, 1)
  l <- geom_line_comparable(d, ggplot2::aes(h, v, colour = var), c("s", "var"), "h")
  expect_length(l, 3)
  b <- ggplot2::ggplot_build(ggplot2::ggplot() + l)
  expect_equal(nrow(b$data[[1]]), 8)
  expect_equal(unique(b$data[[1]]$alpha), 0.3)
  expect_equal(b$data[[2]]$x, c(1, 2, 5, 6))
  expect_equal(length(unique(b$data[[2]]$group)), 2) # hours 1-2 and 5-6 are not joined
  expect_equal(unique(b$data[[2]]$alpha), 1)
})

test_that("pblh_sources_style: complete, unique ids per group, valid colours and line types (spec 09 figure style)", {
  st <- pblh_sources_style
  expect_false(anyDuplicated(paste(st$group, st$id)) > 0)
  expect_true(all(!is.na(st$colour) & !is.na(st$label) & !is.na(st$order)))
  expect_silent(grDevices::col2rgb(st$colour))
  expect_true(all(st$linetype %in% c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "22")))
  # every standard source id and every footprint case of the pipeline has a style row
  src <- c("pblhmod_default", "pblhmod_clim", "pblhmod_gamma001", "pblhmod_gamma0025", "era5_blh_center",
           "era5_blh_bilinear", "icos_oper_qc0", "icos_oper_all", "cloudnet_eps_mlh", "sirta_mlh_chm15k",
           "sirta_mlh_cl31", "meteoswiss_stratfinder_mlh")
  expect_silent(style_of(src))
  expect_silent(style_of(c("pblhmod", "pblhmod_clim", "era5", "era5_bilinear", "measured", "bound", "pblhmod_x0.8",
                           "pblhmod_x1.2", "pblhmod_gamma001", "pblhmod_gamma0025", "oper"), "colour", "footprint_case"))
  expect_equal(unname(style_of(c("pblhmod_default", "pblhmod_clim"), "label")), c("PBLHmod", "PBLHmod-clim"))
  expect_equal(style_ids("source", c("era5_blh_center", "pblhmod_clim", "pblhmod_default")),
               c("pblhmod_default", "pblhmod_clim", "era5_blh_center"))
  expect_error(style_of("nope"), "No style")
  expect_equal(names(source_palette()), names(source_labels()))
  expect_equal(regime_labels()[["convective"]], "convective (B&G)")
})

test_that("figures take colours only from pblh_sources_style and the theme only from theme_pblh()", {
  f <- c(list.files(here::here("R"), "[.]R$", full.names = TRUE), list.files(here::here("scripts"), "[.]R$", full.names = TRUE))
  f <- setdiff(f, here::here("R", "plot_utils.R"))
  txt <- unlist(lapply(f, function(p) paste(basename(p), readLines(p, warn = FALSE))))
  expect_equal(grep("\"#[0-9A-Fa-f]{6}\"", txt, value = TRUE), character())
  # colour names in plotting calls (the legacy port's `fill` option is not a colour)
  expect_equal(grep("(colour|color|fill) = \"(white|black|gr[ae]y|red|blue|green|orange|yellow|purple)", txt,
                    value = TRUE), character())
  expect_equal(grep("theme_(pub|bw|minimal|classic)[(]", txt, value = TRUE), character())
})
