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
