# Source all project functions for the tests (R/ holds functions only).
for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

fixture_path <- function(...) testthat::test_path("fixtures", ...)
