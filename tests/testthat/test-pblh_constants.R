test_that("constants have the specified values", {
  k <- pblh_constants()
  expect_equal(k$k, 0.4)
  expect_equal(k$g, 9.80665)
  expect_equal(k$omega, 7.2921159e-5)
  expect_equal(c(k$bg_A, k$bg_B, k$bg_C), c(0.2, 2.5, 8))
  expect_equal(c(k$nieuwstadt_a, k$nieuwstadt_b), c(3.8, 2.28))
  expect_equal(k$c_n, 0.3)
  # B1 neutral limit is consistent with B2: b / (2 a) = 2.28 / 7.6 = 0.3 = c_n
  expect_equal(k$nieuwstadt_b / (2 * k$nieuwstadt_a), k$c_n)
  # scale height ~ r_d * 288 K / g
  expect_equal(k$scale_height, k$r_d * 288 / k$g, tolerance = 1e-3)
})
