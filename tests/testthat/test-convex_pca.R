test_that("make_weights_works", {
  SE <- matrix(c(0, 1, 2, NA), nrow = 2)
  weights <- make_weights(SE)
  expect_equal(weights, matrix(c(1, 1, 0.25, 0), nrow = 2))
})

test_that("make_weights_max", {
  SE <- matrix(c(1, 0.1, 0.01, 0.001), nrow = 2)
  weights <- make_weights(SE, max_weight = 200)
  expect_equal(weights, matrix(c(1, 100, 200, 200), nrow = 2))
})

test_that("convex_pca_worker_works", {
  Y <- matrix(rnorm(25), nrow = 5)
  Y[1, 2] <- NA
  Y[3, 5] <- NA

  X <- convex_pca(Y, r = 1)
  expect_equal(sum(is.na(X)), 0)
  expect_equal(sum(svd(X, nu = 0, nv = 0)$d), 1, tol = 1e-4)
})

test_that("find_best_r_runs", {
  Y <- matrix(rnorm(25), nrow = 5)
  r <- find_best_r(Y, rs = c(1, 2), its = 3)
  expect_is(r, "numeric")
  expect_equal(length(r), 1)
})

test_that("convex_pca_r", {
  Y <- matrix(rnorm(25), nrow = 5)
  X <- convex_pca(Y, its = 3)
  expect_is(X, "matrix")
  expect_equal(dim(X), c(5, 5))

  X <- convex_pca(Y, r = 1, its = 3)
  expect_is(X, "matrix")
  expect_equal(dim(X), c(5, 5))

  X <- convex_pca(Y, r = c(1, 2), its = 3)
  expect_is(X, "matrix")
  expect_equal(dim(X), c(5, 5))
})

test_that("convex_pca_pcas", {
  Y <- matrix(rnorm(25), nrow = 5)
  res <- convex_pca(Y, n_components = 2, its = 10)
  expect_is(res, "list")
  expect_length(res$d, 2)
  expect_equal(dim(res$u), c(5, 2))
  expect_equal(dim(res$v), c(5, 2))
})
