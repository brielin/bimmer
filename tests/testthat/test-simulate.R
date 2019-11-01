set.seed(123)

test_that("generate_beta_runs", {
  M <- 5
  D <- 2
  p <- 0.5
  result <- generate_beta(M, D, p)
  expect_equal(dim(result$beta), c(M, D))
  expect_length((result$condition), D)
  expect_is(result$beta, "sparseMatrix")
})

test_that("generate_beta_no_pleiotropy", {
  M <- 5
  D <- 2
  p <- 0.5
  result <- generate_beta(M, D, p, pleiotropy = FALSE)
  expect_equal(unname(Matrix::rowSums(abs(result$beta) > 0) <= 1), rep(TRUE, 5))
})

test_that("generate_network_runs", {
  D <- 10
  p <- 0.5
  R <- generate_network(D, p, normalize = FALSE)
  expect_equal(dim(R), c(D, D))
  expect_is(R, "sparseMatrix")
})

test_that("generate_network_normalizes", {
  D <- 10
  p <- 0.5
  R <- generate_network(D, p)
  eigenvalues <- eigen(R, only.values = TRUE)$values
  expect_lte(max(abs(eigenvalues)), 1.0 + 1e-8)
})

test_that("generate_dataset_runs", {
  N <- 10
  D <- 5
  M <- 7
  p_beta <- 0.5
  p_network <- 0.5
  noise_ratio <- 0.5
  result <- generate_dataset(N, M, D, p_beta, p_network, noise_ratio)
  expect_equal(dim(result$Y), c(N, D))
  expect_equal(dim(result$X), c(N, M))
  expect_equal(dim(result$beta), c(M, D))
  expect_equal(dim(result$R), c(D, D))
})
