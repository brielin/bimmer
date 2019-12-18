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

test_that("generate_network_symmetry", {
  D <- 10
  p <- 0.5
  R <- as.matrix(generate_network(D, p, symmetric = TRUE))
  expect_true(all(R == t(R)))
  R <- as.matrix(generate_network(D, p, symmetric = FALSE))
  expect_false(any((R > 0) & (t(R) > 0)))
})

test_that("generate_genotypes_runs", {
  N <- 10
  M <- 5
  X <- generate_genotypes(N, M)
  expect_is(X, "matrix")
  expect_equal(dim(X), c(N, M))
})

test_that("generate_genotypes_whitens", {
  N <- 10
  M <- 5
  X <- generate_genotypes(N, M, whiten = TRUE)
  expect_equal(unname(cov(X)), diag(M))
})

test_that("generate_confounding_works", {
  N <- 10
  D <- 2
  C <- 5

  Ug <- generate_confounding(N, C, D, sigma_g = 1.0)
  expect_equal(dim(Ug), c(N, D))

  sigma_g <- matrix(c(1, 0.5, 0.5, 1), nrow = D)
  Ug <- generate_confounding(N, C, D, sigma_g = sigma_g)
  expect_equal(dim(Ug), c(N, D))
})

test_that("generate_dataset_runs", {
  N <- 10
  D <- 5
  M <- 7
  p_beta <- 0.5
  p_net <- 0.5
  noise_ratio <- 0.5
  result <- generate_dataset(N, M, D,
    p_beta = p_beta, p_net = p_net,
    noise = noise_ratio,
    pleiotropy = TRUE
  )
  expect_equal(dim(result$Y), c(N, D))
  expect_equal(dim(result$X), c(N, M))
  expect_equal(dim(result$beta), c(M, D))
  expect_equal(dim(result$R), c(D, D))
})

test_that("generate_dataset_with_confounding", {
  N <- 10
  D <- 5
  C <- 5
  M <- 7
  p_beta <- 0.5
  p_net <- 0.5
  noise_ratio <- 0.5
  result <- generate_dataset(N, M, D, C,
    p_beta = p_beta, p_net = p_net, noise_ratio,
    pleiotropy = TRUE, conf_ratio = 0.5, sigma_g = 1.0
  )
  expect_equal(dim(result$Y), c(N, D))
  expect_equal(dim(result$X), c(N, M))
  expect_equal(dim(result$beta), c(M, D))
  expect_equal(dim(result$R), c(D, D))
})


test_that("generate_sumstats_works", {
  N <- 10
  D <- 2
  M <- 2
  p_beta <- 1.0
  p_net <- 0.5
  dataset <- generate_dataset(N, M, D, p_beta = p_beta, p_net = p_net, pleiotropy = TRUE)
  sumstats <- generate_sumstats(dataset$X, dataset$Y)
  dataset$Y <- scale(dataset$Y)
  dataset$X <- scale(dataset$X)
  expect_equal(
    sumstats$beta_hat[1, 1],
    unname(lm(dataset$Y[, 1] ~ dataset$X[, 1])$coefficients)[2]
  )
  expect_equal(
    sumstats$beta_hat[1, 2],
    unname(lm(dataset$Y[, 2] ~ dataset$X[, 1])$coefficients)[2]
  )
  expect_equal(
    sumstats$beta_hat[2, 1],
    unname(lm(dataset$Y[, 1] ~ dataset$X[, 2])$coefficients)[2]
  )
  expect_equal(
    sumstats$beta_hat[2, 2],
    unname(lm(dataset$Y[, 2] ~ dataset$X[, 2])$coefficients)[2]
  )
})


test_that("select_snps_oracle_works", {
  N <- 10
  D <- 2
  M <- 10
  p_beta <- 1.0
  p_net <- 0.5
  dataset <- generate_dataset(N, M, D, p_beta = p_beta, p_net = p_net, pleiotropy = FALSE)
  sumstats <- generate_sumstats(dataset$X, dataset$Y)
  snps <- select_snps_oracle(dataset$beta)
  expect_is(snps, "list")
  expect_equal(names(snps), c("P1", "P2"))
  expect_is(get("P1", snps), "list")
  expect_is(get("P2", snps), "list")
  expect_is(get("P2", get("P1", snps)), "logical")
  expect_equal(length(get("P2", get("P1", snps))), M)
})
