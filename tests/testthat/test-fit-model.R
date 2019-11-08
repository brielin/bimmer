set.seed(123)
N <- 1000
D <- 3
dataset <- generate_dataset(
  N, D, D,
  p_beta = 1.0, p_net = 0.5, noise = 0.0, pleiotropy = FALSE, sd_net = 0.5
)
sumstats <- generate_sumstats(dataset$X, dataset$Y)
R <- dataset$R

test_that("get_direct_observed_works", {
  expect_equal(matrix(get_direct(get_observed(R))), matrix(R))
})

test_that("get_tce_fit_exact_works", {
  expect_equal(matrix(fit_exact(get_tce(get_observed(R)))), matrix(R))
})

test_that("fit_regularized_works", {
  R_tce <- get_tce(get_observed(R))
  R_hat <- fit_regularized(R_tce, lambda = 0)
  expect_equal(as.matrix(R), R_hat, tolerance = 0.01)
})

test_that("fit_direct_works", {
  res <- fit_direct(dataset$X, dataset$Y, lambda = 0, niter = 3)
  expect_is(res$R_hat, "matrix")
  expect_is(res$B_hat, "matrix")
  expect_equal(dim(res$R_hat), c(3, 3))
  expect_equal(dim(res$B_hat), c(3, 3))
})

test_that("naive_ma_works", {
  beta_exp <- c(1, 2, 3, 4)
  beta_out <- c(0.5, 2, 1.5, 4)
  se_exp <- c(1, 1, 1, 1)
  se_out <- c(1, 1, 1, 1)
  tce_hat <- 0.75
  se_tce <- 0.144337

  res <- naive_ma(beta_exp, beta_out, se_exp, se_out)
  expect_equal(res$beta.hat, tce_hat)
  expect_equal(res$beta.se, se_tce, tolerance = 1e-6)
})

test_that("select_snps_works", {
  snps <- select_snps(sumstats)
  expect_is(snps, "list")
  expect_equal(names(snps), c("P1", "P2", "P3"))
  expect_is(get("P1", snps), "list")
  expect_is(get("P2", snps), "list")
  expect_is(get("P3", snps), "list")
  expect_is(get("P2", get("P1", snps)), "logical")
  expect_equal(length(get("P2", get("P1", snps))), D)
})

test_that("fit_tce_min_instruments", {
  selected <- select_snps(sumstats)
  R_tce_hat <- fit_tce(sumstats, selected, "mean")
  expect_is(R_tce_hat, "matrix")
  expect_equal(sum(is.na(R_tce_hat)), D * (D - 1))
})

test_that("fit_tce_mean_works", {
  selected <- select_snps(sumstats)
  R_tce_hat <- fit_tce(sumstats, selected, "mean", min_instruments = 1)
  expect_is(R_tce_hat, "matrix")
  expect_false(any(is.na(R_tce_hat)))
})

test_that("fit_tce_raps_works", {
  selected <- select_snps(sumstats)
  R_tce_hat <- fit_tce(sumstats, selected, "ps", min_instruments = 1)
  expect_is(R_tce_hat, "matrix")
  expect_false(any(is.na(R_tce_hat)))
})

test_that("fit_sumstats_exact_works", {
  R_hat <- fit_sumstats(sumstats, sumstats, "mean", "exact",
    min_instruments = 1
  )
  expect_is(R_hat, "matrix")
  expect_false(any(is.na(R_hat)))
})

test_that("fit_error_exactly_zero", {
  M <- 10
  dataset_exact <- generate_dataset(
    N, M, D,
    p_beta = 1.0, p_net = 0.5, noise = 0.0, pleiotropy = FALSE,
    sd_net = 0.2, whiten = TRUE
  )
  Y_sds <- apply(dataset_exact$Y, 2, sd)

  sumstats_exact <- generate_sumstats(dataset_exact$X, dataset_exact$Y)
  selected <- select_snps_oracle(dataset_exact$beta, sumstats_exact$p_value)
  R_tce_hat <- fit_tce(sumstats_exact, selected, "mean", min_instruments = 1)
  R_tce <- as.matrix(get_tce(get_observed(dataset_exact$R), normalize = Y_sds))
  expect_equal(R_tce_hat, R_tce)
})
