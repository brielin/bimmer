set.seed(123)
N <- 100
D <- 3
dataset <- generate_dataset(
  N, D, D,
  p_beta = 1.0, p_net = 0.5, noise = 0.0, pleiotropy = FALSE
)
sumstats <- generate_sumstats(dataset$X, dataset$Y)
R <- dataset$R
R_tce_expected <- matrix(c(
  1, -0.7718971, -0.2976823, -0.7619209, 1, 0.3412960,
  -2.268937, 1.779083, 1
), D, D)

test_that("get_direct_observed_works", {
  expect_equal(matrix(get_direct(get_observed(R))), matrix(R))
})

test_that("get_tce_fit_exact_works", {
  expect_equal(matrix(fit_exact(get_tce(get_observed(R)))), matrix(R))
})

test_that("naive_ma_works", {
  beta_exp <- c(1, 2, 3, 4)
  beta_out <- c(0.5, 2, 1.5, 4)
  se_exp <- c(1, 1, 1, 1)
  se_out <- c(1, 1, 1, 1)
  tce_hat <- 0.75
  se_tce <- 0.144337

  res <- naive_ma(beta_exp, beta_out, se_exp, se_out)
  expect_equal(res$tce_hat, tce_hat)
  expect_equal(res$se_tce, se_tce, tolerance = 1e-6)
})

test_that("fit_tce_mean_works", {
  R_tce_hat <- fit_tce(sumstats, sumstats, "mean")
  expect_is(R_tce_hat, "matrix")
  expect_equal(R_tce_expected, unname(R_tce_hat), tolerance = 1e-6)
})

test_that("select_snps_works", {
  snps <- select_snps(sumstats)
  expect_is(snps, "list")
  expect_equal(names(snps), c("P1", "P2", "P3"))
  expect_equal(names(get("P1", snps)), c("rs1", "rs2"))
  expect_equal(names(get("P2", snps)), c("rs1", "rs2"))
  expect_equal(names(get("P3", snps)), c("rs1", "rs2", "rs3"))
})

test_that("fit_sumstats_exact_works", {
  R_hat <- fit_sumstats(sumstats, sumstats, "mean", "exact")
  expect_equal(unname(fit_exact(R_tce_expected)), unname(R_hat), tolerance = 1e-6)
})

test_that("fit_direct_works", {
  res <- fit_direct(dataset$X, dataset$Y, lambda = 0, niter = 3)
  expect_is(res$R_hat, "matrix")
  expect_is(res$B_hat, "matrix")
  expect_equal(dim(res$R_hat), c(3, 3))
  expect_equal(dim(res$B_hat), c(3, 3))
})
