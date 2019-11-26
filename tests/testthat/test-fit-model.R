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
  expect_equal(matrix(fit_exact(get_tce(get_observed(R)))$R_hat), matrix(R))
})

test_that("fit_exact_inverse", {
  expect_equal(
    matrix(fit_exact(get_tce(get_observed(R)))$R_tce_inv),
    matrix(solve(get_tce(get_observed(R))))
  )
})

test_that("fit_regularized_works", {
  R_tce <- get_tce(get_observed(R))
  fit_res <- fit_regularized(R_tce, lambda = 0)
  R_hat <- fit_res$R_hat
  R_tce_inv <- fit_res$R_tce_inv
  expect_equal(as.matrix(R), R_hat, tolerance = 0.01)
  expect_equal(dim(R_tce_inv), c(D, D), tolerance = 0.01)
})

test_that("fit_regularized_cv_works", {
  R_tce <- get_tce(get_observed(R))
  fit_res <- fit_regularized_cv(R_tce)
  R_hat <- fit_res$R_hat
  R_tce_inv <- fit_res$R_tce_inv
  expect_equal(dim(R_hat), c(D, D))
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

test_that("welch_test_works", {
  mu1 <- 0
  mu2 <- 1
  s1 <- 0.01
  s2 <- 0.01
  n1 <- 1000
  n2 <- 1000
  expect_equal(welch_test(mu1, mu2, s1, s2, n1, n2), 1)
  expect_equal(welch_test(mu2, mu1, s1, s2, n1, n2), 0)
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

test_that("shinkage_works", {
  selected <- select_snps(sumstats)
  tce_res <- fit_tce(sumstats, selected, "mean", min_instruments = 1)
  R_shrunk <- shrink_R(tce_res$R_tce, tce_res$SE_tce, tce_res$N_obs)
  expect_true(all(abs(R_shrunk) <= abs(tce_res$R_tce)))
})

test_that("fit_tce_min_instruments", {
  selected <- select_snps(sumstats)
  R_tce_hat <- fit_tce(sumstats, selected, "mean")$R_tce
  expect_is(R_tce_hat, "matrix")
  expect_equal(sum(is.na(R_tce_hat)), D * (D - 1))
})

test_that("fit_tce_mean_works", {
  selected <- select_snps(sumstats)
  tce_res <- fit_tce(sumstats, selected, "mean", min_instruments = 1)
  R_tce_hat <- tce_res$R_tce
  SE_tce <- tce_res$SE_tce
  N_obs <- tce_res$N_obs
  expect_is(R_tce_hat, "matrix")
  expect_false(any(is.na(R_tce_hat)))
  expect_is(SE_tce, "matrix")
  expect_equal(dim(SE_tce), c(D, D))
  expect_is(N_obs, "matrix")
  expect_equal(dim(N_obs), c(D, D))
})

test_that("fit_tce_raps_works", {
  selected <- select_snps(sumstats)
  R_tce_hat <- fit_tce(sumstats, selected, "ps", min_instruments = 1)$R_tce
  expect_is(R_tce_hat, "matrix")
  expect_false(any(is.na(R_tce_hat)))
})

test_that("resample_cde_works", {
  selected <- select_snps(sumstats)
  tce_res <- fit_tce(sumstats, selected, "ps", min_instruments = 1)
  rs_res <- resample_cde(tce_res$R_tce, tce_res$SE_tce, fit_exact, niter = 10)
  expect_is(rs_res$R_cde, "matrix")
  expect_equal(dim(rs_res$R_cde), dim(tce_res$R_tce))
  expect_is(rs_res$SE_cde, "matrix")
  expect_equal(dim(rs_res$SE_cde), dim(tce_res$SE_tce))
})


test_that("delta_cde_runs", {
  selected <- select_snps(sumstats)
  tce_res <- fit_tce(sumstats, selected, "ps", min_instruments = 1)
  fit_res <- fit_exact(tce_res$R_tce)
  SE <- delta_cde(fit_res$R_tce_inv, tce_res$SE_tce)
  expect_equal(dim(SE), c(3, 3))
})

test_that("fit_sumstats_exact_works", {
  ss_res <- fit_sumstats(sumstats, sumstats, "mean", "exact",
    min_instruments = 1
  )
  expect_is(ss_res$R_cde, "matrix")
  expect_false(any(is.na(ss_res$R_cde)))
})

test_that("fit_sumstats_exact_delta_works", {
  ss_res <- fit_sumstats(sumstats, sumstats, "mean", "exact",
    min_instruments = 1, resample = "delta"
  )
  expect_is(ss_res$R_cde, "matrix")
  expect_is(ss_res$SE_cde, "matrix")
  expect_false(any(is.na(ss_res$R_cde)))
})

test_that("fit_error_exactly_zero", {
  M <- 20
  dataset_exact <- generate_dataset(
    N, M, D,
    p_beta = 1.0, p_net = 0.5, noise = 0.0, pleiotropy = FALSE,
    sd_net = 0.2, whiten = TRUE
  )
  Y_sds <- apply(dataset_exact$Y, 2, sd)

  sumstats_exact <- generate_sumstats(dataset_exact$X, dataset_exact$Y)
  selected <- select_snps_oracle(dataset_exact$beta, sumstats_exact$p_value)
  R_tce_hat <- fit_tce(sumstats_exact, selected, "mean", min_instruments = 1)$R_tce
  R_tce <- as.matrix(get_tce(get_observed(dataset_exact$R), normalize = Y_sds))
  expect_equal(R_tce_hat, R_tce)
})
