load("../testdata/dataset.Rdata")
load("../testdata/dataset_exact.Rdata")
sumstats <- generate_sumstats(dataset$X, dataset$Y)
R <- dataset$R
D <- nrow(R)

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
  lambda_res <- cv_fit_regularized(R_tce)
  expect_equal(lambda_res$lambda, 1e-4)
  expect_equal(length(lambda_res$scores), 7)
})

test_that("welch_test_works", {
  b1 <- 0
  b2 <- 0.1
  s1 <- 0.01
  s2 <- 0.01
  expect_equal(welch_test(b1, s1, b2, s2), -1)
  expect_equal(welch_test(b2, s2, b1, s1), 1)
  expect_equal(welch_test(b1, 10*s1, b2, 10*s2), 0)
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
  selected <- select_snps_oracle(dataset$beta)
  tce_res <- fit_tce(sumstats, selected, "mean", min_instruments = 1)
  R_shrunk <- shrink_R(tce_res$R_tce, tce_res$SE_tce)
  expect_true(all(abs(R_shrunk) <= abs(tce_res$R_tce)))
})

test_that("fit_tce_min_instruments", {
  selected <- select_snps(sumstats)
  R_tce_hat <- fit_tce(sumstats, selected, "mean", min_instruments = 10)$R_tce
  expect_is(R_tce_hat, "matrix")
  expect_equal(sum(is.na(R_tce_hat)), D * (D - 1))
})

test_that("fit_tce_mean_works", {
  selected <- select_snps_oracle(dataset$beta)
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
  selected <- select_snps_oracle(dataset$beta)
  R_tce_hat <- fit_tce(sumstats, selected, "ps", min_instruments = 1)$R_tce
  expect_is(R_tce_hat, "matrix")
  expect_false(any(is.na(R_tce_hat)))
})

test_that("fit_tce_egger_works", {
  selected <- select_snps_oracle(dataset$beta)
  R_tce_hat <- fit_tce(sumstats, selected, "egger", min_instruments = 1)$R_tce
  expect_is(R_tce_hat, "matrix")
  expect_false(any(is.na(R_tce_hat)))
})

test_that("fit_tce_mbe_works", {
  selected <- select_snps_oracle(dataset$beta)
  R_tce_hat <- fit_tce(sumstats, selected, "mbe", min_instruments = 1)$R_tce
  expect_is(R_tce_hat, "matrix")
  expect_false(any(is.na(R_tce_hat)))
})

test_that("resample_cde_works", {
  selected <- select_snps_oracle(dataset$beta)
  tce_res <- fit_tce(sumstats, selected, "ps", min_instruments = 1)
  tce_res$SE_tce[is.na(tce_res$SE_tce)] <- 0.1
  rs_res <- resample_cde(tce_res$R_tce, tce_res$SE_tce, fit_exact, niter = 10)
  expect_is(rs_res$R_cde, "matrix")
  expect_equal(dim(rs_res$R_cde), dim(tce_res$R_tce))
  expect_is(rs_res$SE_cde, "matrix")
  expect_equal(dim(rs_res$SE_cde), dim(tce_res$SE_tce))
})


test_that("delta_cde_runs", {
  selected <- select_snps_oracle(dataset$beta)
  tce_res <- fit_tce(sumstats, selected, "ps", min_instruments = 1)
  fit_res <- fit_exact(tce_res$R_tce)
  SE <- delta_cde(fit_res$R_tce_inv, tce_res$SE_tce)
  expect_equal(dim(SE), c(3, 3))
})

test_that("filter_tce_filters", {
  R_tce <- matrix(c(0.5, 2, 0.5, 0.5), nrow = 2)
  SE_tce <- matrix(c(5, 1, NaN, 0.5), nrow = 2)
  filt_res <- filter_tce(R_tce, SE_tce, max_nan_perc = 1)
  expect_equal(filt_res$R_tce, matrix(c(NA, NA, NA, 0.5), nrow = 2))
  expect_equal(filt_res$SE_tce, matrix(c(NA, NA, NA, 0.5), nrow = 2))
})

test_that("filter_tce_drops", {
  R_tce <- matrix(c(1, 0, 0, NA), nrow = 2)
  SE_tce <- matrix(c(1, 1, 1, 1), nrow = 2)
  filt_res <- filter_tce(R_tce, SE_tce, max_nan_perc = 0.1)
  expect_equal(filt_res$R_tce, 1)
})

test_that("fit_error_exactly_zero", {
  Y_sds <- apply(dataset_exact$Y, 2, sd)
  sumstats_exact <- generate_sumstats(dataset_exact$X, dataset_exact$Y)
  selected <- select_snps_oracle(dataset_exact$beta)
  R_tce_hat <- fit_tce(
    sumstats_exact, selected, "mean", min_instruments = 1)$R_tce
  R_tce <- as.matrix(get_tce(get_observed(dataset_exact$R), normalize = Y_sds))
  expect_equal(R_tce_hat, R_tce)
})
