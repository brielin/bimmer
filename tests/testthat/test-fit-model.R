G <- generate_network(D = 3, pert = T, v_max = 0.3, v = 0.1)
dataset <- generate_dataset(N = 10000000, D = 3, G = G, M = 20, p = 1, h = 0.1)
sumstats <- dataset$sumstats_select
D <- nrow(G)


test_that("fit_exact_inverse", {
  expect_equal(
    matrix(fit_exact(get_tce(get_observed(G)))$R_inv),
    matrix(solve(get_tce(get_observed(G))))
  )
})


test_that("fit_inspre_works", {
  R <- get_tce(get_observed(G))
  fit_res <- fit_inspre(R, lambda = 0, verbose = 0)
  R_hat <- fit_res$R_hat[,,1]
  expect_is(R_hat, "matrix")
  # TODO(brielin): consider returing to an exactness test
  # expect_equal(as.matrix(R), R_hat, tolerance = 0.01)
})


test_that("welch_test_works", {
  b1 <- 0
  b2 <- 0.1
  s1 <- 0.01
  s2 <- 0.01
  expect_equal(welch_test(b1, s1, b2, s2)$pass, -1)
  expect_equal(welch_test(b2, s2, b1, s1)$pass, 1)
  expect_equal(welch_test(b1, 10*s1, b2, 10*s2)$pass, 0)
})


# test_that("select_snps_works", {
#   snps <- select_snps(sumstats)
#   expect_is(snps, "list")
#   expect_equal(names(snps), c("P1", "P2", "P3"))
#   expect_is(get("P1", snps), "list")
#   expect_is(get("P2", snps), "list")
#   expect_is(get("P3", snps), "list")
#   expect_is(get("P2", get("P1", snps)), "numeric")
# })
#
#
# test_that("fit_tce_min_instruments", {
#   selected <- select_snps(sumstats)
#   R_hat <- fit_tce(sumstats, selected, "mean", min_instruments = 100)$R_tce
#   expect_is(R_hat, "matrix")
#   expect_equal(sum(is.na(R_hat)), D * (D - 1))
# })


# test_that("fit_tce_raps_works", {
#   selected <- WWER::select_snps(sumstats)
#   R_hat <- fit_tce(sumstats, selected, "ps", min_instruments = 1)$R_tce
#   expect_is(R_hat, "matrix")
#   expect_false(any(is.na(R_hat)))
# })
#
#
# test_that("fit_tce_egger_works", {
#   selected <- WWER::select_snps(sumstats)
#   R_hat <- fit_tce(sumstats, selected, "egger", min_instruments = 1)$R_tce
#   expect_is(R_hat, "matrix")
#   expect_false(any(is.na(R_hat)))
# })
#
#
# test_that("fit_tce_mbe_works", {
#   selected <- WWER::select_snps(sumstats)
#   R_hat <- fit_tce(sumstats, selected, "mbe", min_instruments = 1)$R_tce
#   expect_is(R_hat, "matrix")
#   expect_false(any(is.na(R_hat)))
# })


test_that("delta_cde_runs", {
  selected <- WWER::select_snps(sumstats)
  tce_res <- fit_tce(sumstats, selected, "ps", min_instruments = 1)
  fit_res <- fit_exact(tce_res$R_tce)
  SE <- delta_cde(fit_res$R_inv, tce_res$SE_tce)
  expect_equal(dim(SE), c(3, 3))
})


test_that("filter_tce_filters", {
  R <- matrix(c(0.5, 2, 0.5, 0.5), nrow = 2)
  SE <- matrix(c(5, 1, NaN, 0.5), nrow = 2)
  filt_res <- filter_tce(R, SE, max_nan_perc = 1)
  expect_equal(filt_res$R_tce, matrix(c(NA, NA, NA, 0.5), nrow = 2))
  expect_equal(filt_res$SE_tce, matrix(c(NA, NA, NA, 0.5), nrow = 2))
})


test_that("filter_tce_drops", {
  R <- matrix(c(rep(0, 8), NA), nrow = 3)
  SE <- matrix(rep(0, 9), nrow = 3)
  filt_res <- filter_tce(R, SE, max_nan_perc = 0.1)
  expect_equal(filt_res$R_tce, matrix(0L, nrow = 2, ncol = 2))
})
