test_that("get_direct_observed_works", {
  R <- generate_network(3, g = 1)
  dataset <- generate_dataset(N = 100, D = 3, R = R, M_total = 20, M_s = 5,
                              M_p = 5, prop_shared = NULL, rho = 0, noise = 0.5,
                              p_net = 1, sd_net = 0.1)
  expect_equal(matrix(get_direct(get_observed(dataset$R_cde))), matrix(dataset$R_cde))
})


# TODO(brielin): Why isn't this exact??
# test_that("get_tce_fit_exact_works", {
#   dataset <- generate_dataset(N = 100, D = 3, M_total = 20, M_s = 5, M_p = 5,
#      prop_shared = NULL, rho = 0, noise = 0.5, p_net = 1, sd_net = 0.1)
#   expect_equal(matrix(fit_exact(get_tce(get_observed(dataset$R_cde)))$R_hat),
#      matrix(dataset$R_cde), tolerance = 0.01)
# })


test_that("generate_beta_pair_runs", {
  beta <- generate_beta_pair(5, 10, 0)
  expect_equal(dim(beta), c(25, 2))
})


test_that("generate_beta_runs", {
  result <- generate_beta(5, 10, 3, 0)
  expect_equal(dim(result), c(35, 3))
  result <- generate_beta(5, 10, 5, 0)
  expect_equal(dim(result), c(60, 5))
})


test_that("generate_network_runs", {
  D <- 100
  for(graph in c("random", "hub", "scale-free")){
    for(orient in c("random", "towards", "away")){
      R <- generate_network(D, graph = graph, orient = orient)
      expect_equal(dim(R), c(D, D))
    }
  }
})

test_that("generate_network_normalizes", {
  D <- 20
  R <- generate_network(D)
  eigenvalues <- eigen(R, only.values = TRUE)$values
  expect_lte(max(abs(eigenvalues)), 1.0 + 1e-8)
})

test_that("generate_dataset_runs", {
  R <- generate_network(3, g = 1)
  dataset <- generate_dataset(
    N = 100, D = 3, R = R, M_total = 30, M_s = 5, M_p = 5, prop_shared = NULL,
    rho = 0, noise = 0.5, p_net = 1, sd_net = 0.1)

  expect_equal(dim(dataset$sumstats_select$beta_hat), c(30,3))
  expect_equal(dim(dataset$sumstats_select$se_hat), c(30,3))
  expect_equal(dim(dataset$sumstats_fit$beta_hat), c(30,3))
  expect_equal(dim(dataset$sumstats_fit$se_hat), c(30,3))
  expect_equal(dim(dataset$R_tce), c(3, 3))
  expect_equal(dim(dataset$R_cde), c(3, 3))
  expect_equal(dim(dataset$beta), c(20, 3))
})


test_that("select_snps_oracle_works", {
  R <- generate_network(3, g = 1)
  dataset <- generate_dataset(
    N = 100, D = 3, R = R,  M_total = 30, M_s = 5, M_p = 5, prop_shared = NULL,
    rho = 0, noise = 0.5, p_net = 1, sd_net = 0.1)

  snps <- select_snps_oracle(dataset$beta)
  expect_is(snps, "list")
  expect_equal(names(snps), c("P1", "P2", "P3"))
  expect_is(get("P1", snps), "list")
  expect_is(get("P2", snps), "list")
  expect_is(get("P2", get("P1", snps)), "logical")
  expect_equal(length(get("P2", get("P1", snps))), 5)
})
