#' Generates sparse effects and Mendelian randomization conditioning statistic.
#'
#' @description
#' Generates an MxD matrix of SNP effects with p% non-zero in expectation. Also
#' returns a vector of MR conditioning, the maximum ratio of the effect size of
#' any SNP on each gene relative to every other gene.
#'
#' @param M Integer. Number of SNPs to simulate.
#' @param D Integer. Number of phenotypes to simulate.
#' @param p Float. Proportion of non-zero effects.
#' @param pleiotropy Bool. TRUE to allow for SNPs to effect multiple phenotypes.
#' @return A list,
#'   beta: sparse float matrix of dimension MxD.
#'   condition: float vector of length D.
#' @example generate_beta(100, 10, 0.2)
generate_beta <- function(M, D, p = 0.1, pleiotropy = FALSE) {
  # TODO(brielin): Separate beta generation from condition analysis?
  beta <- Matrix::Matrix(
    stats::rnorm(M * D) * (stats::runif(M * D) < p), M, D,
    sparse = TRUE
  )
  colnames(beta) <- paste0("P", 1:D)
  rownames(beta) <- paste0("rs", 1:M)
  non_zero_index <- !apply(abs(beta) == 0, 1, all)
  beta_non_zero <- beta[non_zero_index, ]
  beta_snp_max <- apply(abs(beta_non_zero), 1, max)
  beta_over_max <- abs(beta_non_zero) / beta_snp_max
  max_index <- Matrix::which(beta_over_max == 1.0, arr.ind = TRUE)
  if (!pleiotropy) {
    beta_non_zero[Matrix::which(beta_over_max < 1.0, arr.ind = TRUE)] <- 0
    beta[non_zero_index] <- beta_non_zero
  }
  beta_over_max[max_index] <- 0
  second_max <- apply(abs(beta_over_max), 1, max)
  beta_over_max[max_index] <- 1.0 / (second_max)
  condition <- apply(beta_over_max, 2, max)
  list("beta" = beta, "condition" = condition)
}

#' Generates and optionally normalizes causal graph.
#'
#' Generates a DxD sparse matrix with normally distributed edge weights.
#' Optionally normalize the matrix to have eigenvalues with modulus between
#' -1 and 1.
#'
#' @param D Integer. Number of phenotypes to simulate.
#' @param p Float. Proportion of non-zero edges.
#' @param normalize Bool. TRUE to normalize the generated matrix to have
#'   eigenvalues with modulus between -1 and 1.
#' @param epsilon Float. Number to add to maximum eigenvalue to better-condition
#'   normalization.
#' @return A DxD sparse matrix.
#' @example generate_network(10, 0.2)
generate_network <- function(D, p = 0.1, normalize = TRUE, epsilon = 0.1) {
  R <- matrix(stats::rnorm(D * D) * (stats::runif(D * D) < p), D, D)
  diag(R) <- 0.0
  R <- Matrix::Matrix(R, sparse = TRUE)
  colnames(R) <- paste0("P", 1:D)
  rownames(R) <- paste0("P", 1:D)
  if (normalize == TRUE) {
    values <- eigen(R, symmetric = FALSE, only.values = TRUE)$values
    max_value <- max(abs(values))
    if (max_value > 1.0) {
      R <- R / (max_value + epsilon)
    }
  }
  R
}

#' Generates genotypes.
#'
#' @param N Integer. Number of individuals to simulate.
#' @param M Integer. Number of SNPs to simulate.
#' @return and NxM matrix of simulated genotypes.
#' @example generate_genotype(100, 10)
generate_genotypes <- function(N, M) {
  # TODO(brielin): Get real genotypes from Biobank.
  X <- matrix(stats::rnorm(N * M), N, M)
  # vars <- apply(X, 2, stats::var)
  rownames(X) <- paste0("N", 1:N)
  colnames(X) <- paste0("rs", 1:M)
  # t(t(X)/sqrt(vars))
  return(X)
}

#' Generates dataset.
#'
#' Generates a network mendelian randomization dataset with observed phenotypes,
#' genotypes, true network effects and true effect sizes.
#'
#' @param N Integer. Number of individuals to simulate.
#' @param M Integer. Number of SNPs to simulate.
#' @param D Integer. Number of phenotypes to simulate.
#' @param p_beta Float. Proportion of SNPs with non-zero effect size.
#' @param p_net Float. Proportion of non-zero network edges.
#' @param noise Float. Proprotion of variance in phenotype attributable to
#'   noise. 0 for no noise, 1 for all noise with no genetic/network component.
#' @param pleiotropy Bool. TRUE to allow for SNPs to effect multiple phenotypes.
generate_dataset <- function(N, M, D, p_beta = 0.2, p_net = 0.2, noise = 0.0,
                             pleiotropy = FALSE) {
  X <- generate_genotypes(N, M)
  if ((M == D) && pleiotropy == FALSE) {
    beta_res <- list("beta" = diag(D), "condition" = rep(1e5, D))
  } else {
    beta_res <- generate_beta(M, D, p_beta, pleiotropy)
  }
  R <- generate_network(D, p_net)
  network <- X %*% beta_res$beta %*% solve(diag(D) - R)
  network_var <- apply(network, 2, stats::var)
  epsilon <- t(t(matrix(stats::rnorm(N * D), N, D)) * sqrt(network_var))
  Y <- sqrt(1 - noise) * network + sqrt(noise) * epsilon
  list("Y" = as.matrix(Y), "X" = X, "beta" = beta_res$beta, "R" = R)
}

#' Generate summary statistics.
#'
#' @param X N x M matrix of genotypes.
#' @param Y N x D matrix of phenotypes.
#' @return A list,
#'   beta_hat: An M x D matrix of calculated effect sizes.
#'   se_hat: An M x D matrix of corresponding estimated standard errors.
#'   p_value: An M x D matrix of corresponding p-values.
#'   r_squared: An M x D matrix of variances explained.
generate_sumstats <- function(X, Y) {
  # TODO(brielin): Add test for this function.
  N <- dim(X)[1]
  M <- dim(X)[2]
  D <- dim(Y)[2]
  X_no_mean <- t(t(X) - colMeans(X))
  Y_no_mean <- t(t(Y) - colMeans(Y))
  diag_xtxi <- 1 / colSums(X_no_mean^2)
  beta_hat <- (diag_xtxi * t(X_no_mean)) %*% Y_no_mean
  calc_s <- function(X_col, beta_row) {
    eps2sum <- colSums((Y_no_mean - outer(X_col, beta_row))^2)
    X2sum <- sum(X_col^2)
    sqrt(eps2sum / ((N - 2) * X2sum))
  }
  se_hat <- t(mapply(calc_s, data.frame(X), data.frame(t(beta_hat))))
  Z_hat <- beta_hat / se_hat
  p_value <- 2 * (1 - stats::pnorm(abs(Z_hat)))

  Xvar <- apply(X, 2, stats::var)
  Yvar <- apply(Y, 2, stats::var)
  normalizer <- outer(sqrt(Xvar), sqrt(1 / Yvar))
  r_squared <- (normalizer * beta_hat)^2
  list(
    "beta_hat" = beta_hat,
    "se_hat" = se_hat,
    "p_value" = p_value,
    "r_squared" = r_squared
  )
}
