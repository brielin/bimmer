#' Gets matrix of TCE from observed network.
#'
#' @param R_obs D x D matrix of observed effects.
#' @param normalize A length D vector which is used to convert R_tce from the
#'   per-allele to the per-variance scale. Each entry should be the
#'   std dev of the corresponding phenotype. Set to NULL for no normalization.
#' @return D x D matrix of total causal effects.
get_tce <- function(R_obs, normalize = NULL) {
  diag_R_obs <- diag(R_obs)
  R_tce <- (1 / (1 + diag_R_obs)) * R_obs
  diag(R_tce) <- 1
  if (!is.null(normalize)) {
    R_tce <- R_tce * outer(normalize, 1 / normalize)
  }
  return(R_tce)
}

#' Gets observed network from direct effects.
#'
#' @param R D x D matrix of direct effects.
#' @return D x D matrix of observed effects.
get_observed <- function(R) {
  D <- dim(R)[1]
  return(as.matrix(R %*% solve(diag(D) - R)))
}

#' Gets direct network from observed effects.
#'
#' @param R_obs D x D matrix of direct effects.
#' @return D x D matrix of observed effects.
get_direct <- function(R_obs) {
  D <- dim(R_obs)[1]
  return(R_obs %*% solve(diag(D) + R_obs))
}

#' Uses naive meta analysis method for TCE estimation with multiple SNPs.
#'
#' @param b_exp Vector of SNP effects on exposure variable.
#' @param b_out Vector of SNP effects on outcome variable.
#' @param se_exp Vector of SNP effects on exposure variable.
#' @param se_out Vector of SNP effects on outcome variable.
#' @return A list with two elements.
#'   tce_hat: The total causal efffect estimate.
#'   se_tce: The standard error of the estimate.
naive_ma <- function(b_exp, b_out, se_exp, se_out) {
  tce_hat <- mean(b_out / b_exp)
  se_tce <- stats::sd(b_out / b_exp) / sqrt(length(b_exp))
  return(list("beta.hat" = tce_hat, "beta.se" = se_tce))
}

#' Generates sparse effects and Mendelian randomization conditioning statistic.
#'
#' @description
#' Generates an MxD matrix of SNP effects with p% non-zero in expectation.
#'
#' @param M Integer. Number of SNPs to simulate.
#' @param D Integer. Number of phenotypes to simulate.
#' @param p Float. Proportion of non-zero effects.
#' @param sd Float. Standard deviation normal distribution of effect sizes.
#' @param pleiotropy Bool. TRUE to allow for SNPs to effect multiple phenotypes.
generate_beta <- function(M, D, p = 0.1, sd = 1.0, pleiotropy = FALSE) {
  beta <- Matrix::Matrix(
    stats::rnorm(M * D, sd = sd) * (stats::runif(M * D) < p), M, D,
    sparse = TRUE
  )
  colnames(beta) <- paste0("P", 1:D)
  rownames(beta) <- paste0("rs", 1:M)
  non_zero_index <- !apply(abs(beta) == 0, 1, all)
  beta_non_zero <- beta[non_zero_index, ]
  beta_snp_max <- apply(abs(beta_non_zero), 1, max)
  beta_over_max <- abs(beta_non_zero) / beta_snp_max
  if (!pleiotropy) {
    beta_non_zero[Matrix::which(beta_over_max < 1.0, arr.ind = TRUE)] <- 0
    beta[non_zero_index] <- beta_non_zero
  }
  return(beta)
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
#' @param symmetric Bool. TRUE to force the returned matrix to be
#'   symmetric, false for anti-symmetric, and NULL for no symmetry enforcement.
#' @param sd Float. Standard deviation of network edge weights. NA for binary
#'   matrices with equal probability of
#' @return A DxD sparse matrix.
generate_network <- function(D, p = 0.1, normalize = TRUE, epsilon = 0.1,
                             symmetric = NULL, sd = 1.0) {
  if (is.na(sd)) {
    R <- matrix((2 * stats::rbinom(D * D, 1, 0.5) - 1) *
      (stats::runif(D * D) < p), D, D)
  } else {
    R <- matrix(stats::rnorm(D * D, sd = sd) * (stats::runif(D * D) < p), D, D)
  }
  diag(R) <- 0.0
  if (!is.null(symmetric)) {
    if (symmetric) {
      R[lower.tri(R)] <- 0
      R <- R + t(R)
    }
    else if (!symmetric) {
      # TODO(brielin): Fix this so that it evenly sets upper and lower indices
      #   to 0 and preserves fraction of non-zeros.
      mutual_ut <- which(
        (abs(R)[upper.tri(R)] > 0) & (t(abs(R))[upper.tri(R)] > 0)
      )
      R[upper.tri(R)][mutual_ut] <- 0
    }
  }
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
  return(R)
}

#' Generates genotypes.
#'
#' @param N Integer. Number of individuals to simulate.
#' @param M Integer. Number of SNPs to simulate.
#' @param whiten Bool. TRUE to whiten such that the observed covariance is
#'   exactly 0. This should only be used for testing.
#' @return and NxM matrix of simulated genotypes.
generate_genotypes <- function(N, M, whiten = FALSE) {
  # TODO(brielin): Get real genotypes from Biobank.
  X <- matrix(stats::rnorm(N * M), N, M)
  if (whiten) {
    X <- t(t(X) - colMeans(X))
    X_svd <- svd(X)
    X <- sqrt(N - 1) * X_svd$u %*% t(X_svd$v)
  }
  rownames(X) <- paste0("N", 1:N)
  colnames(X) <- paste0("rs", 1:M)
  return(X)
}

#' Generates (optionally correlated) environmental confounding.
#'
#' Note that this is different from the (always uncorrelated) measurement error
#' since it is mixed by the network with the genetic component.
#'
#' @param N Integer, number of individuals.
#' @param C Integer, dimensionality of confounding.
#' @param D Integer, number of phenotypes.
#' @param sigma_g Float or D x D matrix of floats. If a single value it will be
#'   taken as the variance of the generated effect size. If a CxC matrix it will
#'   be the covariance.
generate_confounding <- function(N, C, D, sigma_g) {
  U <- matrix(stats::rnorm(N * C), N, C)
  if (is.null(dim(sigma_g))) {
    gamma <- matrix(stats::rnorm(C * D, sd = sqrt(sigma_g)), C, D)
  } else if (identical(sigma_g, matrix(1L, D, D))) {
    gamma <- matrix(rep(stats::rnorm(C, sd = sqrt(sigma_g)), D), nrow = C)
  } else {
    gamma <- MASS::mvrnorm(C, mu = rep(0, D), Sigma = sigma_g)
  }
  U %*% gamma
}

#' Generates dataset.
#'
#' Generates a network mendelian randomization dataset with observed phenotypes,
#' genotypes, true network effects and true effect sizes.
#'
#' @param N Integer. Number of individuals to simulate.
#' @param M Integer. Number of SNPs to simulate.
#' @param D Integer. Number of phenotypes to simulate.
#' @param C Integer. Number of confounding components to simulate.
#' @param p_beta Float. Proportion of SNPs with non-zero effect size.
#' @param p_net Float. Proportion of non-zero network edges.
#' @param noise Float between 0 and 1. Proportion of variance in phenotype
#'   attributable to noise. 0 for no noise, 1 for all noise with no
#'   genetic/network/confounding component.
#' @param conf_ratio Float between 0 and 1. Ratio of genetic vs confounding
#'   component prior to mixing by R and adding noiose.
#' @param pleiotropy Bool. TRUE to allow for SNPs to effect multiple phenotypes.
#' @param whiten Bool. TRUE to whiten generated genotypes. See
#'   generate_genotypes for more information.
#' @param symmetric Bool or NULL. Whether to enforce a symmetry constraint on
#'   generated network. See generate_network for more.
#' @param sd_net Float. Standard deviation of network edge weights.
#' @param sd_beta Float. Standard deviation of SNP effect sizes.
#' @param sigma_g Float or CxC matrix of floats. Variance or covariance matrix
#'   of confounding effects.
#' @param fix_R Null or DxD matrix. Use to provide R in order to repeatedly
#'   resample from the same model.
#' @param fix_beta Null or MxD matrix. Use to provide beta in order to
#'   repeatedly resample from the same model.
#' @return A list:
#'   Y: N x D matrix of observed phenotypes.
#'   X: N x M matrix of genotypes.
#'   beta: M x D matrix of genotype effect sizes.
#'   R: D x D network effect matrix.
generate_dataset <- function(N, M, D, C = 0, p_beta = 0.2, p_net = 0.2,
                             noise = 0.0, conf_ratio = 0.0, pleiotropy = FALSE,
                             whiten = FALSE, symmetric = NULL, sd_net = 1.0,
                             sd_beta = 1.0, sigma_g = 1.0, fix_R = NULL,
                             fix_beta = NULL) {
  # Generate the various components.
  if (is.null(fix_beta)) {
    beta <- generate_beta(
      M, D, p_beta, pleiotropy = pleiotropy, sd = sd_beta)
  } else {
    beta <- fix_beta
  }
  if (is.null(fix_R)) {
    R <- generate_network(D, p_net, symmetric = symmetric, sd = sd_net)
  } else {
    R <- fix_R
  }
  X <- generate_genotypes(N, M, whiten = whiten)
  genetic <- X %*% beta

  # Normalize the confouding contribution to have equal variance to the genetic.
  if (C > 0) {
    confounding <- generate_confounding(N, C, D, sigma_g = sigma_g)
    genetic_var <- apply(genetic, 2, stats::var)
    confounding_var <- apply(confounding, 2, stats::var)
    confounding <- t(t(confounding) * sqrt(genetic_var / confounding_var))
  }
  else {
    confounding <- 0
  }

  # Calculate the direct and mixed effect plus noise to produce Y.
  direct <- sqrt(1 - conf_ratio) * genetic + sqrt(conf_ratio) * confounding
  network <- direct %*% solve(diag(D) - R)
  network_var <- apply(network, 2, stats::var)
  epsilon <- t(t(matrix(stats::rnorm(N * D), N, D)) * sqrt(network_var))
  Y <- sqrt(1 - noise) * network + sqrt(noise) * epsilon

  return(list("Y" = as.matrix(Y), "X" = X, "beta" = beta, "R" = R))
}

#' Generate summary statistics.
#'
#' Note that this always returns values on the normalized (per-variance) scale.
#'
#' @param X N x M matrix of genotypes.
#' @param Y N x D matrix of phenotypes.
#' @param N_ind Length D vector with each entry specifying the number of samples
#'   from each phenotype to use.
#' @return A list,
#'   beta_hat: An M x D matrix of calculated effect sizes.
#'   se_hat: An M x D matrix of corresponding estimated standard errors.
generate_sumstats <- function(X, Y, N_ind = NULL, M_null = 0) {
  N <- dim(X)[1]
  M <- dim(X)[2]
  D <- dim(Y)[2]

  if(!is.null(N_ind)){
    select_function <- function(y, n){
      drop = sample.int(N, size=(N-n))
      y[drop] = NA
      return(y)
    }
    Y <- as.matrix(purrr::map2_dfr(data.frame(Y), N_ind, select_function))
  }

  Xs <- scale(X)
  Ys <- scale(Y)
  sumstats_one_pheno <- function(y){
    mask <- !is.na(y)
    N_mask <- sum(mask)
    y_mask <- y[mask]
    if(any(!mask)){
      X_mask <- scale(X[mask, ])
    } else {
      X_mask <- Xs
    }
    beta_hat <- c((t(X_mask) %*% y_mask)/(N_mask-1))
    names(beta_hat) <- colnames(X_mask)
    eps2sum <- colSums((y_mask - t(t(X_mask) * beta_hat))**2)
    se_hat <- sqrt(eps2sum/((N_mask - 1)*(N_mask - 2)))
    return(list(beta_hat, se_hat))
  }

  sumstats <- purrr::transpose(purrr::map(data.frame(Ys), sumstats_one_pheno))
  beta_hat <- as.data.frame(sumstats[[1]])
  se_hat <- as.data.frame(sumstats[[2]])
  if(M_null > 0){
    N_null <- rep(N, D)
    if(!is.null(N_ind)){
      N_null <- N_ind
    }
    beta_null <- t(matrix(stats::rnorm(
      M_null*D, mean = 0, sd = 1/sqrt(N_null)), nrow = D))
    se_null <- t(matrix(rep(1/sqrt(N_null), M_null), nrow = D))
    rownames(beta_null) <- paste0("rs", (M+1):(M+M_null))
    rownames(se_null) <- paste0("rs", (M+1):(M+M_null))
    colnames(beta_null) <- colnames(beta_hat)
    colnames(se_null) <- colnames(se_hat)
    beta_hat <- rbind(beta_hat, beta_null)
    se_hat <- rbind(se_hat, se_null)
  }
  return(list("beta_hat" = beta_hat, "se_hat" = se_hat))
}

#' Selects SNPs based on true effect sizes.
#'
#' @param beta_true M x D matrix of true effect sizes.
select_snps_oracle <- function(beta_true) {
  select <- function(beta) {
    mask <- abs(beta) > 0
    res <- purrr::map(colnames(beta_true), function(x) {
      return(mask)
    })
    names(res) <- colnames(beta_true)
    res <- purrr::prepend(res, list("names" = rownames(beta_true)))
    return(res)
  }
  return(purrr::map(data.frame(as.matrix(beta_true)), select))
}


#' Fits L1-regularized approximate inverse while finding best lambda (cheating).
#'
#' @param R_tce DxD matrix of "total causal effects".
#' @param R DxD matrix of true "causal direct effects".
fit_regularized_cheat <- function(R_tce, R) {
  mae <- Inf
  best_R <- NULL
  for (lambda in c(0.0, 0.0001, 0.001, 0.01, 0.1)) {
    tce_res <- tryCatch(
      fit_regularized(R_tce, lambda),
      error = function(cond) {
        return(NA)
      }
    )
    lam_mae <- mean(abs(R - tce_res$R_hat))
    if (is.na(lam_mae)) {
      lam_mae <- Inf
    }
    if (lam_mae < mae) {
      mae <- lam_mae
      best_R <- tce_res
    }
  }
  return(best_R)
}

#' A simple helper function to count the number of instuments for each pair.
#'
#' @param selected A list of lists, output of `select_snps`.
count_instruments <- function(selected) {
  return(do.call(cbind, purrr::map(selected, function(pheno) {
    return(purrr::map(pheno[-1], sum))
  })))
}

instrument_metrics <- function(selected, beta){
  phenos <- names(selected)
  beta <- as.matrix(beta)
  all_snps <- unique(unlist(purrr::map(selected, function(pheno){pheno$names})))
  null_snps <- all_snps[!(all_snps %in% rownames(beta))]
  beta_null <- matrix(0L, nrow = length(null_snps), ncol = ncol(beta), dimnames = list(null_snps, colnames(beta)))
  beta <- rbind(beta, beta_null)
  purrr::imap(selected, function(select_pheno, pheno1){
    names <- select_pheno$names
    all_p1 <- sum(abs(beta[, pheno1]) > 0)
    purrr::imap(select_pheno[-1], function(use, pheno2){
      dir_eff <- abs(beta[names, , drop=FALSE][use, pheno1]) > 0
      other_eff <- purrr::map_lgl(abs(beta[names, ,drop=FALSE][use, !(colnames(beta) %in% c(pheno1)), drop=FALSE])>0, any)
      correct <- sum(dir_eff & !other_eff)
      reversed <- sum(!dir_eff & other_eff)
      total <- length(dir_eff)
      # Precision, recall, reversed
      return(c(correct/total, correct/all_p1, reversed/total))
    })
  })
}
