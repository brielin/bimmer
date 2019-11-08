#' Fits exact model to data.
#'
#' @param R_tce D x D matrix of "total causal effects".
#' @return D x D matrix with zero diagonal of deconvoluted direct effects.
fit_exact <- function(R_tce) {
  D <- dim(R_tce)[1]
  R_tce_inv <- solve(R_tce)
  diag(D) - t((1 / diag(R_tce_inv)) * t(R_tce_inv))
}

#' Fits L1-regularized approximate inverse model.
#'
#' @param R_tce D x D matrix of "total causal effects".
#' @param lambda Float. Regularization strength.
fit_regularized <- function(R_tce, lambda = 0.01) {
  D <- dim(R_tce)[1]
  resp <- diag(D)
  R_tce_inv <- matrix(0L, nrow = D, ncol = D)
  for (d in 1:D) {
    fit <- glmnet::glmnet(R_tce, resp[, d],
      alpha = 1.0, lambda = lambda,
      standardize = FALSE, intercept = FALSE
    )
    R_tce_inv[, d] <- matrix(fit$beta)
  }
  R_hat <- diag(D) - t((1 / diag(R_tce_inv)) * t(R_tce_inv))
  colnames(R_hat) <- colnames(R_tce)
  rownames(R_hat) <- rownames(R_tce)
  return(R_hat)
}

#' Fit simple inspre model
#'
#' @param X N x M matrix of genotypes.
#' @param Y N x D matrix of response.
#' @param lambda regularization strength, passed to glmnet.
#' @param niter max iterations.
#' @return A list:
#'   R_hat: D x D matrix of deconvoluted direct effects.
#'   B_hat: M x D matrix of inferred genotype effect sizes.
#'   progress_tibble: A tible containing information on convergence.
#'
#' TODO(brielin): Figure out what to do with the progress of this and why it's
#'   performance is so poor.
fit_direct <- function(X, Y, lambda = 1., niter = 20, true_R = NULL, true_B = NULL) {
  N <- nrow(X)
  M <- ncol(X)
  D <- ncol(Y)
  Y_hat <- matrix(0L, nrow = N, ncol = D)
  R_hat <- matrix(0L, nrow = D, ncol = D)
  for (d in 1:D) {
    fit <- glmnet::glmnet(X, Y[, d],
      alpha = 1.0, lambda = lambda,
      standardize = FALSE, intercept = FALSE
    )
    Y_hat[, d] <- glmnet::predict.glmnet(fit, X)
  }

  `%do%` <- foreach::`%do%`
  progress_tibble <- foreach::foreach(i = 1:niter, .combine = dplyr::bind_rows) %do% {
    R_hat_next <- matrix(0L, nrow = D, ncol = D)
    Y_hat_next <- matrix(0L, nrow = N, ncol = D)
    B_hat <- matrix(0L, nrow = M, ncol = D)
    for (d in 1:D) {
      fit <- glmnet::glmnet(cbind(Y_hat[, -d], X), Y[, d],
        alpha = 1.0,
        lambda = lambda, standardize = FALSE,
        intercept = FALSE
      )
      Y_hat_next[, d] <- glmnet::predict.glmnet(fit, cbind(Y_hat[, -d], X))
      if (d == 1) { # I fucking hate R.
        Rdnext <- c(0, fit$beta[1:(D - 1)])
      } else if (d == D) {
        Rdnext <- c(fit$beta[1:(D - 1)], 0)
      }
      else {
        Rdnext <- c(fit$beta[1:(d - 1)], 0, fit$beta[d:(D - 1)])
      }
      R_hat_next[, d] <- Rdnext
      B_hat[, d] <- fit$beta[D:(D + M - 1)]
    }
    delta_Y <- sum((Y_hat - Y_hat_next)^2)
    delta_R <- sum((R_hat - R_hat_next)^2)
    Y_hat <- Y_hat_next
    R_hat <- R_hat_next

    tibble::tibble(
      iteration = i,
      mae_recover_R = ifelse(is.null(true_R), NA, mean(abs(true_R - R_hat_next))),
      mae_recover_B = ifelse(is.null(true_B), NA, mean(as.matrix(abs(true_B - B_hat)))),
      delta_Y = delta_Y,
      delta_R = delta_R
    )
  }

  list(R_hat = R_hat, B_hat = B_hat, progress_tibble = progress_tibble)
}


#' Gets matrix of TCE from observed network.
#'
#' @param R_obs D x D matrix of observed effects.
#' @param normalize A length D vector which is used to convert R_tce from the
#'   per-allele to the per-variance scale. Each entry should be the
#'   std dev of the corresponding phenotype. Set to NULL for no normalization.
#' @return D x D matrix of total causal effects.
get_tce <- function(R_obs, normalize = NULL) {
  diag_R_obs <- Matrix::diag(R_obs)
  R_tce <- (1 / (1 + diag_R_obs)) * R_obs
  Matrix::diag(R_tce) <- 1
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
  R %*% solve(diag(D) - R)
}

#' Gets direct network from observed effects.
#'
#' @param R_obs D x D matrix of direct effects.
#' @return D x D matrix of observed effects.
get_direct <- function(R_obs) {
  D <- dim(R_obs)[1]
  R_obs %*% solve(diag(D) + R_obs)
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
  list("beta.hat" = tce_hat, "beta.se" = se_tce)
}

#' Simple implementation of a welch test.
#'
#' This tests the null hypothesis mu1 = mu2 against the one-sided alternative
#' mu1 > mu2.
#'
#' @param mu1 Float, mean of the first sample.
#' @param mu2 Float, mean of the second sample.
#' @param s1 Float, SD of estimate of mu1.
#' @param s2 Float, SD of estimate of mu2.
#' @param n1 Integer, number of samples in dataset 1.
#' @param n2 Integer, number of samples in dataset 2.
#' @returns Float, the p-value corresponding to the test.
welch_test <- function(mu1, mu2, s1, s2, n1, n2) {
  t_val <- (mu2 - mu1) / sqrt(s1^2 + s2^2)
  nu <- (s1^2 + s2^2)^2 / (s1^4 / ((n1 - 1)) + s2^4 / ((n2 - 1)))
  stats::pt(t_val, round(nu))
}

#' Selects SNPs for inclusion in MR by comparing per-variance effect sizes.
#'
#' Notes: sumstats passed to this function must be computed on the per-variance
#' scale. This function does twice as much work as necessary.
#'
#' @param sumstats_select List with elements "beta_hat", "se_hat", "n_mat", and
#'   "p_value", each M x D matrices.
#' @param p_thresh Float, p-value threshold to use for SNP inclusion.
#' @param welch_thresh FLOAT, p-value threshold for Welch test of equal betas.
#' @return A list of lists. The outer list is indexed by the phenotype names.
#'   Each inner list is also indexed by the phenotype names, and the value of
#'   result$P1$P2 is a boolean vector of length M corresponding to the SNPs to
#'   use in the estimation of TCE of P1 on P2.
select_snps <- function(sumstats_select, p_thresh = 1e-5, welch_thresh = 0.01) {
  select_function <- function(beta, se, p_val, n_val) {
    run_test <- function(beta_other, se_other, n_other) {
      index <- p_val < p_thresh
      welch_res <- mapply(
        welch_test,
        abs(beta[index]),
        abs(beta_other[index]),
        se[index],
        se_other[index],
        n_val[index],
        n_other[index]
      )
      index[index == TRUE][welch_res > welch_thresh] <- FALSE
      return(index)
    }
    mapply(run_test,
      data.frame(sumstats_select$beta_hat),
      data.frame(sumstats_select$se_hat),
      data.frame(sumstats_select$n_mat),
      SIMPLIFY = FALSE
    )
  }
  result <- mapply(select_function,
    data.frame(sumstats_select$beta_hat),
    data.frame(sumstats_select$se_hat),
    data.frame(sumstats_select$p_value),
    data.frame(sumstats_select$n_mat),
    SIMPLIFY = FALSE
  )
  result
}

#' Calculates matrix of total causal effects using a specified method.
#'
#' @param sumstats_fit List representing summary statistics from the second
#'   dataset. Must include entries "beta_hat" and "se_hat".
#' @param selected A list of lists with names of each equal to the phenotype
#'   names.
#' @param p_thresh Float. p-value threshold for inclusion.
#' @param mr_method String, one of c("mean", "raps"). Method to use for TCE
#'   estimate between every pair.
#' @param ... Additional parameters to pass to mr_method.
#' @returns
fit_tce <- function(sumstats_fit,
                    selected_snps,
                    mr_method = c("mean", "ps", "aps", "raps"),
                    p_thresh = 1e-5,
                    min_instruments = 5,
                    ...) {
  mr_method_func <- switch(mr_method,
    mean = naive_ma,
    ps = mr.raps::mr.raps,
    aps = function(...) {
      mr.raps::mr.raps(over.dispersion = TRUE, ...)
    },
    raps = function(...) {
      mr.raps::mr.raps(over.dispersion = TRUE, loss.function = "huber", ...)
    }
  )

  run_one_pheno <- function(snp_list, beta_exp, stderr_exp) {
    run_paired_pheno <- function(snps_to_use, beta_out, stderr_out) {
      # TODO(brielin): Do something with the SE of this estimate
      if (sum(snps_to_use) < min_instruments) {
        return(NA)
      }
      else {
        mr_method_func(
          b_exp = beta_exp[snps_to_use],
          b_out = beta_out[snps_to_use],
          se_exp = stderr_exp[snps_to_use],
          se_out = stderr_out[snps_to_use],
          ...
        )$beta.hat
      }
    }
    mapply(
      run_paired_pheno,
      snp_list,
      data.frame(sumstats_fit$beta_hat),
      data.frame(sumstats_fit$se_hat)
    )
  }
  R_tce <- t(mapply(
    run_one_pheno,
    selected_snps,
    data.frame(sumstats_fit$beta_hat),
    data.frame(sumstats_fit$se_hat)
  ))
  diag(R_tce) <- 1.0
  return(R_tce)
}

#' Fits network mendelian randomization model to data.
#'
#' @param sumstats_select List representing summary statistics from the first
#'   dataset. Must include "p-value" and "r-squared".
#' @param sumstats_fit List representing summary statistics from the second
#'   dataset. Must include entries "beta_hat" and "se_hat".
#' @param mr_method String, one of c("mean", "raps"). Method to use for TCE
#'   estimate between every pair. Note that this will be called with default
#'   arguments. For more detailed control, call fit_tce() directly.
#' @param fit_method String, one of c("exact", "regularized"). Method to use for
#'   network optimization. Note that this will  be called with the default
#'   arguments. For more detailed control, call the fit method directly.
#' @param p_thresh Float. p-value threshold for inclusion.
#' @returns A modified version of sumtats_fit or dataset_fit with only SNPs
#'   that are selected for further analysis.
fit_sumstats <- function(sumstats_select,
                         sumstats_fit,
                         mr_method = c("mean", "ps", "aps", "raps"),
                         fit_method = c("exact", "regularized"),
                         p_thresh = 1e-5,
                         min_instruments = 5) {
  fit_method_func <- switch(fit_method,
    exact = fit_exact,
    regularized = fit_regularized
  )
  selected <- select_snps(sumstats_select)
  R_tce_hat <- fit_tce(
    sumstats_fit, selected, mr_method, p_thresh,
    min_instruments
  )
  fit_method_func(R_tce_hat)
}
