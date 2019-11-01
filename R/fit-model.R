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
fit_regularized <- function(R_tce, lambda = 0.1) {
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
  diag(D) - t((1 / diag(R_tce_inv)) * t(R_tce_inv))
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
fit_direct <- function(X, Y, lambda = 1., niter = 20, true_R = NULL, true_B = NULL) {
  N <- nrow(X)
  P <- ncol(X)
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
    B_hat <- matrix(0L, nrow = P, ncol = D)
    for (d in 1:D) {
      fit <- glmnet::glmnet(cbind(Y_hat[, -d], X), Y[, d],
        alpha = 1.0, lambda = lambda,
        standardize = FALSE, intercept = FALSE
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
      B_hat[, d] <- fit$beta[D:(D + P - 1)]
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
#' @return D x D matrix of total causal effects.
get_tce <- function(R_obs) {
  diag_R_obs <- Matrix::diag(R_obs)
  R_tce <- (1 / (1 + diag_R_obs)) * R_obs
  Matrix::diag(R_tce) <- 1
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
#' @param beta_exp Vector of SNP effects on exposure variable.
#' @param beta_out Vector of SNP effects on outcome variable.
#' @param se_exp Vector of SNP effects on exposure variable.
#' @param se_out Vector of SNP effects on outcome variable.
#' @return A list with two elements.
#'   tce_hat: The total causal efffect estimate.
#'   se_tce: The standard error of the estimate.
naive_ma <- function(beta_exp, beta_out, se_exp, se_out) {
  tce_hat <- mean(beta_out / beta_exp)
  se_tce <- stats::sd(beta_out / beta_exp) / sqrt(length(beta_exp))
  list("tce_hat" = tce_hat, "se_tce" = se_tce)
}

#' Selects SNP sets for each phenotype.
#'
#' @param sumstats_select List with elements "r_squared" and "p_value", each
#'   M x D matrices.
#' @param p_thresh Float, p-value threshold to use for SNP inclusion.
#' @return A list with names from the columns of r_squared (phenotypes), where
#'   each entry is a vector with names from the rows (SNPs) consististing of the
#'   r_squared value for each SNP that has a significant effect on that
#'   phenotype accoring to p_thresh.
select_snps <- function(sumstats_select, p_thresh = 1e-5) {
  snp_names <- rownames(sumstats_select$r_squared)
  select_one_pheno <- function(p_value, r_squared) {
    sig_index <- p_value < p_thresh
    names(r_squared) <- snp_names
    r_squared[sig_index]
  }
  mapply(select_one_pheno,
    data.frame(sumstats_select$p_value),
    data.frame(sumstats_select$r_squared),
    SIMPLIFY = FALSE
  )
}

#' Calculates matrix of total causal effects using a specified method.
#'
#' @param sumstats_select List representing summary statistics from the first
#'   dataset. Must include "p-value" and "r-squared".
#' @param sumstats_fit List representing summary statistics from the second
#'   dataset. Must include entries "beta_hat" and "se_hat".
#' @param mr_method String, one of c("mean", "raps"). Method to use for TCE
#'   estimate between every pair.
#' @param p_thresh Float. p-value threshold for inclusion.
#' @returns A modified version of sumtats_fit or dataset_fit with only SNPs
#'   that are selected for further analysis.
fit_tce <- function(sumstats_select,
                    sumstats_fit,
                    mr_method = c("mean", "raps"),
                    p_thresh = 1e-5) {
  mr_method_func <- switch(mr_method,
    mean = naive_ma,
    raps = mr.raps::mr.raps
  )

  # Get SNP list and corresponding r-squareds.
  selected <- select_snps(sumstats_select, p_thresh)
  pheno_names <- colnames(sumstats_select$beta_hat)

  run_mr_method_on_exp <- function(exp_name) {
    exp_snps <- get(exp_name, selected)
    run_mr_method <- function(out_name) {
      out_snps <- get(out_name, selected)
      # Use only SNPs that explain more variance in exposure than outcome.
      # TODO(brielin): Consider ways to force that this difference be large, or
      #   do some weighting based on the difference, or have an exclusion set.
      snps_to_use <- exp_snps > out_snps[names(exp_snps)]
      snps_to_use[is.na(snps_to_use)] <- TRUE # Use any SNPs not in out_snps.
      snps_to_use <- names(snps_to_use[snps_to_use])
      beta_exp <- sumstats_fit$beta_hat[snps_to_use, exp_name, drop = FALSE]
      se_exp <- sumstats_fit$se_hat[snps_to_use, exp_name, drop = FALSE]
      beta_out <- sumstats_fit$beta_hat[snps_to_use, out_name, drop = FALSE]
      se_out <- sumstats_fit$se_hat[snps_to_use, out_name, drop = FALSE]
      # TODO(brielin): Do something with the SE of this estimate
      mr_method_func(beta_exp, beta_out, se_exp, se_out)$tce_hat
    }
    row_res <- sapply(pheno_names, run_mr_method)
    row_res[exp_name] <- 1
    return(row_res)
  }
  t(sapply(pheno_names, run_mr_method_on_exp))
}

#' Fits network mendelian randomization model to data.
#'
#' @param sumstats_select List representing summary statistics from the first
#'   dataset. Must include "p-value" and "r-squared".
#' @param sumstats_fit List representing summary statistics from the second
#'   dataset. Must include entries "beta_hat" and "se_hat".
#' @param mr_method String, one of c("mean", "raps"). Method to use for TCE
#'   estimate between every pair.
#' @param fit_method String, one of c("exact", "regularized"). Method to use for
#'   network optimization.
#' @param p_thresh Float. p-value threshold for inclusion.
#' @returns A modified version of sumtats_fit or dataset_fit with only SNPs
#'   that are selected for further analysis.
fit_sumstats <- function(sumstats_select,
                         sumstats_fit,
                         mr_method = c("mean", "raps"),
                         fit_method = c("exact", "regularized"),
                         p_thresh = 1e-5) {
  fit_method_func <- switch(fit_method,
    exact = fit_exact,
    regularized = fit_regularized
  )
  R_tce_hat <- fit_tce(sumstats_select, sumstats_fit, mr_method, p_thresh)
  fit_method_func(R_tce_hat)
}

#' Fits model using individual level genotypes and sumstats for SNP selection.
#'
#' @param X N x M matrix of genotypes.
#' @param Y N x D matrix of phenotypes.
#' @param sumstats_select List with elements "r_squared" and "p_value", each
#'   M x D matrices.
#' @param p_thresh p-value threshold for inclusion of a SNP in the fit.
#' @param lambda Regularization strength.
#' @param niter Number of iterations.
fit_ind_level <- function(X, Y, sumstats_select, p_thresh = 1e-5,
                          lambda = 0.1, niter = 20) {
  selected <- select_snps(sumstats_select, p_thresh)
  snps_to_use <- unique(unlist(lapply(selected, names), use.names = FALSE))
  fit_direct(X[, snps_to_use], Y, lambda = lambda, niter = niter)
}
