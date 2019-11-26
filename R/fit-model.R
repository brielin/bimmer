#' Fits exact model to data.
#'
#' @param R_tce D x D matrix of "total causal effects".
#' @return D x D matrix with zero diagonal of deconvoluted direct effects.
fit_exact <- function(R_tce) {
  D <- dim(R_tce)[1]
  R_tce_inv <- solve(R_tce)
  R_hat <- diag(D) - t((1 / diag(R_tce_inv)) * t(R_tce_inv))
  list("R_hat" = R_hat, "R_tce_inv" = R_tce_inv)
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
  list("R_hat" = R_hat, "R_tce_inv" = R_tce_inv)
}

#' Fits L1-regularized approximate inverse model.
#'
#' @param R_tce D x D matrix of "total causal effects".
#' @param weights Length D vector of per-row weights. Default is 1 for each observation.
#' @param k Number of zero-values to drop during CV.
#' @param nfolds Number of folds. NULL for k=1 to use leave-one-out CV.
#' @param lambda Vector of lambda values to try.
fit_regularized_cv <- function(R_tce, weights = NULL, k = 1, nfolds = NULL,
                               lambda = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)) {
  D <- dim(R_tce)[1]
  if (k == 1) {
    nfolds <- D - 1
  }
  if (is.null(weights)) {
    weights <- rep(1, D)
  }

  resp <- diag(D)
  error <- matrix(0L, nrow = D * nfolds, ncol = length(lambda))
  for (d in 1:D) {
    # The d'th entry in response is 1.0, which we cannot drop.
    for (fold in 1:nfolds) {
      if (k == 1) {
        index <- if (fold < d) fold else fold + 1
      } else {
        index <- sample(c(1:D)[-d], k)
      }
      fit <- glmnet::glmnet(R_tce[-index, ], resp[-index, d],
        alpha = 1.0, lambda = lambda, weights = weights[-index],
        standardize = FALSE, intercept = FALSE
      )
      pred <- glmnet::predict.glmnet(fit, R_tce[index, , drop = FALSE])
      pred_one <- glmnet::predict.glmnet(fit, R_tce[d, , drop = FALSE])
      # Weight such that the zero-rows contribute D-1 entries to score.
      error[d * fold, ] <- colSums(pred**2) * round((D - 1) / k) + (1 - pred_one)**2
    }
  }
  scores <- colSums(error)
  best_lambda <- lambda[which.min(scores)]
  fit_regularized(R_tce, best_lambda)
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
#' @return Float, the p-value corresponding to the test.
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

#' Shrinks the estimate of the total causal effect matrix.
#'
#' @param R_tce A matrix of floats with diagonals equal to 1.0. The estimate
#'   of R_tce to be shunk.
#' @param SE_tce A matrix of floats with diagonals equal to 0.0. Standard errors
#'   of the entries in R_tce.
#' @param N_obs A matrix of ints with diagonals equal to 0.0. The number of
#'   samples (instruments) used in the estimate of R_tce for each entry.
#' @param lambda Float less than 1.0 or NULL. Amount of shrinkage. NULL
#'   computes shrinkage automatically.
shrink_R <- function(R_tce, SE_tce, N_obs, lambda = NULL) {
  if (is.null(lambda)) {
    # TODO(brielin): Something is wrong with this formula.
    # Diagonal of SE_tce should be 0.
    var_tce <- N_obs * (SE_tce)^2
    numerator <- sum(var_tce, na.rm = TRUE)
    # Diagonal of R_tce should be 1, so technically I should subtract D
    # here but that gives bad results.
    denominator <- sum(R_tce^2, na.rm = TRUE)
    lambda <- numerator / denominator
    lambda_s <- max(0, min(1, lambda))
  } else {
    lambda_s <- lambda
  }
  R_tce <- (1 - lambda_s) * R_tce
  diag(R_tce) <- 1.0
  return(R_tce)
}

#' Calculates matrix of total causal effects using a specified method.
#'
#' @param sumstats_fit List representing summary statistics from the second
#'   dataset. Must include entries "beta_hat" and "se_hat".
#' @param selected_snps A list of lists with names of each equal to the
#'   phenotype names. The inner lists are boolean vectors of length equal to
#'   the number of SNPs and TRUE indicating to use that SNP.
#' @param p_thresh Float. p-value threshold for inclusion.
#' @param mr_method String, one of c("mean", "raps"). Method to use for TCE
#'   estimate between every pair.
#' @param min_instruments Integer. Return NA if there are less than
#'   this many instruments for a pair of phenotypes.
#' @param shrink Boolean. True to shrink estimates of R_tce to 0.
#' @param ... Additional parameters to pass to mr_method.
fit_tce <- function(sumstats_fit,
                    selected_snps,
                    mr_method = c("mean", "ps", "aps", "raps"),
                    p_thresh = 1e-5,
                    min_instruments = 5,
                    shrink = FALSE,
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
      n_instruments <- sum(snps_to_use)
      if (n_instruments < min_instruments) {
        list(NA, NA, n_instruments)
      }
      else {
        mr_res <- mr_method_func(
          b_exp = beta_exp[snps_to_use],
          b_out = beta_out[snps_to_use],
          se_exp = stderr_exp[snps_to_use],
          se_out = stderr_out[snps_to_use],
          ...
        )
        list(mr_res$beta.hat, mr_res$beta.se, n_instruments)
      }
    }
    mapply(
      run_paired_pheno,
      snp_list,
      data.frame(sumstats_fit$beta_hat),
      data.frame(sumstats_fit$se_hat)
    )
  }
  res <- t(mapply(
    run_one_pheno,
    selected_snps,
    data.frame(sumstats_fit$beta_hat),
    data.frame(sumstats_fit$se_hat)
  ))
  # Fuck this garbage programming language.
  R_tce <- matrix(as.numeric(res[, c(TRUE, FALSE, FALSE)]),
    nrow = nrow(res),
    dimnames = list(rownames(res), rownames(res))
  )
  SE_tce <- matrix(as.numeric(res[, c(FALSE, TRUE, FALSE)]),
    nrow = nrow(res),
    dimnames = list(rownames(res), rownames(res))
  )
  N_obs <- matrix(as.numeric(res[, c(FALSE, FALSE, TRUE)]),
    nrow = nrow(res),
    dimnames = list(rownames(res), rownames(res))
  )
  diag(R_tce) <- 1.0
  diag(SE_tce) <- 0.0
  diag(N_obs) <- 0.0
  if (shrink) {
    if (is.logical(shrink)) {
      R_tce <- shrink_R(R_tce, SE_tce, N_obs)
    } else if (is.numeric(shrink)) {
      R_tce <- shrink_R(R_tce, SE_tce, N_obs, lambda = shrink)
    }
  }
  list("R_tce" = R_tce, "SE_tce" = SE_tce, "N_obs" = N_obs)
}

#' Resample CDE estimate using observed standard errors of TCE.
#'
#' @param R_tce DxD matrix of floats. Estimates of R_tce to resample.
#' @param SE_tce DxD matrix of floats. Standard errors of entries in R_tce.
#' @param fit_method_func Function. Method for inferring R_cde given R_tce.
#'   Must return a list with entry "R_hat".
#' @param niter Integer. Number of resampling iterations.
resample_cde <- function(R_tce, SE_tce, fit_method_func, niter = 100) {
  run_sum <- rep(0, length(R_tce))
  run_sum_sq <- rep(0, length(R_tce))
  SE_cde <- rep(0, length(R_tce))
  R_cde <- rep(0, length(R_tce))
  for (i in 1:niter) {
    R_tce_i <- stats::rnorm(length(R_tce),
      mean = as.vector(R_tce),
      sd = as.vector(SE_tce)
    )
    R_tce_i <- matrix(R_tce_i, nrow = nrow(R_tce))
    R_cde_i <- fit_method_func(R_tce_i)$R_hat
    run_sum <- run_sum + as.vector(R_cde_i)
    run_sum_sq <- run_sum_sq + as.vector(R_cde_i)^2
    R_cde <- run_sum / i
    SE_cde_next <- sqrt((run_sum_sq / i - R_cde^2) * (i / (i - 1)))
    eps <- mean(abs(SE_cde - SE_cde_next))
    # print(c(i, eps))
    SE_cde <- SE_cde_next
  }
  R_cde <- matrix(R_cde, nrow = nrow(R_tce))
  SE_cde <- matrix(SE_cde, nrow = nrow(R_tce))
  list("R_cde" = R_cde, "SE_cde" = SE_cde)
}

#' Uses delta method to calculate CDE standard error.
#'
#' @param R_tce_inv DxD matrix. Inverse or approximate inverse of the TCE.
#' @param SE_tce DxD matrix. Standard error of the TCE.
delta_cde <- function(R_tce_inv, SE_tce) {
  delta_one_entry <- function(Ri_row_k, Ri_col_j) {
    sum((SE_tce**2) * outer(Ri_row_k**2, Ri_col_j**2))
  }
  delta_one_row <- function(Ri_col_j) {
    apply(R_tce_inv, 1, delta_one_entry, Ri_col_j = Ri_col_j)
  }
  var_cde <- apply(R_tce_inv, 2, delta_one_row)
  diag(var_cde) <- 0
  sqrt(t(t(var_cde) / diag(R_tce_inv)))
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
#' @param shrink Boolean. True to shrink estimates of R_tce to 0.
#' @param min_instruments Integer. Return NA if there are less than
#'   this many instruments for a pair of phenotypes.
#' @param resample Non-negative integer or "delta". The number of resample iterations.
fit_sumstats <- function(sumstats_select,
                         sumstats_fit,
                         mr_method = c("mean", "ps", "aps", "raps"),
                         fit_method = c("exact", "regularized"),
                         p_thresh = 1e-5,
                         shrink = FALSE,
                         min_instruments = 5,
                         resample = 0) {
  fit_method_func <- switch(fit_method,
    exact = fit_exact,
    regularized = fit_regularized_cv
  )
  selected <- select_snps(sumstats_select)
  tce_res <- fit_tce(
    sumstats_fit, selected, mr_method, p_thresh,
    min_instruments,
    shrink = shrink
  )
  if (resample == "delta") {
    fit_res <- fit_method_func(tce_res$R_tce)
    list(
      "R_cde" = fit_res$R_hat,
      "SE_cde" = delta_cde(fit_res$R_tce_inv, tce_res$SE_tce)
    )
  }
  else if (resample) {
    resample_cde(tce_res$R_tce, tce_res$SE_tce, fit_method_func, resample)
  } else {
    list("R_cde" = fit_method_func(tce_res$R_tce)$R_hat, "SE_cde" = NA)
  }
}
