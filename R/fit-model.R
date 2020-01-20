#' Fits exact model to data.
#'
#' @param R_tce D x D matrix of "total causal effects".
#' @return D x D matrix with zero diagonal of deconvoluted direct effects.
fit_exact <- function(R_tce) {
  D <- dim(R_tce)[1]
  R_tce_inv <- solve(R_tce)
  R_hat <- diag(D) - t((1 / diag(R_tce_inv)) * t(R_tce_inv))
  return(list("R_hat" = R_hat, "R_tce_inv" = R_tce_inv))
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
  return(list("R_hat" = R_hat, "R_tce_inv" = R_tce_inv))
}

#' Fits L1-regularized approximate inverse model.
#'
#' @param R_tce D x D matrix of "total causal effects".
#' @param weights Length D vector of per-row weights. Default is 1 for each
#'   observation.
#' @param k Number of zero-values to drop during CV.
#' @param nfolds Number of folds. NULL for k=1 to use leave-one-out CV.
#' @param lambda Vector of lambda values to try.
cv_fit_regularized <- function(R_tce, weights = NULL, k = 1, nfolds = NULL,
                               lambda = NULL) {
  if(is.null(lambda)){
    lambda = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
  }

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
      fold_err <- colSums(pred**2) * round((D - 1) / k) + (1 - pred_one)**2
      error[d * fold, ] <- fold_err
    }
  }
  scores <- colSums(error)
  best_lambda <- lambda[which.min(scores)]
  return(fit_regularized(R_tce, best_lambda))
}


#' Simple implementation of a welch test.
#'
#' This tests the null hypothesis abs(beta1) = abs(beta2) against the one-sided
#' alternative abs(beta1) > abs(beta2).
#'
#' @param beta1 Float or vector of floats, mean of the first sample.
#' @param beta2 Float or vector of floats, mean of the second sample.
#' @param se1 Float or vector of floats, SD of estimate of mu1.
#' @param se2 Float or vector of floats, SD of estimate of mu2.
#' @param welch_thresh Float, p_value threshold for significance.
#' @return Float
welch_test <- function(beta1, se1, beta2, se2, welch_thresh = 0.05) {
  n1 <- 1 / (se1**2)
  n2 <- 1 / (se2**2)
  t_val <- (abs(beta2) - abs(beta1)) / sqrt(se1^2 + se2^2)
  nu <- (se1^2 + se2^2)^2 / (se1^4 / (n1 - 1) + se2^4 / (n2 - 1))
  res <- stats::pt(t_val, round(nu)) < welch_thresh
  res[is.na(t_val)] <- FALSE
  return(res)
}

#' Wrapper for `welch_test` that sets SNPs significant for both to FALSE.
#'
#' @param beta1 Float or vector of floats, mean of the first sample.
#' @param beta2 Float or vector of floats, mean of the second sample.
#' @param se1 Float or vector of floats, SD of estimate of mu1.
#' @param se2 Float or vector of floats, SD of estimate of mu2.
#' @param sig2 Bool or vector of bools indicating whether this SNP is also
#'   significant for phenotype two.
#' @param welch_thresh Float, p_value threshold for significance.
welch_filter_both_sig <- function(beta1, se1, beta2, se2, sig2,
                                  welch_thresh = 0.05) {
  welch_res <- welch_test(beta1, se1, beta2, se2, welch_thresh)
  welch_res[sig2] <- FALSE
  return(welch_res)
}

#' Selects SNPs for inclusion in MR by comparing per-variance effect sizes.
#'
#' Notes: sumstats passed to this function must be computed on the per-variance
#' scale. This function does twice as much work as necessary.
#'
#' @param sumstats List with elements "beta_hat", "se_hat", "n_mat", and
#'   "p_value", each M x D matrices.
#' @param snps_to_use List or NULL. A list named by phenotypes where each list
#'   entry is a list of SNPs that can be used for that phenotype. Usually
#'   the result of clumping to avoid correlated SNPs.
#' @param p_thresh Float, p-value threshold to use for SNP inclusion.
#' @param welch_thresh FLOAT, p-value threshold for Welch test of equal betas.
#' @param verbose Bool. If true, print phenotype label during iteration.
select_snps <- function(sumstats, snps_to_use = NULL, p_thresh = 1e-4,
                        welch_thresh = 0.05, verbose = FALSE) {
  z_scores <- as.matrix(abs(sumstats$beta_hat / sumstats$se_hat))
  p_vals <- 2 * (1 - stats::pnorm(z_scores))
  sig_p_vals <- data.frame(p_vals < p_thresh)
  run_pheno <- function(pheno) {
    if (verbose) {
      print(pheno)
    }
    snp_mask <- sig_p_vals[, pheno]
    if (!is.null(snps_to_use)) {
      pheno_snps <- get(pheno, snps_to_use)
      snp_mask <- (rownames(sumstats$beta_hat) %in% pheno_snps) & snp_mask
    }
    snp_mask[is.na(snp_mask)] <- FALSE
    beta_sub <- sumstats$beta_hat[snp_mask, ]
    se_sub <- sumstats$se_hat[snp_mask, ]
    sig_pv_sub <- sig_p_vals[snp_mask, ]
    beta <- beta_sub[, pheno]
    se <- se_sub[, pheno]
    snps <- purrr::pmap(list(beta_sub, se_sub, sig_pv_sub),
                        welch_filter_both_sig,
                        beta1 = beta, se1 = se, welch_thresh = welch_thresh)
    snps$names <- rownames(beta_sub)
    return(snps)
  }
  snps_to_use <- purrr::map(colnames(sumstats$beta_hat), run_pheno)
  names(snps_to_use) <- colnames(sumstats$beta_hat)
  return(snps_to_use)
}


#' Shrinks the estimate of the total causal effect matrix.
#'
#' @param R_tce A matrix of floats with diagonals equal to 1.0. The estimate
#'   of R_tce to be shunk.
#' @param SE_tce A matrix of floats, standard errors of the entries in `R_tce`.
#' @param lambda Float less than 1.0 or NULL. Amount of shrinkage. NULL
#'   computes shrinkage automatically.
shrink_R <- function(R_tce, SE_tce, lambda = NULL) {
  if (is.null(lambda)) {
    # Diagonal of SE_tce should be 0.
    var_tce <- SE_tce^2
    numerator <- sum(var_tce, na.rm = TRUE)
    denominator <- sum(R_tce^2, na.rm = TRUE) - dim(R_tce)[1]
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
#' @param sumstats List representing summary statistics from the second
#'   dataset. Must include entries "beta_hat" and "se_hat".
#' @param selected_snps A list of lists with names of each equal to the
#'   phenotype names. The inner lists are boolean vectors of length equal to
#'   the number of SNPs and TRUE indicating to use that SNP.
#' @param mr_method String, one of c("mean", "raps"). Method to use for TCE
#'   estimate between every pair.
#' @param min_instruments Integer. Return NA if there are less than
#'   this many instruments for a pair of phenotypes.
#' @param shrink Boolean. True to shrink estimates of R_tce to 0.
#' @param verbose Bpplean. True to print progress.
#' @param ... Additional parameters to pass to mr_method.
fit_tce <- function(sumstats, selected_snps, mr_method = c("raps", "ps", "aps"),
                    min_instruments = 5, shrink = FALSE, verbose = FALSE, ...) {
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

  # all.equal returns a STRING if they dont have the same length??
  if (!isTRUE(all.equal(names(sumstats$beta_hat), names(selected_snps)))) {
    common_phenotypes <- intersect(
      names(sumstats$beta_hat), names(selected_snps))
    sumstats$beta_hat <- dplyr::select(
      sumstats$beta_hat, dplyr::one_of(common_phenotypes))
    sumstats$se_hat <- dplyr::select(
      sumstats$se_hat, dplyr::one_of(common_phenotypes))
    selected_snps <- purrr::map(selected_snps[common_phenotypes], function(x) {
      x[c(common_phenotypes, "names")]
    })
  }

  run_tce_row <- function(snps_to_use, exp) {
    if (verbose) {
      print(exp)
    }
    beta_sub <- sumstats$beta_hat[snps_to_use$names, ]
    se_sub <- sumstats$se_hat[snps_to_use$names, ]
    beta_exp <- beta_sub[, exp]
    se_exp <- se_sub[, exp]
    run_tce_entry <- function(beta_out, se_out, out) {
      snp_mask <- get(out, snps_to_use) & !is.na(beta_out) & !is.na(beta_exp)
      n_instruments <- sum(snp_mask)
      if (n_instruments < min_instruments) {
        list("R" = NA, "SE" = NA, "N" = n_instruments)
      }
      else {
        tryCatch(
          {
            mr_res <- mr_method_func(
              b_exp = beta_exp[snp_mask],
              b_out = beta_out[snp_mask],
              se_exp = se_exp[snp_mask],
              se_out = se_out[snp_mask],
              ...
            )
            return(list("R" = mr_res$beta.hat, "SE" = mr_res$beta.se,
                        "N" = n_instruments))
          },
          error = function(cond) {
            message(c("Error when processing ", exp, " ", out))
            message(cond)
            list("R" = NA, "SE" = NA, "N" = n_instruments)
          }
        )
      }
    }
    return(purrr::pmap_dfr(
      list(beta_sub, se_sub, colnames(beta_sub)), run_tce_entry, .id = "out"))
  }
  tce_res <- purrr::imap_dfr(selected_snps, run_tce_row, .id = "exp")

  R_tce <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "R", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  SE_tce <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "SE", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  N_obs <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "N", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  diag(R_tce) <- 1.0
  diag(SE_tce) <- 0.0
  diag(N_obs) <- 0.0

  if (shrink) {
    if (is.logical(shrink)) {
      R_tce <- shrink_R(R_tce, SE_tce)
    } else if (is.numeric(shrink)) {
      R_tce <- shrink_R(R_tce, SE_tce, lambda = shrink)
    }
  }
  return(list("R_tce" = as.matrix(R_tce), "SE_tce" = as.matrix(SE_tce),
              "N_obs" = as.matrix(N_obs)))
}

#' Resample TCE estimates using observed SEs while imputing missing values.
#'
#' @param R_tce Matrix.
#' @param SE_tce Matrix, same dimensnions as `R_tce`.
#' @param niter Integer, number of resampling iterations.
#' @param impute_function Function to use to impute values of R_tce.
#' @param rmse_target Float. Target RMSE change per iteration to stop.
#' @param verbose Bool. True to print progress.
resample_tce <- function(R_tce, SE_tce, niter = 100, impute_function = NULL,
                         rmse_target = 1e-4, verbose = FALSE) {
  if (any(is.na(R_tce)) && is.null(impute_function)) {
    stop("There are NA values in the R_tce, matrix.
         You must provide an imputation function.")
  }
  run_sum <- rep(0, length(R_tce))
  run_sum_sq <- rep(0, length(R_tce))
  SE <- rep(0, length(R_tce))
  R <- rep(0, length(R_tce))
  for (i in 1:niter) {
    R_tce_i <- stats::rnorm(length(R_tce),
      mean = as.vector(R_tce),
      sd = as.vector(SE_tce)
    )
    R_tce_i <- matrix(R_tce_i, nrow = nrow(R_tce))
    if (!is.null(impute_function)) {
      R_tce_i <- impute_function(R_tce_i)
    }
    run_sum <- run_sum + as.vector(R_tce_i)
    run_sum_sq <- run_sum_sq + as.vector(R_tce_i)^2

    R <- run_sum / i
    if (i > 1) {
      SE_next <- sqrt((run_sum_sq / i - R^2) * (i / (i - 1)))
      rmse_change <- sqrt(mean((SE - SE_next)^2, na.rm = TRUE))
      if (rmse_change < rmse_target){
        break
      }
      SE <- SE_next
    }
    if (verbose) {
      print(c(i, rmse_change))
    }
  }
  R <- matrix(R, nrow = nrow(R_tce))
  SE <- matrix(SE, nrow = nrow(R_tce))
  return(list("R_tce" = R, "SE_tce" = SE))
}

#' Resample CDE estimate using observed standard errors of TCE.
#'
#' @param R_tce DxD matrix of floats. Estimates of R_tce to resample.
#' @param SE_tce DxD matrix of floats. Standard errors of entries in R_tce.
#' @param fit_method_func Function. Method for inferring R_cde given R_tce.
#'   Must return a list with entry "R_hat".
#' @param impute_function Function. Method to use for imputing missing values in
#'   `R_tce`.
#' @param niter Integer. Maximum number of resampling iterations.
#' @param rmse_target Float. Stop when change in SE of CDE is below this.
#' @param verbose Boolean. True to print convergence progress.
resample_cde <- function(R_tce, SE_tce, fit_method_func, impute_function = NULL,
                         niter = 100, rmse_target = 1e-3, verbose = FALSE) {
  if (any(is.na(R_tce)) && is.null(impute_function)) {
    stop("There are NA values in the R_tce, matrix.
         You must provide an imputation function.")
  }
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
    if (!is.null(impute_function)) {
      R_tce_i <- impute_function(R_tce_i)
    }
    R_cde_i <- fit_method_func(R_tce_i)$R_hat
    run_sum <- run_sum + as.vector(R_cde_i)
    run_sum_sq <- run_sum_sq + as.vector(R_cde_i)^2

    R_cde <- run_sum / i
    if (i > 1) {
      SE_cde_next <- sqrt((run_sum_sq / i - R_cde^2) * (i / (i - 1)))
      rmse_change <- sqrt(mean((SE_cde - SE_cde_next)^2, na.rm = TRUE))
      if (rmse_change < rmse_target){
        break
      }
      SE_cde <- SE_cde_next
    }
    if (verbose) {
      print(c(i, rmse_change))
    }
  }
  R_cde <- matrix(R_cde, nrow = nrow(R_tce))
  SE_cde <- matrix(SE_cde, nrow = nrow(R_tce))
  return(list("R_cde" = R_cde, "SE_cde" = SE_cde))
}

#' Uses delta method to calculate CDE standard error.
#'
#' @param R_tce_inv DxD matrix. Inverse or approximate inverse of the TCE.
#' @param SE_tce DxD matrix. Standard error of the TCE.
delta_cde <- function(R_tce_inv, SE_tce) {
  delta_one_entry <- function(Ri_row_k, Ri_col_j) {
    return(sum((SE_tce**2) * outer(Ri_row_k**2, Ri_col_j**2)))
  }
  delta_one_row <- function(Ri_col_j) {
    return(apply(R_tce_inv, 1, delta_one_entry, Ri_col_j = Ri_col_j))
  }
  var_cde <- apply(R_tce_inv, 2, delta_one_row)
  diag(var_cde) <- 0
  return(sqrt(t(t(var_cde) / diag(R_tce_inv))))
}

#' Helper function to do basic filtering of the TCE matrix.
#'
#' Large values of R, entries with a high SE, and row/columns with many nans
#' can be removed.
#'
#' @param R_tce Matrix or data.frame. Estimates of TCE.
#' @param SE_tce Matrix or data.frame. Standard errors of the entries in R_tce.
#' @param R_tce_true Matrix or data.frame. For simulation and testing, remove
#'   columns/rows from the true R_tce as well as our estimate.
#' @param max_R Float. Set all entries where `abs(R_tce) > max_R` to `NA`.
#' @param max_SE Float. Set all entries where `SE_tce > max_SE` to `NA`.
#' @param max_nan_perc Float. Remove columns and rows that are more than
#'   `max_nan_perc` NAs.
filter_tce <- function(R_tce, SE_tce, R_tce_true = NULL, max_R = 1.5,
                       max_SE = 3, max_nan_perc = 0.5) {
  R_tce[is.nan(SE_tce)] <- NA
  SE_tce[is.nan(SE_tce)] <- NA

  R_too_large <- abs(R_tce) > max_R
  R_tce[R_too_large] <- NA
  SE_tce[R_too_large] <- NA

  SE_too_large <- SE_tce > max_SE
  R_tce[SE_too_large] <- NA
  SE_tce[SE_too_large] <- NA

  drop_rows <- rowMeans(is.na(R_tce)) > max_nan_perc
  R_tce <- R_tce[!drop_rows, !drop_rows, drop = FALSE]
  SE_tce <- SE_tce[!drop_rows, !drop_rows, drop = FALSE]
  drop_cols <- colMeans(is.na(R_tce)) > max_nan_perc

  if (!is.null(R_tce_true)) {
    R_tce_true <- R_tce_true[!drop_rows, !drop_rows]
    R_tce_true <- R_tce_true[!drop_cols, !drop_cols]
    return(list("R_tce" = R_tce[!drop_cols, !drop_cols],
                "SE_tce" = SE_tce[!drop_cols, !drop_cols],
                "R_tce_true" = R_tce_true))
  }
  else {
    return(list("R_tce" = R_tce[!drop_cols, !drop_cols],
                "SE_tce" = SE_tce[!drop_cols, !drop_cols]))
  }
}
