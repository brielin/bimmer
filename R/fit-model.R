#' Fits inverse sparse regression model.
#'
#' See also inspre::inspre() for more details.
#'
#' @param R_tce D x D matrix of  "total causal effects".
#' @param W DxD Matrix of weights.
#' @param rho Float. Initial learning rate for ADMM.
#' @param lambda Float, sequence of floats of NULL. L1 regularization strength
#'   on inverse of X. If NULL, a logarithmicallly spaced set of values between
#'   the maximimum absolute off diagonal element of X and lambda_min_ratio
#'   times this value will be used.
#' @param lambda_min_ratio Float, ratio of maximum lambda to minimum lambda.
#' @param nlambda Integer. Number of lambda values to try.
#' @param alpha Float between 0 and 1 or NULL. If > 0, the model will be fit
#'   once with gamma = 0 to find L0, then all subsequent fits will use
#'   gamma = alpha * L0 / D. Set to NULL to provide gamma directly.
#' @param gamma Float or sequence of nlambda floats or NULL. Determinant
#'   regularization strength to use (for each lambda value). It is recommended
#'   to set alpha rather than setting this directly.
#' @param its Integer. Maximum number of iterations.
#' @param delta_target Float. Target change in solution.
#' @param verbose 0, 1 or 2. 2 to print convergence progress for each lambda,
#'   1 to print convergence result for each lambda, 0 for no output.
#' @param train_prop Float between 0 and 1. Proportion of data to use for
#'   training in cross-validation.
#' @param cv_folds Integer. Number of cross-validation folds to perform.
#' @param mu rho modification parameter for ADMM. Rho will be
#'   increased/decreased when the dual constrant and primal constraint are off
#'   by a factor of > mu.
#' @param tau rho modification parameter for ADMM. When called for, rho will be
#'   increased/decreased by the factor tau.
#' @param solve_its Integer, number of iterations of bicgstab/lasso to run
#'   for each U and V update.
#' @param ncores Integer, number of cores to use.
#' @param warm_start Logical. Whether to use previous lambda value result as
#'   starting point for next fit.
#' @export
fit_inspre <- function(R_tce, W = NULL, rho = 10.0, lambda = NULL,
                       lambda_min_ratio = 1e-2, nlambda = 20, alpha = 0,
                       gamma = NULL, its = 100, delta_target = 1e-4,
                       verbose = 1, train_prop = 0.8,
                       cv_folds = 0, mu = 10, tau = 2, solve_its = 3,
                       ncores = 1, warm_start = TRUE){
  D <- dim(R_tce)[1]
  inspre_res <- inspre::inspre(
    X = R_tce, W = W, rho = rho, lambda = lambda,
    lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, alpha = alpha,
    gamma = gamma, its = its, delta_target = delta_target, symmetrize = FALSE,
    verbose = verbose, train_prop = train_prop, cv_folds = cv_folds, mu = mu,
    tau = tau, solve_its = solve_its, ncores = ncores, warm_start = warm_start)
  inspre_res$R_hat <- array(0L, dim = dim(inspre_res$V))
  for(i in 1:length(inspre_res$lambda)){
    inspre_res$R_hat[ , , i] <-
      diag(D) - inspre_res$V[ , , i] / diag(inspre_res$V[ , , i])
  }
  dimnames(inspre_res$R_hat) <- list(rownames(R_tce), colnames(R_tce), inspre_res$lambda)
  return(inspre_res)
}


#' Simple implementation of a welch test.
#'
#' This tests the null hypothesis abs(beta1) = abs(beta2) against the two
#' alternatives abs(beta1) > abs(beta2) and abs(beta1) < abs(beta2).
#'
#' @param beta1 Float or vector of floats, mean of the first sample.
#' @param beta2 Float or vector of floats, mean of the second sample.
#' @param se1 Float or vector of floats, SD of estimate of mu1.
#' @param se2 Float or vector of floats, SD of estimate of mu2.
#' @param welch_thresh Float, p_value threshold for significance.
welch_test <- function(beta1, se1, beta2, se2, welch_thresh = 0.05) {
  n1 <- 1 / (se1**2)
  n2 <- 1 / (se2**2)
  t_val <- (abs(beta2) - abs(beta1)) / sqrt(se1^2 + se2^2)
  nu <- (se1^2 + se2^2)^2 / (se1^4 / (n1 - 1) + se2^4 / (n2 - 1))
  res_12 <- stats::pt(t_val, round(nu))
  res_21 <- 1-res_12
  sig_12 <- res_12 < welch_thresh
  sig_21 <- res_21 < welch_thresh
  res <- sig_12 - sig_21
  res[is.na(t_val)] <- 0
  return(list("pass" = res, "t" = t_val, "nu" = nu) )
}


egger_w <- function(b_exp, b_out, se_exp, se_out, weights){
  lm_res <- summary(stats::lm(
    sign(b_exp)*b_out ~ abs(b_exp),
    weights = weights/(mean(weights)*(se_out**2))))
  beta_hat <- lm_res$coefficients[2, 1]
  beta_se <- lm_res$coefficients[2, 2]/min(lm_res$sigma, 1)
  beta_p <- 2 * (1 - stats::pnorm(abs(beta_hat / beta_se)))
  return(
    list("beta.hat" = beta_hat, "beta.se" = beta_se, "beta.p.value" = beta_p))
}


egger <- function(b_exp, b_out, se_exp, se_out, weights){
  lm_res <- summary(stats::lm(
    sign(b_exp)*b_out ~ abs(b_exp),
    weights = 1/se_out**2))
  beta_hat <- lm_res$coefficients[2, 1]
  beta_se <- lm_res$coefficients[2, 2]/min(lm_res$sigma, 1)
  beta_p <- 2 * (1 - stats::pnorm(abs(beta_hat / beta_se)))
  return(
    list("beta.hat" = beta_hat, "beta.se" = beta_se, "beta.p.value" = beta_p))
}


#' Selects SNPs for inclusion in MR by comparing per-variance effect sizes.
#'
#' Notes: sumstats passed to this function must be computed on the per-variance
#' scale.
#'
#' @param sumstats List with elements "beta_hat", "se_hat", both M x D matrices.
#' @param snps_to_use List or NULL. A list named by phenotypes where each list
#'   entry is a list of SNPs that can be used for that phenotype. Usually
#'   the result of clumping to avoid correlated SNPs.
#' @param p_thresh Float, p-value threshold to use for SNP inclusion.
#' @param exclusive Bool. True to only use SNPs significant for one phenotype
#'   but *not* the other.
#' @param weight Bool. True to store welch-test weights for regression.
#' @param filter Double of NULL. If not NULL, filter variants with welch
#'   statistic less than filter.
#' @param verbose Bool. If true, print phenotype label during iteration.
# @export
# select_snps <- function(sumstats, snps_to_use = NULL, p_thresh = 5e-8,
#                         exclusive = FALSE, weight = TRUE, filter = 1.6, verbose = FALSE) {
#   z_scores <- as.matrix(abs(sumstats$beta_hat / sumstats$se_hat))
#   p_vals <- 2 * (1 - stats::pnorm(z_scores))
#   sig_p_vals <- dplyr::as_tibble(p_vals < p_thresh)
#
#   selected_snps <- list()
#   phenos <- colnames(sumstats$beta_hat)
#   D <- length(phenos)
#   snps <- rownames(sumstats$beta_hat)
#
#   for(index in 1:D){
#     pheno1 <- phenos[index]
#     if(verbose){
#       print(pheno1)
#     }
#     mask1 <- rep(TRUE, length(snps))
#     if (!is.null(snps_to_use)) {
#       p1_snps <- get(pheno1, snps_to_use)
#       mask1 <- (snps %in% p1_snps)
#     }
#     sig1 <- dplyr::pull(sig_p_vals, pheno1)
#     candidate1 <- sig1 & mask1
#     candidate1[is.na(candidate1)] <- FALSE
#     selected_snps[[pheno1]]$names <- snps[candidate1]
#     for(pheno2 in phenos[index:D]){
#       mask2 <- rep(TRUE, length(snps))
#       if (!is.null(snps_to_use)) {
#         p2_snps <- get(pheno2, snps_to_use)
#         mask2 <- (snps %in% p2_snps)
#       }
#       sig2 <- dplyr::pull(sig_p_vals, pheno2)
#       candidate2 <- sig2 & mask2
#       candidate2[is.na(candidate2)] <- FALSE
#       selected_snps[[pheno2]]$names <- snps[candidate2]
#
#       keep <- candidate1 | candidate2
#       b1 <- sumstats$beta_hat[keep, pheno1]
#       b2 <- sumstats$beta_hat[keep, pheno2]
#       s1 <- sumstats$se_hat[keep, pheno1]
#       s2 <- sumstats$se_hat[keep, pheno2]
#
#       welch_res <- welch_test(b1, s1, b2, s2)
#       weights_12 <- -welch_res$t[candidate1[keep]]
#       weights_21 <- welch_res$t[candidate2[keep]]
#
#       if(!is.null(filter)){
#         weights_12 <- ifelse(weights_12 > filter, weights_12, 0.0)
#         weights_21 <- ifelse(weights_21 > filter, weights_21, 0.0)
#       }
#       if(isFALSE(weight)){
#         weights_12 <- (weights_12 != 0)
#         weights_21 <- (weights_21 != 0)
#       }
#       if(isTRUE(exclusive)){
#         weights_12 <- weights_12 * as.numeric(!sig2[candidate1])
#         weights_21 <- weights_21 * as.numeric(!sig2[candidate1])
#       }
#
#       selected_snps[[pheno1]][[pheno2]] <- weights_12
#       selected_snps[[pheno2]][[pheno1]] <- weights_21
#     }
#   }
#   return(selected_snps)
# }

#' Calculates matrix of total causal effects using a specified method.
#'
#' @importFrom dplyr %>%
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
#' @param verbose Bpplean. True to print progress.
#' @param ... Additional parameters to pass to mr_method.
#' @export
fit_tce <- function(sumstats, selected_snps, mr_method = "egger_w",
                    min_instruments = 5, verbose = FALSE, ...) {
  mr_method_func <- switch(mr_method,
    ps = function(b_exp, b_out, se_exp, se_out, weight){
      suppressWarnings(mr.raps::mr.raps(
        b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out))
    },
    aps = function(b_exp, b_out, se_exp, se_out, weight) {
      suppressWarnings(mr.raps::mr.raps(
        b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out,
        over.dispersion = TRUE))
    },
    raps = function(b_exp, b_out, se_exp, se_out, weight) {
      suppressWarnings(mr.raps::mr.raps(
        b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out,
        over.dispersion = TRUE, loss.function = "huber"))
    },
    egger_p = function(b_exp, b_out, se_exp, se_out, ...) {
      input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
      egger_res <- MendelianRandomization::mr_egger(input, robust = FALSE, penalized = TRUE)
      return(list("beta.hat" = egger_res$Estimate, "beta.se" = egger_res$StdError.Est, "beta.p.value" = egger_res$Pvalue.Est))
    },
    egger = egger,
    mbe = function(b_exp, b_out, se_exp, se_out, ...){
      input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
      mbe_res <- MendelianRandomization::mr_mbe(input, seed = NA, stderror = "simple")
      return(list("beta.hat" = mbe_res$Estimate, "beta.se" = mbe_res$StdError, "beta.p.value" = mbe_res$Pvalue))
    },
    median = function(b_exp, b_out, se_exp, se_out, ...){
      input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
      median_res <- MendelianRandomization::mr_median(input)
      return(list("beta.hat" = median_res$Estimate, "beta.se" = median_res$StdError, "beta.p.value" = median_res$Pvalue))
    },
    ivw = function(b_exp, b_out, se_exp, se_out, ...){
      input <- MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
      ivw_res <- MendelianRandomization::mr_ivw(input)
      return(list("beta.hat" = ivw_res$Estimate, "beta.se" = ivw_res$StdError, "beta.p.value" = ivw_res$Pvalue))
    },
    mr_presso = function(b_exp, b_out, se_exp, se_out, ...){
      input <- data.frame("b_exp" = b_exp, "b_out" = b_out, "se_exp" = se_exp, "se_out" = se_out)
      mr_presso_res <- MRPRESSO::mr_presso(data = input, BetaOutcome = "b_out", BetaExposure = "b_exp", SdOutcome = "se_out", SdExposure = "se_exp",
                                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000)$`Main MR results`
      result_index = 2
      if(is.na(mr_presso_res$`Causal Estimate`[[result_index]])){
        result_index = 1
      }
      return(list("beta.hat" = mr_presso_res$`Causal Estimate`[[result_index]], "beta.se" = mr_presso_res$Sd[[result_index]], "beta.p.value" = mr_presso_res$`P-value`[[result_index]]))
    },
    mr_mix = function(b_exp, b_out, se_exp, se_out, ...){
      # TODO(brielin): This seems to flip the result?? Double check this.
      mrmix_res <- MRMix::MRMix(b_exp, -b_out, se_exp, se_out)
      return(list("beta.hat" = mrmix_res$theta, "beta.se" = mrmix_res$SE_theta, "beta.p.value" = mrmix_res$pvalue_theta))
    },
    # TODO(brielin): current implementation (probably) won't work on real data
    # because the global SNP matrix is not LD pruned (just per-phenotype).
    # The SE is also asuming the posterior is normal which is probably wrong.
    cause = function(X, variants){
      params <- suppressWarnings(cause::est_cause_params(X, X$snp))
      cause_res <- suppressWarnings(cause::cause(X=X, variants = variants, param_ests = params, force = TRUE))
      summary_cause <- summary(cause_res)
      quants <- summary_cause$quants[[2]]
      beta.hat <- quants[1, "gamma"]
      beta.se <- (quants[3, "gamma"] - quants[2, "gamma"])/(2*1.96)
      beta.p.value <- summary_cause$p
      return(list("beta.hat" = beta.hat, "beta.se" = beta.se, "beta.p.value" = beta.p.value))
    },
    egger_w = egger_w
  )

  # all.equal returns a STRING if they dont have the same length??
  if (!isTRUE(all.equal(names(sumstats$beta_hat), names(selected_snps)))) {
    common_phenotypes <- intersect(
      colnames(sumstats$beta_hat), names(selected_snps))
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
    run_tce_entry <- function(beta_out, se_out, out){
      mask_or_weight = get(out, snps_to_use)
      snp_mask <- (mask_or_weight > 0) & !is.na(mask_or_weight) & !is.na(beta_out) & !is.na(beta_exp)
      n_instruments <- sum(snp_mask)
      if ((n_instruments < min_instruments) | (exp == out)){
        list("R" = NA, "SE" = NA, "N" = n_instruments, "p" = NA)
      } else {
        tryCatch(
          {
            if(mr_method != "cause"){
              mr_res <- mr_method_func(
                b_exp = beta_exp[snp_mask],
                b_out = beta_out[snp_mask],
                se_exp = se_exp[snp_mask],
                se_out = se_out[snp_mask],
                weight = mask_or_weight[snp_mask],
                ...
              )
            } else {
              X <- tibble::as_tibble(tibble::rownames_to_column(sumstats$beta_hat[, c(exp, out)], var = "snp"))
              X <- X %>% dplyr::rename(beta_hat_1 = exp, beta_hat_2 = out)
              X$seb1 <- sumstats$se_hat[,exp]
              X$seb2 <- sumstats$se_hat[,out]
              X$A1 <- "A"
              X$A2 <- "G"
              variants <- dplyr::filter(X, snp %in% snps_to_use$names[snp_mask])$snp
              mr_res <- mr_method_func(X, variants)
            }

            return(list("R" = mr_res$beta.hat, "SE" = mr_res$beta.se,
                        "N" = n_instruments, "p" = mr_res$beta.p.value))
          },
          error = function(cond) {
            message(c("Error when processing ", exp, " ", out))
            message(cond)
            list("R" = NA, "SE" = NA, "N" = n_instruments, "p" = NA)
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
  p_val <- tce_res %>%
    tidyr::pivot_wider(
      names_from = "out", values_from = "p", id_cols = "exp") %>%
    tibble::column_to_rownames("exp")
  diag(R_tce) <- 1.0
  diag(SE_tce) <- 0.0
  diag(N_obs) <- 0.0
  diag(p_val) <- 1.0

  return(list("R_tce" = as.matrix(R_tce), "SE_tce" = as.matrix(SE_tce),
              "N_obs" = as.matrix(N_obs), "p_val" = as.matrix(p_val)))
}

#' Uses delta method to calculate CDE standard error.
#'
#' @param R_tce_inv DxD matrix. Inverse or approximate inverse of the TCE.
#' @param SE_tce DxD matrix. Standard error of the TCE.
#' @param na.rm Logical. Pass through to sum if SE_tce matrix has NA values.
#' @export
delta_cde <- function(R_tce_inv, SE_tce, na.rm = FALSE) {
  delta_one_entry <- function(Ri_row_k, Ri_col_j) {
    return(sum((SE_tce**2) * outer(Ri_row_k**2, Ri_col_j**2), na.rm = na.rm))
  }
  delta_one_row <- function(Ri_col_j) {
    return(apply(R_tce_inv, 1, delta_one_entry, Ri_col_j = Ri_col_j))
  }
  var_cde <- apply(R_tce_inv, 2, delta_one_row)
  diag(var_cde) <- 0
  return(sqrt(t(t(var_cde) / diag(R_tce_inv))))
}

#' Filters estimate of the graph G based on edge inclusion in CV.
#'
#' @param G D x D matrix
#' @param xi D x D matrix
#' @param thresh Entries in G with corresponding xi < thresh are set to 0
#' @export
filter_G_on_xi <- function(G, xi, thresh = 0.45){
  drop = (xi < thresh) & (xi > 0)
  G[drop] = 0
  return(G)
}

#' Helper function to do basic filtering of the TCE matrix.
#'
#' Large values of R, entries with a high SE, and row/columns with many nans
#' can be removed.
#'
#' @param R_tce Matrix or data.frame. Estimates of TCE.
#' @param SE_tce Matrix or data.frame. Standard errors of the entries in R_tce.
#' @param max_R Float. Set all entries where `abs(R_tce) > max_R` to `NA`.
#' @param max_SE Float. Set all entries whwere `SE > max_SE` tp `NA`.
#' @param max_nan_perc Float. Remove columns and rows that are more than
#'   `max_nan_perc` NAs.
#' @export
filter_tce <- function(R_tce, SE_tce, max_R = 1, max_SE = 0.5, max_nan_perc = 0.5) {
  R_tce[is.nan(SE_tce)] <- NA
  SE_tce[is.nan(SE_tce)] <- NA

  R_too_large <- abs(R_tce) > max_R
  R_tce[R_too_large] <- NA
  SE_tce[R_too_large] <- NA
  SE_too_large <- SE_tce > max_SE
  R_tce[SE_too_large] <- NA
  SE_tce[SE_too_large] <- NA

  row_nan_perc <- rowMeans(is.na(R_tce))
  col_nan_perc <- colMeans(is.na(R_tce))
  max_row_nan = max(row_nan_perc)
  max_col_nan = max(col_nan_perc)
  while((max_row_nan > max_nan_perc) | (max_col_nan > max_nan_perc)){
    if(max_row_nan >= max_col_nan){
      which_max_row_nan <- which.max(row_nan_perc)
      R_tce <- R_tce[-which_max_row_nan, -which_max_row_nan]
      SE_tce <- SE_tce[-which_max_row_nan, -which_max_row_nan]
    } else{
      which_max_col_nan <- which.max(col_nan_perc)
      R_tce <- R_tce[-which_max_col_nan, -which_max_col_nan]
      SE_tce <- SE_tce[-which_max_col_nan, -which_max_col_nan]
    }
    row_nan_perc <- rowMeans(is.na(R_tce))
    col_nan_perc <- colMeans(is.na(R_tce))
    max_row_nan = max(row_nan_perc)
    max_col_nan = max(col_nan_perc)
  }
  return(list("R_tce" = R_tce, "SE_tce" = SE_tce))
}


#' Creates an igraph from CDE matrix.
#'
#' @param R_cde Matrix of causal effects.
#' @param min_edge_value Minimum edge strength for pruning.
#' @param max_edge_value Set edges above this number to this.
#' @export
make_igraph <- function(R_cde, min_edge_value = 0.01, max_edge_value = 0.999){
  adj_matrix <- R_cde
  adj_matrix[abs(adj_matrix) < min_edge_value] = 0
  adj_matrix[abs(adj_matrix) > max_edge_value] = max_edge_value
  zeros <- adj_matrix == 0
  adj_matrix <- -log(abs(adj_matrix))
  adj_matrix[zeros] = 0
  return(igraph::graph_from_adjacency_matrix(
    adj_matrix, mode = "directed", weighted = TRUE))
}
