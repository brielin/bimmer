#' Fits exact model to data.
#'
#' @param R_tce D x D matrix of "total causal effects".
#' @return D x D matrix with zero diagonal of deconvoluted direct effects.
fit_exact <- function(R_tce) {
  D <- dim(R_tce)[1]
  R_tce[is.na(R_tce)] <- 0
  R_tce_inv <- solve(R_tce)
  R_hat <- diag(D) - R_tce_inv / diag(R_tce_inv)
  return(list("R_hat" = R_hat, "R_tce_inv" = R_tce_inv))
}


#' Fit's MVMR model to data.
#'
#' @importFrom foreach %do%
#'
#' @param sumstats_fit List with elements beta_hat and se_hat.
#' @param snps_to_use List of per-phenotype pair SNP to use.
fit_mvmr <- function(sumstats_fit, snps_to_use){
  snps_to_use <- purrr::transpose(snps_to_use)
  snp_names <- snps_to_use$names
  snps_to_use$names <- NULL
  pheno_names <- names(snps_to_use)

  out = NULL
  R_hat <- as.matrix(foreach::foreach(out = pheno_names, .combine = dplyr::bind_rows) %do% {
    exp_snps <- dplyr::setdiff(unique(unlist(snp_names[names(snp_names) != out])),  snp_names[[out]])
    exp_phenos <- pheno_names[pheno_names != out]
    betaX <- as.matrix(sumstats_fit$beta_hat[exp_snps, exp_phenos])
    betaXse <- as.matrix(sumstats_fit$se_hat[exp_snps, exp_phenos])
    betaY <- sumstats_fit$beta_hat[exp_snps, out]
    names(betaY) <- rownames(betaX)
    betaYse <- sumstats_fit$se_hat[exp_snps, out]
    names(betaY) <- rownames(betaXse)


    mrinput <- MendelianRandomization::mr_mvinput(bx = betaX, bxse = betaXse, by = betaY, byse = betaYse)
    mr_res <- as.vector(MendelianRandomization::mr_mvegger(mrinput)$Estimate)
    names(mr_res) <- exp_phenos
    self <- c(0.0)
    names(self) = out
    return(c(self, mr_res))
  })

  rownames(R_hat) <- colnames(R_hat)
  return(list("R_hat" = R_hat))
}


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
  return(solve(diag(D) - R, R))
}


#' Gets direct network from observed effects.
#'
#' @param R_obs D x D matrix of direct effects.
#' @return D x D matrix of observed effects.
get_direct <- function(R_obs) {
  D <- dim(R_obs)[1]
  return(solve(diag(D) + R_obs, R_obs))
}


#' Generates pairs of optionally correleated SNP effects.
#'
#' @param M_s Integer. Number of shared SNPs.
#' @param M_p Integer or sequence of two integers. Number of private SNPs to
#'  simulate per phenotype pair.
#' @param rho Float [-1, 1]. Correlation of shared SNPs.
generate_beta_pair <- function(M_s, M_p, rho){
  beta_s <- mvtnorm::rmvnorm(M_s, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  if(length(M_p) > 1){
    beta_p1 <- stats::rnorm(M_p[1])
    beta_p2 <- stats::rnorm(M_p[2])
    beta <- rbind(cbind(beta_p1, rep(0, M_p[1])), cbind(rep(0, M_p[2]), beta_p2), beta_s)
  } else{
    beta_p1 <- stats::rnorm(M_p)
    beta_p2 <- stats::rnorm(M_p)
    beta <- rbind(cbind(beta_p1, rep(0, M_p)), cbind(rep(0, M_p), beta_p2), beta_s)
  }
  return(beta)
}


#' Generates many pairs of optionally correlated SNP effects
#'
#'
#' @param M_s Integer. Number of shared direct effects to simulate per phenotype.
#' @param M_p Integer or sequence of two integers. Number of private SNPs to
#'  simulate per phenotype pair.
#' @param D Integer >= 2. Number of phenotypes to simulate.
#' @param rho Float in [-1, 1]. Correlation of each pair of phenotypes.
generate_beta <- function(M_s, M_p, D, rho) {
  beta <- generate_beta_pair(M_s, M_p, rho)
  if(D > 3){
    for( i in 1:(floor(D/2)-1) ){
      beta_pair <- generate_beta_pair(M_s, M_p, rho)
      beta <- rbind(cbind(beta, matrix(0L, nrow = nrow(beta), ncol = 2)),
                    cbind(matrix(0L, nrow = nrow(beta_pair), ncol = ncol(beta)), beta_pair))
    }
  }
  if(D%%2==1){
    beta <- rbind(cbind(beta, matrix(0L, nrow = nrow(beta), ncol = 1)),
                  cbind(matrix(0L, nrow = M_p[1], ncol = ncol(beta)), stats::rnorm(M_p[1])))
  }
  colnames(beta) <- paste0("P", 1:ncol(beta))
  rownames(beta) <- paste0("rs", 1:nrow(beta))
  return(beta)
}


#' Generates and optionally normalizes causal graph.
#'
#' Generates a DxD sparse matrix with normally distributed edge weights.
#' Optionally normalize the matrix to have eigenvalues with modulus between
#' -1 and 1.
#'
#' @param D Integer. Number of phenotypes to simulate.
generate_network <- function(D, graph, prob, g, v, orient = "random"){
  R <- huge.generator(n = D, d = D, graph = graph, prob = prob, g = g, v = v)$omega
  diag(R) = 0
  adj <- matrix(abs(R) > 1e-8, nrow = D)
  node_degree <- rowSums(adj)
  for( i in 1:(D-1) ){
    for(j in (i+1):D){
      if(R[i, j] > 0){
        if(orient == "random"){
          if(runif(1) > 0.5){
            R[i, j] = 0
          } else{
            R[j, i] = 0
          }
        } else if(orient == "away"){
          if(runif(1) > node_degree[i]/(node_degree[i] + node_degree[j])){
            R[i, j] = 0
          } else{
            R[j, i] = 0
          }
        } else if(orient == "towards"){
          if(runif(1) > node_degree[j]/(node_degree[i] + node_degree[j])){
            R[i, j] = 0
          } else{
            R[j, i] = 0
          }
        }
      }
    }
  }

  values <- eigen(R, symmetric = FALSE, only.values = TRUE)$values
  max_value <- max(abs(values))
  if (max_value > 1.0) {
      R <- R / (max_value + 0.1)
  }
  return(R)
}

#' Generates sumstats.
#'
#' @param beta M x D matrix of effect sizes.
#' @param M_null Integer. Adds M_null additional SNPs with null effects.
#' @param N Integer or array of integers of length D. Number of samples for each
#'   phenotype.
generate_sumstats <- function(beta, M_null, N){
  M <- nrow(beta)
  D <- ncol(beta)
  beta_eff <- t(matrix(stats::rnorm(length(beta), mean = t(beta), sd = 1/sqrt(N)), nrow = D))
  beta_null <- t(matrix(stats::rnorm(M_null * D, mean = 0, sd = 1/sqrt(N)), nrow = D))
  beta_hat <- rbind(beta_eff, beta_null)
  se_hat <- t(matrix(rep(1/sqrt(N), D/length(N)*(M + M_null)), nrow = D))

  colnames(beta_hat) <- paste0("P", 1:ncol(beta_hat))
  rownames(beta_hat) <- paste0("rs", 1:nrow(beta_hat))
  colnames(se_hat) <- paste0("P", 1:ncol(se_hat))
  rownames(se_hat) <- paste0("rs", 1:nrow(se_hat))
  return(list("beta_hat" = data.frame(beta_hat), "se_hat" = data.frame(se_hat)))
}

#' Generates dataset.
#'
#' Generates a network mendelian randomization dataset with observed phenotypes,
#' genotypes, true network effects and true effect sizes.
#'
#' @param N Integer or array of integers of length D. Number of individuals to simulate.
#' @param D Integer. Number of phenotypes to simulate.
#' @param M_total Integer. Total number of SNPs to simulate.
#' @param M_s Integer. Number of shared SNPs to simulate per pair.
#' @param M_p Integer or sequence of two integers. Number of private SNPs to
#'   simulate per phenotype pair.
#' @param prop_shared NULL, float between 0 and 1 or sequence thereof of length D.
#'   proportion of the variance explained by shared SNPs.
#' @param rho Float. Correlation of shared SNPs.
#' @param p_net Float. Proportion of non-zero network edges.
#' @param sd_net Float. Standard deviation of network edge weights.
#' @param noise Float. Proportion of variance explained by environmental and
#'   measurement noise.
#' @param fix_R Null or DxD matrix. Use to provide R in order to repeatedly
#'   resample from the same model.
#' @param fix_beta Null or MxD matrix. Use to provide beta in order to
#'   repeatedly resample from the same model.
generate_dataset <- function(N, D, M_total, M_s, M_p, prop_shared, rho, noise, p_net = NULL, sd_net = NULL,
                             fix_R = NULL, fix_beta = NULL) {
  if(M_total < D*M_p + floor(D/2)*M_s){
    stop("Total number of SNPs less than combined shared and private SNPs.")
  }

  # Generate the various components.
  beta <- fix_beta
  if (is.null(beta)) {
    beta <- generate_beta(M_s = M_s, M_p = M_p, D = D, rho = rho)
  }

  R <- fix_R
  if (is.null(R)) {
    stop(":-(")
    # R <- generate_network(D = D, p = p_net, sd = sd_net)
  }

  mix_mat <- solve(diag(D) - R)
  beta_obs <- beta %*% mix_mat

  if(!is.null(prop_shared)){
    for(i in seq(1, 2*floor(D/2), 2)){
      prop_1 = prop_shared[i]
      prop_2 = prop_shared[i+1]
      eff1 <- (abs(beta[, i]) > 0)
      eff2 <- (abs(beta[, i+1]) > 0)
      shared <- eff1 & eff2
      private_1 = eff1 & !eff2
      private_2 = eff2 & !eff1

      M_s1_eff <- sum(beta_obs[shared, i]**2)
      M_s2_eff <- sum(beta_obs[shared, i+1]**2)
      M_p1_eff <- sum(beta_obs[private_1, i]**2)
      M_p2_eff <- sum(beta_obs[private_2, i+1]**2)
      v1s = M_s1_eff/(M_s1_eff + M_p1_eff)
      v1p = M_p1_eff/(M_s1_eff + M_p1_eff)
      v2s = M_s2_eff/(M_s2_eff + M_p2_eff)
      v2p = M_p2_eff/(M_s2_eff + M_p2_eff)

      beta_obs[shared, i] = sqrt(prop_1/v1s)*beta_obs[shared, i]
      beta_obs[private_1, i] = sqrt((1-prop_1)/v1p)*beta_obs[private_1, i]
      beta_obs[shared, i+1] = sqrt(prop_2/v2s)*beta_obs[shared, i+1]
      beta_obs[private_2, i+1] = sqrt((1-prop_2)/v2p)*beta_obs[private_2, i+1]
    }
  }

  Y_var <- colSums(beta_obs**2)
  beta_scale <- t(t(beta_obs) / sqrt(Y_var)) * sqrt(1-noise)

  colnames(beta_scale) = colnames(beta)
  sumstats_select <- generate_sumstats(beta_scale, M_total - nrow(beta), N)
  sumstats_fit <- generate_sumstats(beta_scale, M_total - nrow(beta), N)

  R_tce = get_tce(get_observed(R), normalize=sqrt(Y_var))
  R_normed <- fit_exact(R_tce)$R_hat
  return(list("sumstats_select" = sumstats_select,
              "sumstats_fit" = sumstats_fit,
              "R_cde" = R_normed,
              "R_tce" = R_tce,
              "beta" = beta))
}


#' Selects SNPs based on true effect sizes.
#'
#' @param beta M x D matrix of true effect sizes.
select_snps_oracle <- function(beta) {
  non_pleio_snps <- (rowSums(abs(beta) > 0) == 1)
  selected <- purrr::imap(as.data.frame(beta), function(b, name){
    mask <- (abs(b) > 0) & (non_pleio_snps)
    res <- purrr::map(colnames(beta), function(x) {
      if(x == name){
        return(rep(FALSE, sum(mask)))
      } else{
        return(rep(TRUE, sum(mask)))
      }
    })
    names(res) <- colnames(beta)
    res <- purrr::prepend(res, list("names" = rownames(beta)[mask]))
    return(res)
  })
  return(selected)
}


#' A simple helper function to count the number of instuments for each pair.
#'
#' @param selected A list of lists, output of `select_snps`.
count_instruments <- function(selected) {
  return(do.call(cbind, purrr::map(selected, function(pheno) {
    return(purrr::map(pheno[-1], function(x){sum(x>0)}))
  })))
}

# TODO(brielin): Not updated for pleitropic SNPs, fix this if we want to
# revisit these metrics.
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

#' Calculates metrics for model evaluation.
#'
#' @param X DxD matrix of predicted parameters.
#' @param X_true DxD matrix of true parameters
#' @param eps float. Absolute values below eps will be considered 0.
calc_metrics <- function(X, X_true, eps = 1e-10) {
  D <- ncol(X)
  X <- X[!diag(D)]
  X_true <- X_true[!diag(D)]
  X[abs(X) < eps] <- 0
  X_true[abs(X_true) < eps] <- 0

  rmse <- sqrt(mean((X - X_true)^2))
  mae <- mean(abs(X - X_true))

  sign_X <- sign(X)
  sign_Xt <- sign(X_true)
  TN <- sum(!(abs(sign_X) | abs(sign_Xt)))
  TS <- sum(sign_X == sign_Xt) - TN
  FN <- sum((1 - abs(sign_X)) & abs(sign_Xt))
  FS <- sum(sign_X != sign_Xt) - FN

  acc <- mean(sign_X == sign_Xt)
  N_pos <- sum(X_true > 0)
  N_neg <- sum(X_true < 0)
  N_zero <- sum(X_true == 0)
  weights <- sign_Xt
  weights[sign_Xt > 0] <- 1 / N_pos
  weights[sign_Xt < 0] <- 1 / N_neg
  weights[sign_Xt == 0] <- 1 / N_zero
  weight_acc <- sum((sign_X == sign_Xt) * weights) / sum(weights)
  precision = TS / (TS + FS)
  recall = TS / (TS + FN)
  F1 <- 2*precision*recall/(precision + recall)
  F1 <- dplyr::if_else(!is.na(F1), F1, 0.0)

  return(list("precision" = precision, "recall" = recall,
              "F1" = F1,
              "rmse" = rmse, "mae" = mae, "acc" = acc,
              "weight_acc" = weight_acc))
}
