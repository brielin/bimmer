#' Fits exact model to data.
#'
#' @param R D x D matrix of "total causal effects".
#' @return D x D matrix with zero diagonal of deconvoluted direct effects.
fit_exact <- function(R) {
  D <- dim(R)[1]
  R[is.na(R)] <- 0
  R_inv <- solve(R)
  # This is really G
  R_hat <- diag(D) - R_inv / diag(R_inv)
  return(list("R_hat" = R_hat, "R_inv" = R_inv))
}


#' Fit's MVMR model to data.
#'
#' @importFrom foreach %do%
#'
#' @param sumstats_fit List with elements beta_hat and se_hat.
#' @param snps_to_use List of per-phenotype pair SNP to use.
#' @param max_iter Integer. Number of BMA iterations to run, default 1000.
#' @param verbose Boolean. True to print progress.
fit_mrbma <- function(sumstats_fit, snps_to_use, max_iter = 1000,
                      verbose = FALSE){
  snps_to_use <- purrr::transpose(snps_to_use)
  snp_names <- snps_to_use$names
  snps_to_use$names <- NULL
  pheno_names <- names(snps_to_use)

  out = NULL
  R_hat <- as.matrix(
    foreach::foreach(out = pheno_names, .combine = dplyr::bind_rows) %do% {
      if(verbose){
        print(out)
      }
      exp_snps <- dplyr::setdiff(
        unique(unlist(snp_names[names(snp_names) != out])),  snp_names[[out]])
      exp_phenos <- pheno_names[pheno_names != out]
      betaX <- as.matrix(sumstats_fit$beta_hat[exp_snps, exp_phenos])
      betaXse <- as.matrix(sumstats_fit$se_hat[exp_snps, exp_phenos])
      betaY <- sumstats_fit$beta_hat[exp_snps, out]
      names(betaY) <- rownames(betaX)
      betaYse <- sumstats_fit$se_hat[exp_snps, out]
      names(betaY) <- rownames(betaXse)

      mrinput <- methods::new(
        "mvMRInput", betaX = betaX, betaY = as.matrix(betaY),
        betaXse = betaXse, betaYse = as.matrix(betaYse))
      mr_res <- summarymvMR_SSS(
        mrinput, kmin = 1, kmax = length(pheno_names) - 1, max_iter = max_iter,
        print = verbose)@BestModel_Estimate
      names(mr_res) <- exp_phenos
      self <- c(0.0)
      names(self) = out
      return(c(self, mr_res))
    })

  rownames(R_hat) <- colnames(R_hat)
  return(list("R_hat" = R_hat))
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
  R_hat <- as.matrix(
    foreach::foreach(out = pheno_names, .combine = dplyr::bind_rows) %do% {
      exp_snps <- dplyr::setdiff(
        unique(unlist(snp_names[names(snp_names) != out])),  snp_names[[out]])
      exp_phenos <- pheno_names[pheno_names != out]
      betaX <- as.matrix(sumstats_fit$beta_hat[exp_snps, exp_phenos])
      betaXse <- as.matrix(sumstats_fit$se_hat[exp_snps, exp_phenos])
      betaY <- sumstats_fit$beta_hat[exp_snps, out]
      names(betaY) <- rownames(betaX)
      betaYse <- sumstats_fit$se_hat[exp_snps, out]
      names(betaY) <- rownames(betaXse)

      mrinput <- MendelianRandomization::mr_mvinput(
        bx = betaX, bxse = betaXse, by = betaY, byse = betaYse)
      mr_res <- as.vector(MendelianRandomization::mr_mvegger(mrinput)$Estimate)
      names(mr_res) <- exp_phenos
      self <- c(0.0)
      names(self) = out
      return(c(self, mr_res))
    })

  rownames(R_hat) <- colnames(R_hat)
  return(list("R_hat" = R_hat))
}


#' Fit's MVMR model to data.
#'
#' @importFrom foreach %do%
#'
#' @param sumstats_fit List with elements beta_hat and se_hat.
#' @param snps_to_use List of per-phenotype pair SNP to use.
fit_glmnet <- function(sumstats_fit, snps_to_use){
  snps_to_use <- purrr::transpose(snps_to_use)
  snp_names <- snps_to_use$names
  snps_to_use$names <- NULL
  pheno_names <- names(snps_to_use)

  out = NULL
  R_hat <- as.matrix(
    foreach::foreach(out = pheno_names, .combine = dplyr::bind_rows) %do% {
      exp_snps <- dplyr::setdiff(
        unique(unlist(snp_names[names(snp_names) != out])),  snp_names[[out]])
      exp_phenos <- pheno_names[pheno_names != out]
      betaX <- as.matrix(sumstats_fit$beta_hat[exp_snps, exp_phenos])
      betaXse <- as.matrix(sumstats_fit$se_hat[exp_snps, exp_phenos])
      betaY <- sumstats_fit$beta_hat[exp_snps, out]
      names(betaY) <- rownames(betaX)
      betaYse <- sumstats_fit$se_hat[exp_snps, out]
      names(betaY) <- rownames(betaXse)

      glm_coef <- stats::coef(glmnet::cv.glmnet(
        betaX, betaY, alpha = 0.5), s = "lambda.1se")
      glm_coef <- as.vector(glm_coef[2:length(glm_coef)])

      names(glm_coef) <- exp_phenos
      self <- c(0.0)
      names(self) = out
      return(c(self, glm_coef))
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
#' @param G D x D matrix of direct effects.
#' @return D x D matrix of observed effects.
get_observed <- function(G) {
  D <- dim(G)[1]
  return(solve(diag(D) - G, G))
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
  beta_s <- mvtnorm::rmvnorm(
    M_s, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  if(length(M_p) > 1){
    beta_p1 <- stats::rnorm(M_p[1])
    beta_p2 <- stats::rnorm(M_p[2])
    beta <- rbind(cbind(beta_p1, rep(0, M_p[1])),
                  cbind(rep(0, M_p[2]), beta_p2), beta_s)
  } else{
    beta_p1 <- stats::rnorm(M_p)
    beta_p2 <- stats::rnorm(M_p)
    beta <- rbind(cbind(beta_p1, rep(0, M_p)),
                  cbind(rep(0, M_p), beta_p2), beta_s)
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
                    cbind(matrix(0L, nrow = nrow(beta_pair), ncol = ncol(beta)),
                          beta_pair))
    }
  }
  if(D%%2==1){
    beta <- rbind(cbind(beta, matrix(0L, nrow = nrow(beta), ncol = 1)),
                  cbind(matrix(0L, nrow = M_p[1], ncol = ncol(beta)),
                        stats::rnorm(M_p[1])))
  }
  colnames(beta) <- paste0("P", 1:ncol(beta))
  rownames(beta) <- paste0("rs", 1:nrow(beta))
  return(beta)
}


#' Generates and optionally normalizes causal graph.
#'
#' This is a wrapper around the huge.generator function from that package. After
#' generating the adjacency matrix, this orients the sign randomly and then the
#' edges according to the specified method.
#'
#' @param D Integer. Number of phenotypes to simulate.
#' @param graph String, one of "random", "hub" or "scale-free". Specifies the
#'   kind of graph to be generated.
#' @param orient String, one of "random", "towards" or "away". Specifies the
#'   edge orientation strategy. Randomly, towards high degree nodes, or away
#'   from high degree nodes.
#' @param v Float between 0 and 1. Roughly corresponds to edge weight, see
#'   huge.generator for more information. Default 0.3.
#' @param v_max Float greater than v. Used as max value if pert=T. Default 0.9.
#' @param v_min Float greater than 0. Minimim value of pert=T. Default 0.01.
#' @param pert Boolean. TRUE to replace huge produced value with one sampled
#'   from the betaPERT distribution with mode v and max v_max.
#' @export
generate_network <- function(D, graph = "random", orient = "random", v = 0.3,
                             v_max = 0.9, v_min = 0.01, pert = F){
  R <- huge::huge.generator(n = D, d = D, graph = graph, v = v, verbose = FALSE)$omega
  diag(R) = 0
  adj <- matrix(abs(R) > 1e-8, nrow = D)
  node_degree <- rowSums(adj)
  for( i in 1:(D-1) ){
    for(j in (i+1):D){
      if(R[i, j] > 1e-10){
        if(pert) {
          R_ij <- mc2d::rpert(1, min = v_min, mode = v, max = v_max)
        } else {
          R_ij <- R[i, j]
        }
        R_ij <- dplyr::if_else(stats::runif(1) > 0.5, R_ij, -R_ij)
        if(orient == "random"){
          if(stats::runif(1) > 0.5){
            R[i, j] = 0
            R[j, i] = R_ij
          } else{
            R[j, i] = 0
            R[i, j] = R_ij
          }
        } else if(orient == "away"){
          if(stats::runif(1) > node_degree[i]/(node_degree[i] + node_degree[j])){
            R[i, j] = 0
            R[j, i] = R_ij
          } else{
            R[j, i] = 0
            R[i, j] = R_ij
          }
        } else if(orient == "towards"){
          if(stats::runif(1) > node_degree[j]/(node_degree[i] + node_degree[j])){
            R[i, j] = 0
            R[j, i] = R_ij
          } else{
            R[j, i] = 0
            R[i, j] = R_ij
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



generate_network_dsfg <- function(D, alpha=0.1, beta=0.5, gamma=0.4,
                                  delta_in=0.0, delta_out=0.2, v_mode = 0.3,
                                  v_max = 0.9, v_min = 0.1){
  A = matrix(c(0, 0, 1, 0), nrow = 2)
  while(nrow(A) < D){
    d_in <- colSums(A)
    d_out <- rowSums(A)
    choice <- stats::runif(1)
    D_curr <- nrow(A)
    if(choice < alpha){ # Edge from new to existing
      w <- sample(D_curr, 1, prob = d_in + delta_in)
      a_next <- rep(0, D_curr)
      a_next[w] <- 1
      A <- cbind(rbind(A, a_next, deparse.level = 0), rep(0, D_curr + 1))
    } else if(choice < alpha + beta){ # Edge between existing
      v <- sample(D_curr, 1, prob = d_out + delta_out)
      w <- sample(D_curr, 1, prob = d_in + delta_in)
      if(v != w) A[v, w] = 1
    } else{ # Edge from existing to new
      v <- sample(D_curr, 1, prob = d_out + delta_out)
      a_next <- rep(0, D_curr)
      a_next[v] <- 1
      A <- rbind(cbind(A, a_next, deparse.level = 0), rep(0, D_curr + 1))
    }
  }
  G = A
  G[A != 0] <- sample(c(1, -1), sum(A), replace = T, prob = c(0.6, 0.4)) * mc2d::rpert(sum(A), min = v_min, mode = v_mode, max = v_max)
  return(G)
}

#' Turns a fully observed network into a partially observed network.
#'
#' This turns a fully observed network matrix into a partially observed one
#' where TCEs of missing nodes are integrated into pseudo-DCEs for the observed
#' nodes.
#'
#' @param G D x D matrix representing direct causal graph.
#' @param p Float 0 to 1, proportion of nodes to hide.
#' @export
censor_network <- function(G, p){
  R <- get_tce(get_observed(G))
  D <- nrow(R)
  keep_cols <- sort(sample.int(D, size = round((1-p)*D)))
  R_cens <- R[keep_cols, keep_cols]
  G_cens <- fit_exact(R_cens)$R_hat
  return(G_cens)
}


#' Generates sumstats.
#'
#' @param beta M x D matrix of effect sizes.
#' @param N Integer or array of integers of length D. Number of samples for each
#'   phenotype.
generate_sumstats <- function(beta, N){
  M <- nrow(beta)
  D <- ncol(beta)
  beta_hat <- t(matrix(stats::rnorm(length(beta), mean = t(beta), sd = 1/sqrt(N)), nrow = D))
  se_hat <- t(matrix(rep(1/sqrt(N), D/length(N)*M), nrow = D))

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
#' @param G D x D matrix. Causal graph to simulate.
#' @param M Integer. Number of SNPs to simulate.
#' @param p Float between 0 and 1. Per-SNP effect probability.
#' @param h Float between 0 and 1. Phenotype heritability.
#' @param censor_prob Float between 0 and 1. Hide phenotypes to induce
#'   non-causal genetic correlation. Proportion of phenotype to hide.
#' @param beta M x D matrix of SNP effects or NULL to generate the effects based
#'   on the provided parameters.
#' @export
generate_dataset <- function(N, D, G, M, p, h, censor_prob = 0, beta = NULL) {
  if(is.null(beta)){
    beta <- purrr::map(1:D, function(.x){
      beta_mat <- matrix(0L, M, D)
      beta_mat[, .x] <- stats::rnorm(n = M, mean = 0, sd = sqrt(h/(M*p))) * stats::rbinom(n = M, size = 1, prob = p)
      return(beta_mat)
    })
    beta <- do.call(rbind, beta)
    # TODO(brielin): the below code simulates random pleiotropy, which we'll
    # have plenty of from the network already.
    # beta <- matrix(stats::rnorm(n = M*D, mean = 0, sd = sqrt(h/(M*p))) *
    #                  rbinom(n = M*D, size = 1, prob = p), nrow = M, ncol = D)
  }
  mix_mat <- solve(diag(D) - G)
  beta_obs <- beta %*% mix_mat
  if(censor_prob > 0){
    R <- get_tce(get_observed(G))
    keep_cols <- as.logical(rbinom(nrow(R), 1, 1-p))
    R_cens <- R[keep_cols, keep_cols]
    G <- fit_exact(R_cens)$R_hat
  }

  h_obs <- colSums(beta_obs**2)

  sumstats_select <- generate_sumstats(beta_obs, N)
  sumstats_fit <- generate_sumstats(beta_obs, N)
  R = get_tce(get_observed(G))

  return(list("sumstats_select" = sumstats_select,
              "sumstats_fit" = sumstats_fit,
              "R" = R,
              "beta" = beta,
              "h_obs" = h_obs))
}


#' Generates dataset.
#'
#' Generates a network mendelian randomization dataset with observed phenotypes,
#' genotypes, true network effects and true effect sizes.
#'
#' @param N Integer or array of integers of length D. Number of individuals to simulate.
#' @param D Integer. Number of phenotypes to simulate.
#' @param R D x D matrix. Causal graph to simulate.
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
#' @param fix_beta Null or MxD matrix. Use to provide beta in order to
#'   repeatedly resample from the same model.
# generate_dataset <- function(N, D, R, M_total, M_s, M_p, prop_shared, rho,
#                              noise, p_net = NULL, sd_net = NULL,
#                              fix_beta = NULL) {
#   if(M_total < D*M_p + floor(D/2)*M_s){
#     stop("Total number of SNPs less than combined shared and private SNPs.")
#   }
#
#   # Generate the various components.
#   beta <- fix_beta
#   if (is.null(beta)) {
#     beta <- generate_beta(M_s = M_s, M_p = M_p, D = D, rho = rho)
#   }
#
#   mix_mat <- solve(diag(D) - R)
#   beta_obs <- beta %*% mix_mat
#
#   if(!is.null(prop_shared)){
#     for(i in seq(1, 2*floor(D/2), 2)){
#       prop_1 = prop_shared[i]
#       prop_2 = prop_shared[i+1]
#       eff1 <- (abs(beta[, i]) > 0)
#       eff2 <- (abs(beta[, i+1]) > 0)
#       shared <- eff1 & eff2
#       private_1 = eff1 & !eff2
#       private_2 = eff2 & !eff1
#
#       M_s1_eff <- sum(beta_obs[shared, i]**2)
#       M_s2_eff <- sum(beta_obs[shared, i+1]**2)
#       M_p1_eff <- sum(beta_obs[private_1, i]**2)
#       M_p2_eff <- sum(beta_obs[private_2, i+1]**2)
#       v1s = M_s1_eff/(M_s1_eff + M_p1_eff)
#       v1p = M_p1_eff/(M_s1_eff + M_p1_eff)
#       v2s = M_s2_eff/(M_s2_eff + M_p2_eff)
#       v2p = M_p2_eff/(M_s2_eff + M_p2_eff)
#
#       beta_obs[shared, i] = sqrt(prop_1/v1s)*beta_obs[shared, i]
#       beta_obs[private_1, i] = sqrt((1-prop_1)/v1p)*beta_obs[private_1, i]
#       beta_obs[shared, i+1] = sqrt(prop_2/v2s)*beta_obs[shared, i+1]
#       beta_obs[private_2, i+1] = sqrt((1-prop_2)/v2p)*beta_obs[private_2, i+1]
#     }
#   }
#
#   Y_var <- colSums(beta_obs**2)
#   beta_scale <- t(t(beta_obs) / sqrt(Y_var)) * sqrt(1-noise)
#
#   colnames(beta_scale) = colnames(beta)
#   sumstats_select <- generate_sumstats(beta_scale, N)
#   sumstats_fit <- generate_sumstats(beta_scale, N)
#
#   R_tce = get_tce(get_observed(R), normalize=sqrt(Y_var))
#   R_normed <- fit_exact(R_tce)$R_hat
#   return(list("sumstats_select" = sumstats_select,
#               "sumstats_fit" = sumstats_fit,
#               "R_cde" = R_normed,
#               "R_tce" = R_tce,
#               "beta" = beta))
# }


# TODO(brielin): This is obsolete now that shared effects come from hidden
#  phenotypes. Stage for removal.
#' Selects SNPs based on true effect sizes.
#'
#' @param beta M x D matrix of true effect sizes.
# select_snps_oracle <- function(beta) {
#   non_pleio_snps <- (rowSums(abs(beta) > 0) == 1)
#   selected <- purrr::imap(as.data.frame(beta), function(b, name){
#     mask <- (abs(b) > 0) & (non_pleio_snps)
#     res <- purrr::map(colnames(beta), function(x) {
#       if(x == name){
#         return(rep(FALSE, sum(mask)))
#       } else{
#         return(rep(TRUE, sum(mask)))
#       }
#     })
#     names(res) <- colnames(beta)
#     res <- purrr::prepend(res, list("names" = rownames(beta)[mask]))
#     return(res)
#   })
#   return(selected)
# }


#' A simple helper function to count the number of instuments for each pair.
#'
#' @param selected A list of lists, output of `select_snps`.
count_instruments <- function(selected) {
  return(do.call(cbind, purrr::map(selected, function(pheno) {
    return(purrr::map(pheno[-1], function(x){sum(x>0, na.rm=T)}))
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
#' @export
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
