#' Calculates weights to use for convex PCA from observed SEs.
#'
#' This uses standard inverse variance weighting for 0 for missing entries.
#' If there are entries with SE 0, they will receive the same weight as the
#' largest weighted non-zero entry, rather than Infinity.
#'
#' @param SE A matrix of floats, with missing entries (NA) allowed.
#' @param max_weight The maximum allowable weight to use. Entries with 1/SE^2
#'   above `max_weight` will get weight `max_weight`.
make_weights <- function(SE, max_weight = NULL) {
  weights <- 1 / SE^2
  weights[is.na(SE)] <- 0

  if (is.null(max_weight)) {
    infs <- is.infinite(weights)
    weights[infs] <- 0
    max_weight <- max(weights)
    weights[infs] <- max_weight
  }

  weights[weights > max_weight] <- max_weight
  return(weights)
}

#' Worker function for convex PCA. See `convex_pca` for more.
#'
#' @param Y Matrix. See `convex_pca`.
#' @param r Float. Maximum nuclear norm of the returned solution.
#' @param weights Matrix. See `convex_pca`.
#' @param its Integer. See `convex_pca`.
#' @param rmse_target Float. See `convex_pca`.
#' @param warm_start Matrix. See `convex_pca`
#' @param verbose Boolean. See `convex_pca`
convex_pca_worker <- function(Y, r = NULL, weights = NULL, its = 1000,
                              rmse_target = 1e-4, warm_start = NULL,
                              verbose = FALSE) {
  if (is.null(warm_start)) {
    X <- matrix(0, nrow(Y), ncol(Y))
  } else {
    X <- warm_start
  }
  if (is.null(weights)) {
    weights <- !is.na(Y)
  }
  class(weights) <- "numeric"

  old_rmse <- Inf
  for (it in 0:its) {
    tol <- max(1e-1 / (it + 1)^2, 1e-6)
    g <- weights * (X - Y)
    # Y may contain NA's and 0*NA = NA.
    g[weights == 0] <- 0

    svd_temp <- irlba::irlba(-g, 1, 1, tol = tol)
    ruv <- (r * svd_temp$u) %*% t(svd_temp$v)
    erv <- (X - ruv) * weights
    step_size <- sum(g * erv) / sum(erv * erv)
    if (step_size < 0) {
      warning(c("Warning: step size is ", step_size, "\n"))
    }
    step_size <- min(step_size, 1)
    X <- (1.0 - step_size) * X + step_size * ruv
    er <- weights * (X - Y)
    # Again, Y may contain NA's and 0*NA = NA.
    er[weights == 0] <- 0

    rmse <- sqrt(sum(er^2) / sum(weights))
    rmse_delta <- abs(rmse - old_rmse)
    if (verbose) cat(it, rmse, step_size, rmse_delta, "\n")
    if (rmse_delta < rmse_target) {
      break
    }
    old_rmse <- rmse
  }
  return(X)
}

#' Cross-validation wrapper for convex PCA worker.
#'
#' @param Y See `convex_pca`.
#' @param rs A list of r values to try.
#' @param weights See `convex_pca`.
#' @param its See `convex_pca`.
#' @param rmse_target See `convex_pca`.
#' @param verbose Print progress.
find_best_r <- function(Y, rs = NULL, weights = NULL, its = 1000,
                        rmse_target = 1e-3, verbose = FALSE) {
  if (is.null(weights)) {
    weights <- !is.na(Y)
  }
  class(weights) <- "numeric"

  if (is.null(rs)) {
    rs <- 10^seq(0, log10(sum(Y^2, na.rm = TRUE) + 1), length.out = 10)
  }
  # Break the matrix into training and test sets, equally and at random
  rand <- matrix(stats::runif(nrow(Y) * ncol(Y)), nrow(Y))
  train <- rand < 2 / 3
  test <- !train & (weights > 0)

  # training data
  train_weights <- weights
  train_weights[!train] <- 0

  test_weights <- weights
  test_weights[!test] <- 0

  # try a range of regularisation parameters
  train_errors <- numeric(length(rs))
  test_errors <- numeric(length(rs))
  for (ri in 1:length(rs)) {
    r <- rs[ri]
    if (verbose) {
      cat(sprintf("Trying regularisation parameter r = %f", r))
    }
    if (ri > 1) {
      warm_start <- X
    } else {
      warm_start <- NULL
    }
    X <- convex_pca_worker(Y, r,
      weights = train_weights, its = its,
      rmse_target = rmse_target, warm_start = warm_start,
      verbose = verbose
    )
    train_errors[ri] <- sqrt(mean(train_weights * (X - Y)^2, na.rm = TRUE))
    test_errors[ri] <- sqrt(mean(test_weights * (X - Y)^2, na.rm = TRUE))
    if (verbose) {
      cat(sprintf("Error is %f", test_errors[ri]))
    }
  }
  return(rs[which.min(test_errors)])
}

#' Calculates the convex relaxation of principal component analysis
#'
#' The convex relaxation of PCA is the optimization problem
#'   || W(X - Y) ||_F ^ 2 subject to ||X||_* < r
#' Where W is a weight matrix and ||X||_* is the nuclear (trace) norm of X.
#'
#' This can be ued to generate a nuclear-norm regularizaed approximation
#' to Y, even if Y has missing values.
#'
#' @param Y Matrix. Missing values are allowed.
#' @param r Float, sequence of floats, or NULL. Maximum nuclear norm of the
#'   returned solution.
#' @param n_components Integer or `NULL`. Number of approximate principal
#'   components to keep. If `NULL`, the entire matrix is returned.
#' @param weights Matrix with weights for the loss. If `NULL`, will default
#'   to 1 for every observed value and 0 for every missing value.
#' @param its Integer. Maximum number of iterations to run.
#' @param rmse_target Float. Stop when rmse change in solution is less than
#'   `rmse_target`.
#' @param verbose Boolean. `TRUE` to print convergance progress.
convex_pca <- function(Y, r = NULL, weights = NULL, n_components = NULL,
                       its = 1000, rmse_target = 1e-4, verbose = FALSE) {
  rs <- NULL
  if (!is.null(r)) {
    if (length(r) > 1) {
      rs <- r
    }
  }
  else {
    rs <- 10^seq(0, log10(sum(Y^2, na.rm = TRUE) + 1), length.out = 10)
  }

  if (!is.null(rs)) {
    r <- find_best_r(Y, rs,
      weights = weights, its = its,
      rmse_target = rmse_target, verbose = verbose
    )
    if (verbose) {
      print(sprintf("Using regulariation parameter %f", r))
    }
  }
  X <- convex_pca_worker(Y, r,
    weights = weights, its = its,
    rmse_target = rmse_target, verbose = verbose
  )
  if (!is.null(n_components)) {
    return(irlba::irlba(X, n_components, n_components))
  } else {
    return(X)
  }
}
