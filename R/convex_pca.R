#' Calculates weights to use for convex PCA from observed SEs.
#'
#' This uses standard inverse variance weighting for 0 for missing entries.
#' If there are entries with SE 0, they will receive the same weight as the
#' largest weighted non-zero entry, rather than Infinity.
#'
#' @param SE. A matrix of floats, with missing entries (NA) allowed.
make_weights <- function(SE, max_weight = NULL){
  weights <- 1/SE^2
  weights[is.na(SE)] <- 0

  if(is.null(max_weight)){
    infs <- is.infinite(weights)
    weights[infs] <- 0
    max_weight = max(weights)
    weights[infs] <- max_weight
  }

  weights[weights > max_weight] = max_weight
  return(weights)
}

#' Worker function for convex PCA. See `convex_pca` for more.
#'
#' @param Y Matrix. See `convex_pca`.
#' @param r Float. Maximum nuclear norm of the returned solution.
#' @param weights Matrix. See `convex_pca`.
#' @param rmse_target Float. See `convex_pca`.
#' @param warm_start Matrix. See `convex_pca`
#' @param verbose Boolean. See `convex_pca`
convex_pca_worker  <- function(Y, r=NULL, weights=NULL, its=1000, rmse_target=1e-4, warm_start=NULL, verbose=FALSE){
  P=nrow(Y)
  N=ncol(Y)
  X=if (is.null(warm_start)) matrix(0,P,N) else warm_start
  if (is.null(weights)) weights=!is.na(Y)
  class(weights)="numeric"

  oldRmse=Inf
  for (it in 0:its){
    tol=max( 1e-1/(it+1)^2, 1e-6 )
    g=weights * (X-Y)
    g[weights == 0] = 0  # Remove NA.
    svdTemp=irlba::irlba(-g,1,1,tol=tol)
    ruv=(r * svdTemp$u) %*% t(svdTemp$v)
    erv=(X-ruv) * weights
    stepSize=sum(g*erv)/sum(erv*erv)
    if (stepSize<0)
      cat('Warning: step size is',stepSize,'\n')
    stepSize=min(stepSize,1)
    X = (1.0-stepSize)*X + stepSize* ruv
    er=weights * (X-Y)
    er[weights == 0] = 0  # Remove NA.

    rmse=sqrt(sum(er^2)/sum(weights))
    rmseDelta=abs(rmse-oldRmse)
    if (verbose)
      cat(it,rmse,stepSize,rmseDelta,'\n')
    if ( rmseDelta < rmse_target)
      break
    oldRmse=rmse
  }
  X
}

#' Cross-validation wrapper for convex PCA worker.
#'
#' @param Y See `convex_pca`.
#' @param rs A list of r values to try.
#' @param Weights See `convex_pca`.
#' @param its See `convex_pca`.
#' @param rmse_target See `convex_pca`.
find_best_r <- function(Y, rs, weights=NULL,its=1000,rmse_target=1e-3, verbose=FALSE) {
  if (is.null(weights)) weights=!is.na(Y)
  class(weights)="numeric"
  # break the matrix into training and test sets, equally and at random
  rand=matrix( runif(nrow(Y)*ncol(Y)),nrow(Y))
  train=rand < 4/5
  test= !train & (weights>0)

  # training data
  train_weights=weights
  train_weights[!train]=0

  test_weights=weights
  test_weights[!test]=0

  # try a range of regularisation parameters
  trainErrors=numeric(length(rs))
  testErrors=numeric(length(rs))
  times=numeric(length(rs))
  for (ri in 1:length(rs)){
    r=rs[ri]
    if(verbose) sprintf("Trying regularisation parameter r = %f", r)
    if(ri > 1){
      warm_start = X
    } else {
      warm_start = NULL
    }
    times[ri]=system.time( X<-convex_pca_worker(Y, r, weights=train_weights,its=its,rmse_target=rmse_target,warm_start=warm_start ))[1]
    trainErrors[ri]=sqrt(mean( train_weights * (X-Y)^2,na.rm=TRUE))
    testErrors[ri]=sqrt(mean( test_weights * (X-Y)^2, na.rm=TRUE))
    if(verbose) sprintf("Error is %f",testErrors[ri])
  }
  rs[which.min(testErrors)]
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
#' @param r Float, sequence of floats, or NULL. Maximum nuclear norm of the returned solution.
#' @param n_components Integer or `NULL`. Number of approximate principal components to keep.
#'   If `NULL`, the entire matrix is returned.
#' @param weights Matrix with weights for the loss. If `NULL`, will default
#'   to 1 for every observed value and 0 for every missing value.
#' @param its Integer. Maximum number of iterations to run.
#' @param rmse_target Float. Stop when rmse change in solution is less than `rmse_target`.
#' @param warm_start Matrix or `NULL`. If non-`NULL`, use this matrix as initial solution.
#' @param verbose Boolean. `TRUE` to print convergance progress.
convex_pca <- function(Y, r=NULL, weights=NULL, n_components = NULL, its=1000, rmse_target=1e-4, verbose=FALSE){
  rs <- NULL
  if(!is.null(r)){
    if(length(r) > 1){
      rs <- r
    }
  }
  else{
    rs <- 10^seq(0,log10(sum(Y^2,na.rm=TRUE)+1),length.out = 20)
  }

  if(!is.null(rs)){
    r <- find_best_r(Y, rs, weights = weights, its = its, rmse_target = rmse_target)
    if(verbose) sprintf("Using regulariation parameter %f", r)
  }
  X <- convex_pca_worker(Y, r, weights = weights, its = its, rmse_target = rmse_target, verbose = verbose)
  if(!is.null(n_components)){
    irlba::irlba(X, n_components, n_components)
  } else{
    X
  }
}
