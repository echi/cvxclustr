#' Convex Clustering via ADMM
#' 
#' \code{cvxclust_admm} performs convex clustering via ADMM. This is an R wrapper function around Fortran code.
#' Dimensions of various arguments are as follows:
#' \itemize{
#' \item{q is the number of data points}
#' \item{p is the number of features}
#' \item{k is the number non-zero weights.}
#' }
#' 
#' @param X The q-by-p data matrix whose columns are to be clustered.
#' @param Lambda The q-by-k matrix of Lagrange multipliers.
#' @param V The q-by-k matrix of centroid differences.
#' @param ix The k-by-2 matrix of index pairs.
#' @param w The vector of k positive weights.
#' @param gamma The regularization parameter controlling the amount of shrinkage.
#' @param nu A positive penalty parameter for quadratic deviation term.
#' @param tol The convergence tolerance.
#' @param max_iter The maximum number of iterations.
#' @param type An integer indicating the norm used: 1 = 1-norm, 2 = 2-norm.
#' @param accelerate If \code{TRUE} (the default), acceleration is turned on.
#' @export
#' @author Eric C. Chi
#' @useDynLib cvxclustr
cvxclust_admm = function(X,Lambda,V,ix,w,gamma,nu=1,type=2,max_iter=1e2,tol=1e-4,accelerate=TRUE) {

  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  if (!is.null(type) && !(type %in% c(1,2)))
    stop("type must be 1, 2, or NULL. Only 1-norm and 2-norm penalties are currently supported.")
  nK = as.integer(ncol(Lambda))
  gamma = as.double(gamma)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  storage.mode(X) = "double"
  w = as.double(w)
  storage.mode(Lambda) = "double"
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  nu = as.double(nu)
  type = as.integer(type)
  primal = double(max_iter)
  dual = double(max_iter)
  if (accelerate) {
    fxname = 'convex_cluster_admm_acc'
  } else {
    fxname = 'convex_cluster_admm'
  }
  sol = .Fortran(fxname,X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 primal=primal,dual=dual,max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}
