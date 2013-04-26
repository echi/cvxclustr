#' Convex Clustering via AMA
#' 
#' \code{cvxclust_ama} performs convex clustering via AMA. Let
#' q denote the number of data points and p denote the number of covariates. Let k denote
#' the number non-zero weights.
#' 
#' @param X q-by-p data matrix
#' @param Lambda q-by-k matrix of Lagrange multipliers
#' @param ix k-by-2 matrix of index pairs
#' @param M1
#' @param M2
#' @param s1
#' @param s2
#' @param w vector of k positive weights
#' @param gamma regularization parameter controlling the amount of shrinkage
#' @param eta step decrement
#' @param nu positive penalty parameter for quadratic deviation term
#' @param type integer indicating the norm used
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param accelerate boolean indicating whether to use acceleration
#' @useDynLib cvxclustr
cvxclust_ama = function(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,eta=2,type=2,max_iter=1e2,tol=1e-4,accelerate=TRUE) {
  q = as.integer(nrow(X))
  p = as.integer(ncol(X))
  nK = as.integer(ncol(Lambda))
  mix1 = as.integer(nrow(M1))
  mix2 = as.integer(nrow(M2))
  storage.mode(X) = "double"
  storage.mode(Lambda) = "double"
  U = matrix(0,q,p)
  storage.mode(U) = "double"
  V = matrix(0,q,nK)
  storage.mode(V) = "double"
  storage.mode(ix) = "integer"
  storage.mode(M1) = "integer"
  storage.mode(M2) = "integer"
  s1 = as.integer(s1)
  s2 = as.integer(s2)
  w = as.double(w)
  gamma = as.double(gamma)
  nu = as.double(nu)
  eta = as.double(eta)
  type = as.integer(type)
  max_iter = as.integer(max_iter)
  tol = as.double(tol)
  primal = double(max_iter)
  dual = double(max_iter)
  if (accelerate) {
    fxname = 'convex_cluster_ama_fista_backtrack'
#    fxname = 'convex_cluster_ama_fista'
  } else {
    fxname = 'convex_cluster_ama_backtrack'
#    fxname = 'convex_cluster_ama'
  }  
  sol = .Fortran(fxname,X=X,Lambda=Lambda,U=U,V=V,q=q,p=p,nK=nK,ix=ix,w=w,gamma=gamma,nu=nu,
                 eta=eta,s1=s1,s2=s2,M1=M1,M2=M2,mix1=mix1,mix2=mix2,primal=primal,dual=dual,
                 max_iter=max_iter,iter=integer(1),tol=tol,type=type)
  return(list(U=sol$U,V=sol$V,Lambda=sol$Lambda,nu=sol$nu,
              primal=sol$primal[1:sol$iter],dual=sol$dual[1:sol$iter],iter=sol$iter))
}
