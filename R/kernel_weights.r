#' Compute Gaussian Kernel Weights
#' 
#' \code{kernel_weights} computes Gaussian kernel weights
#' 
#' @param X q-by-p data matrix
#' @param mu nonnegative parameter that controls the scale of kernel weights
#' @export
kernel_weights = function(X,mu) {
  p = ncol(X)
  k = 1
  w = matrix(0,p*(p-1)/2,1)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      w[k] = exp(-mu*norm(as.matrix(X[,i,drop=FALSE]-X[,j,drop=FALSE]),'f')^2)
      k = k+1
    }
  }
  return(weights=w)
}
