#' Proximal Mapping
#'
#' \code{prox} computes the proximal mapping of norms and indicator functions
#' over norm balls.
#' @param x input vector
#' @param tau positive double scaling factor
#' @param type integer to select proximal mapping type
#' @useDynLib cvxclustr
#' @export
#' @author Eric C. Chi
#' @examples
#' set.seed(12345)
#' x = rnorm(10)
#' prox(x,1)
prox = function(x,tau,type=2) {
  x = as.double(x)
  tau = as.double(tau)
  type = as.integer(type)
  n = as.integer(length(x))
  y = double(n)
  sol=.Fortran("proxA",x=x,n=n,y=y,tau=tau,type=type)
  return(sol$y)
}
