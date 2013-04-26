#' Convex Clustering Path via AMA
#' 
#' \code{cvxclust_path_ama} estimates the convex clustering path via the Alternating Minimization Algorithm.
#' Required inputs include a data matrix \code{X} (rows are features; columns are samples), a vector of weights
#' \code{w}, and a sequence of regularization parameters \code{gamma}.
#' Two penalty norms are currently supported: 1-norm and 2-norm.
#' AMA is performing proximal gradient ascent on the dual function, and therefore can be accelerated with FISTA.
#' Additional speed up can be had by employing backtracking. Both speed-ups are employed by default.
#' 
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param w A vector of nonnegative weights. The ith entry w[i] denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param gamma A sequence of regularization parameters.
#' @param nu The initial step size parameter when backtracking is applied. Otherwise it is a fixed step size in which case there are no guarantees of convergence if it exceeds \code{2/ncol(X)}.
#' @param tol The convergence tolerance.
#' @param max_iter The maximum number of iterations.
#' @param type An integer indicating the norm used: 1 = 1-norm, 2 = 2-norm.
#' @param accelerate If \code{TRUE} (the default), acceleration is turned on.
#' @param backtracking If \code{TRUE} (the default), backtracking is used.
#' @export
#' @seealso \code{\link{cvxclust_path_admm}} for estimating the clustering path with ADMM. \code{\link{kernel_weights}} and \code{\link{knn_weights}} compute useful weights.
#' @examples
#' ## Create a small set of points to cluster.
#' set.seed(12345)
#' p = 10
#' q = 2
#' X = matrix(rnorm(p*q),q,p)
#' 
#' ## Pick some weights and a sequence of regularization parameters.
#' w = kernel_weights(X,2)
#' gamma = c()
#' gamma[1] = 0
#' gamma[2] = 1.01
#' i = 3
#' repeat {
#'   g = 1.05*gamma[i-1]
#'   if (g > 50) break
#'   gamma[i] = g
#'   i = i + 1
#' }
#' nGamma = i-1
#' gamma = gamma[1:nGamma]
#' 
#' ## Perform clustering
#' sol = cvxclust_path_ama(X,w,gamma)
#' 
#' ## Plot the cluster path
#' library(ggplot2)
#' df.paths = data.frame(x=c(),y=c(), group=c())
#' for (j in 1:nGamma) {
#'   x = sol$U[[j]][1,]
#'   y = sol$U[[j]][2,]
#'   df = data.frame(x=x, y=y, group=1:p)
#'   df.paths = rbind(df.paths,df)
#' }
#' data_plot = ggplot(data=df.paths,aes(x=x,y=y))
#' data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
#' data_plot = data_plot + geom_point(data=data.frame(x=X[1,],y=X[2,]),aes(x=x,y=y))
#' data_plot + theme_bw()
cvxclust_path_ama = function(X,w,gamma,nu=10,tol=1e-4,max_iter=1e4,type=2,accelerate=TRUE,
                             backtracking=TRUE) {
  call = match.call()
  if (!is.null(type) && !(type %in% c(1,2)))
    stop("type must be 1, 2, or NULL. Only 1-norm and 2-norm penalties are currently supported.")
  nGamma = length(gamma)
  p = ncol(X)
  q = nrow(X)
  edge_info = compactify_edges(w,p)
  nK = length(which(w > 0))
  Lambda = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  ix = edge_info$ix  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  iter_vec = integer(nGamma)
  print("gamma    its | primal obj       dual obj        gap      ")
  print("---------------------------------------------------------")  
  for (ig in 1:nGamma) {
    gam = gamma[ig]
    cc = cvxclust_ama(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gam,nu,max_iter=1e4,tol=tol,type=type,
                      accelerate=accelerate,backtracking=backtracking)
    iter_vec[ig] = cc$iter
    nu = cc$nu
    Lambda = cc$Lambda
    list_U[[ig]] = cc$U
    list_V[[ig]] = cc$V
    list_Lambda[[ig]] = Lambda 
    print(sprintf("%5d  %5d | %5f        %5f      %7f", ig, cc$iter, signif(cc$primal[cc$iter],4),
                  signif(cc$dual[cc$iter],4),
                  signif(cc$primal[cc$iter]-cc$dual[cc$iter],4)))
    #        if (norm(cc$V,'1')==0) {
    #          print('Single cluster')
    #          break
    #        }
  }
  cvxclust_obj <- list(U=list_U,V=list_V,Lambda=list_Lambda,nGamma=ig,iters=iter_vec,call=call)
  class(cvxclust_obj) <- "cvxclustobject"
  return(cvxclust_obj)
}
