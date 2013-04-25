#' Convex Clustering via AMA
#' 
#' \code{cvxclust_path_ama} estimates the convex clustering path via AMA using warm starts.
#' q denote the number of data points and p denote the number of covariates. Let k denote
#' the number non-zero weights.
#' 
#' @param X q-by-p data matrix
#' @param w vector of k positive weights
#' @param gamma sequence of regularization parameters
#' @param nu0 positive penalty parameter for quadratic deviation term
#' @param tol convergence tolerance
#' @param max_iter maximum number of iterations
#' @param type integer indicating the norm used
#' @param accelerate boolean indicating whether to use acceleration
#' @export
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
#'   x = sol$UHx[[j]][1,]
#'   y = sol$UHx[[j]][2,]
#'   df = data.frame(x=x, y=y, group=1:p)
#'   df.paths = rbind(df.paths,df)
#' }
#' data_plot = ggplot(data=df.paths,aes(x=x,y=y))
#' data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
#' data_plot = data_plot + geom_point(data=data.frame(x=X[1,],y=X[2,]),aes(x=x,y=y))
#' data_plot + theme_bw()
cvxclust_path_ama = function(X, w, gamma,nu0=10,tol=1e-4,max_iter=1e6,type=2,accelerate=TRUE) {
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
  nu = nu0
  ix = edge_info$ix  
  M1 = edge_info$M1
  M2 = edge_info$M2
  s1 = edge_info$s1
  s2 = edge_info$s2
  iter_vec = integer(nGamma)
  for (ig in 1:nGamma) {
    gam = gamma[ig]
    cc = cvxclust_ama(X,Lambda,ix,M1,M2,s1,s2,w[w>0],gam,nu,max_iter=1e4,tol=tol,type=type,accelerate=accelerate)
    iter_vec[ig] = cc$iter
    nu = cc$nu
    Lambda = cc$Lambda
    list_U[[ig]] = cc$U
    list_V[[ig]] = cc$V
    list_Lambda[[ig]] = Lambda 
    print(paste0("gamma: ",ig,"| itn: ", cc$iter,"| primal: ", signif(cc$primal[cc$iter],4),
                 "| dual: ", signif(cc$dual[cc$iter],4),
                 "| gap: ", signif(cc$primal[cc$iter]-cc$dual[cc$iter],4)))
    #        if (norm(cc$V,'1')==0) {
    #          print('Single cluster')
    #          break
    #        }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig,iters=iter_vec,call=call))
}
