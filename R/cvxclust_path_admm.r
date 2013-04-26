#' Convex Clustering Path via ADMM
#' 
#' \code{cvxclust_path_admm} estimates the convex clustering path via ADMM.
#' 
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param w A vector of nonnegative weights. The ith entry w[i] denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param gamma A sequence of regularization parameters.
#' @param nu A positive penalty parameter for quadratic deviation term.
#' @param tol The convergence tolerance.
#' @param max_iter The maximum number of iterations.
#' @param type An integer indicating the norm used: 1 = 1-norm, 2 = 2-norm.
#' @param accelerate If \code{TRUE} (the default), acceleration is turned on.
#' @export
#' @seealso \code{\link{cvxclust_path_ama}} for estimating the clustering path with AMA. \code{\link{kernel_weights}} and \code{\link{knn_weights}} compute useful weights.
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
#' sol = cvxclust_path_admm(X,w,gamma)
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
#' 
#' ## Iris Data
#' X = as.matrix(iris[,1:4])
#' X = t(scale(X,center=TRUE,scale=FALSE))
#' dupix = which(duplicated(t(X)))
#' X = X[,-dupix]
#' p = ncol(X)
#' 
#' ## Pick some weights and a sequence of regularization parameters.
#' kgamma = 4
#' knbs = 5
#' w = kernel_weights(X,kgamma)
#' w = knn_weights(w,knbs,p)
#' gamma = c()
#' gamma[1] = 0
#' gamma[2] = 1.01
#' i = 3
#' repeat {
#'   g = 1.05*gamma[i-1]
#'   if (g > 21) break
#'   gamma[i] = g
#'   i = i + 1
#' }
#' gamma = gamma[1:(i-1)]
#' 
#' ## Perform clustering
#' sol = cvxclust_path_admm(X,w,gamma,nu=1,tol=1e-1)
#' 
#' ## Plot the cluster path
#' library(ggplot2)
#' svdX = svd(X)
#' pc = svdX$u[,1:2,drop=FALSE]
#' pc.df = as.data.frame(t(pc)%*%X)
#' nGamma = sol$nGamma
#' df.paths = data.frame(x=c(),y=c(), group=c())
#' for (j in 1:nGamma) {
#'   pcs = t(pc)%*%sol$U[[j]]
#'   x = pcs[1,]
#'   y = pcs[2,]
#'   df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
#'   df.paths = rbind(df.paths,df)
#' }
#' X_data = as.data.frame(t(X)%*%pc)
#' colnames(X_data) = c("x","y")
#' X_data$Species = iris[-dupix,5]
#' data_plot = ggplot(data=df.paths,aes(x=x,y=y))
#' data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
#' data_plot = data_plot + geom_point(data=X_data,aes(x=x,y=y,colour=Species,shape=Species))
#' data_plot = data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
#' data_plot + theme_bw()
cvxclust_path_admm = function(X,w,gamma,nu=1,tol=1e-4,max_iter=1e4,type=2,accelerate=TRUE) {
  call = match.call()
  if (!is.null(type) && !(type %in% c(1,2)))
    stop("type must be 1, 2, or NULL. Only 1-norm and 2-norm penalties are currently supported.")
  nGamma = length(gamma)
  p = ncol(X)
  q = nrow(X)
  nK = p*(p-1)/2
  Lambda = matrix(0,q,nK)
  V = matrix(0,q,nK)
  list_U = vector(mode="list",length=nGamma)
  list_V = vector(mode="list",length=nGamma)
  list_Lambda = vector(mode="list",length=nGamma)
  ix = vec2tri(1:nK,p)  
  iter_vec = integer(nGamma)  
  print("gamma    its | primal res      dual res      max res     ")
  print("---------------------------------------------------------")    
  for (ig in 1:nGamma) {
    gam = gamma[ig]
    cc = cvxclust_admm(X,Lambda,V,ix,w,gam,nu=nu,type=type,max_iter=max_iter,tol=tol,accelerate=accelerate)
    iter_vec[ig] = cc$iter
    Lambda = cc$Lambda
    V = cc$V
    list_U[[ig]] = cc$U
    list_V[[ig]] = V
    list_Lambda[[ig]] = Lambda 
    print(sprintf("%5d  %5d | %5f        %5f      %7f", ig, cc$iter, signif(cc$primal[cc$iter],4),
                  signif(cc$dual[cc$iter],4),
                  signif(max(cc$primal[cc$iter],cc$dual[cc$iter]),4)))    
    #    if (norm(cc$V,'1')==0) {
    #      print('Single cluster')
    #      break
    #    }
  }
#  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig,call=call))
  cvxclust_obj <- list(U=list_U,V=list_V,Lambda=list_Lambda,nGamma=ig,iters=iter_vec,call=call)
  class(cvxclust_obj) <- "cvxclustobject"
  return(cvxclust_obj)  
}
