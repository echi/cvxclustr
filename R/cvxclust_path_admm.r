#' Convex Clustering via ADMM
#' 
#' \code{cvxclust_path_admm} estimates the convex clustering path via ADMM using warm starts.
#' q denote the number of data points and p denote the number of covariates. Let k denote
#' the number non-zero weights.
#' 
#' @param X q-by-p data matrix
#' @param w vector of k positive weights
#' @param gamma sequence of regularization parameters
#' @param nu positive penalty parameter for quadratic deviation term
#' @param tol convergence tolerance
#' @param max_iter maximum number of iterations
#' @param type integer indicating the norm used
#' @param accelerate boolean indicating whether to use acceleration
#' @useDynLib cvxclustr
#' @export
#' @examples
#' # Create a small set of points to cluster.
#' set.seed(12345)
#' p = 10
#' q = 2
#' X = matrix(rnorm(n*q),q,p)
#' 
#' # Pick some weights and a sequence of regularization parameters.
#' w = kernel_weights(X,3)
#' gamma =c()
#' gamma[1] = 0
#' gamma[2] = 1.01
#' i = 3
#' repeat {
#'  gam = 1.05*gamma[i-1]
#'  if (gam > 50) {
#'    break
#'  }
#'  gamma[i] = gam
#'  i = i + 1
#' }
#' nGamma = i-1
#' gamma = gamma[1:nGamma]
#' 
#' # Perform clustering
#' sol = cvxclust_path_admm(X,w,gamma)
#' 
#' # Plot the cluster path
#' library(ggplot2)
#' df.paths = data.frame(x=c(),y=c(), group=c())
#' for (j in 1:nGamma) {
#' x = sol$UHx[[j]][1,]
#' y = sol$UHx[[j]][2,]
#' df = data.frame(x=x, y=y, group=1:p)
#' df.paths = rbind(df.paths,df)
#' }
#' data_plot = ggplot(data=df.paths,aes(x=x,y=y))
#' data_plot = data_plot + geom_path(aes(group=group),colour='grey60',alpha=0.5)
#' data_plot = data_plot + geom_point(data=data.frame(x=X[1,],y=X[2,]),aes(x=x,y=y))
#' data_plot + theme_bw()
cvxclust_path_admm = function(X, w, gamma,nu=1,tol=1e-4,max_iter=1e6,type=2,accelerate=TRUE) {
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
  
  for (ig in 1:nGamma) {
    gam = gamma[ig]
    cc = cvxclust_admm(X,Lambda,V,ix,w,gam,nu=nu,type=type,max_iter=max_iter,tol=tol,accelerate=accelerate)
    Lambda = cc$Lambda
    V = cc$V
    list_U[[ig]] = cc$U
    list_V[[ig]] = V
    list_Lambda[[ig]] = Lambda 
    print(paste0("iters: ", cc$iter,"| primal: ", signif(cc$primal[cc$iter],6),
                     "| dual: ", signif(cc$dual[cc$iter],6),
                     "| max: ", signif(max(cc$primal[cc$iter],cc$dual[cc$iter]),6)))
        print(paste("Completed gamma",ig))
    #    if (norm(cc$V,'1')==0) {
    #      print('Single cluster')
    #      break
    #    }
  }
  return(list(UHx=list_U,VHx=list_V,LambdaHx=list_Lambda,nGamma=ig))
}
