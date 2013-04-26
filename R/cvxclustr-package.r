#' Convex Clustering via Splitting Methods
#' 
#' Clustering is a fundamental problem in science and engineering. Many classic methods such as $k$-means,
#' Gaussian mixture models, and hierarchical clustering, however, employ greedy algorithms which can be
#' entrapped in local minima, sometimes drastical suboptimal ones at that. Recently introduced convex relaxations
#' of $k$-means and hierarchical clustering shrink cluster centroids toward one another and ensure a unique global minimizer. 
#' This package provides two variable splitting methods
#' \itemize{
#' \item{Alternating Method of Multipliers (ADMM)}
#' \item{Alternating Minimization Algorithm (AMA)}
#' }
#' for solving this convex formulation of the clustering problem. We seek the centroids u_i that minimize
#' \deqn{
#' \frac{1}{2} \sum_i || x_i - u_i||_2^2 + \gamma \sum_l w_{l} ||u_{l1} - u_{l2} ||
#' }
#' Two penalty norms are currently supported: 1-norm and 2-norm.
#'
#' The two main functions are \code{\link{cvxclust_path_admm}} and \code{\link{cvxclust_path_ama}} which compute the cluster paths using
#' the ADMM and AMA methods respectively. The function \code{\link{cvxclust}} is a wrapper function that calls either 
#' \code{cvxclust_path_admm} or \code{cvxclust_path_ama} (the default) to perform the computation.
#' 
#' The functions \code{\link{kernel_weights}} and \code{\link{knn_weights}} can be used in sequence
#' to compute weights that can improve the quality of the clustering paths.
#' 
#' The typical usage consists of three steps
#' \itemize{
#' \item Compute weights \code{w}
#' \item Generate a geometrically increasing regularization parameter sequence. Unfortunately a closed form expression for the minimum amount of penalization to get complete coalescence is currently unknown.
#' \item Call \code{\link{cvxclust}} using the data \code{X}, weights \code{w}, and regularization parameter sequence \code{gamma}.
#' }
#'
#' @author Eric C. Chi
#' @references Eric C. Chi, Kenneth Lange (2013). Splitting Methods
#'   for convex clustering. arXiv: 1304.0499 [stat.ML].
#'   \url{http://arxiv.org/abs/1304.0499}.
#' @name cvxclustr
#' @docType package
NULL
