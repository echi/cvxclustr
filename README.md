#cvxclustr

Clustering is a fundamental problem in science and engineering. Many classic methods such as $k$-means, Gaussian mixture models, and hierarchical clustering, however, employ greedy algorithms which can be
entrapped in local minima, sometimes drastical suboptimal ones at that. Recently introduced convex relaxations
of $k$-means and hierarchical clustering shrink cluster centroids toward one another and ensure a unique global minimizer. 
This R package provides two variable splitting methods

	* Alternating Method of Multipliers (ADMM)
	* Alternating Minimization Algorithm (AMA)

for solving this convex formulation of the clustering problem. We seek the centroids $u_i$
that minimize

$\frac{1}{2} \sum_i || x_i - u_i||_2^2 + \gamma \sum_l w_{l} ||u_{l1} - u_{l2} ||$

Two penalty norms are currently supported: 1-norm and 2-norm.

## Remarks

Details on the algorithms are in the paper [Splitting Methods for Convex Clustering] (http://arxiv.org/abs/1304.0499) by Chi and Lange.

Previous takes on this formulation of clustering:

[Just Relax and Come Clustering! A Convexification of k-Means Clustering](http://www.control.isy.liu.se/research/reports/2011/2992.pdf)
by Lindsten, Ohlsson, and Ljung.

[Clusterpath: An Algorithm for Clustering using Convex Fusion Penalties](http://www.icml-2011.org/papers/419_icmlpaper.pdf)
	by Hocking, Joulin, Bach, and Vert.
