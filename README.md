#cvxclustr

Clustering is a fundamental problem in science and engineering. Many classic methods such as $k$-means, Gaussian mixture models, and hierarchical clustering, however, employ greedy algorithms which can be
entrapped in local minima, sometimes drastical suboptimal ones at that. Recently introduced convex relaxations
of $k$-means and hierarchical clustering shrink cluster centroids toward one another and ensure a unique global minimizer. 
This package provides two variable splitting methods

	* Alternating Method of Multipliers (ADMM)
	* Alternating Minimization Algorithm (AMA)

for solving this convex formulation of the clustering problem. We seek the centroids $u_i$
that minimize
$\frac{1}{2} \sum_i || x_i - u_i||_2^2 + \gamma \sum_l w_{l} ||u_{l1} - u_{l2} ||$
Two penalty norms are currently supported: 1-norm and 2-norm.

The R package `cvxclustr` implements two variable splitting methods (ADMM and AMA) for solving a convex formulation of the clustering
problem.

## Remarks

Details on the algorithms are in the paper [Splitting Methods for Convex Clustering] (http://arxiv.org/abs/1304.0499) by Chi and Lange.
