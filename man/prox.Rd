\name{prox}
\alias{prox}
\title{Proximal Mapping}
\usage{
  prox(x, tau, type = 2)
}
\arguments{
  \item{x}{input vector}

  \item{tau}{positive double scaling factor}

  \item{type}{integer to select proximal mapping type}
}
\description{
  \code{prox} computes the proximal mapping of norms and
  indicator functions over norm balls.
}
\examples{
set.seed(12345)
x = rnorm(10)
prox(x,1)
}
\author{
  Eric C. Chi
}

