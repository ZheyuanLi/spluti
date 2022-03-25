\name{solve.PP}

\alias{solve.PP}
\alias{solve.PPA}

\title{
"PP" and "PPA" methods for generic function solve()
}

\description{
Solve an inverse problem \eqn{f^{(m)}(x) = b}, where \eqn{f(x)} is a "PP" or "PPA" object.
}

\usage{
\method{solve}{PP}(a, b = 0, deriv = 0, \dots)

\method{solve}{PPA}(a, b = 0, deriv = 0, \dots)
}

\arguments{
  \item{a}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
  \item{b}{The \eqn{b} value in an inverse problem.}
  \item{deriv}{The \eqn{m} value, i.e., the number of derivatives to take in an inverse problem. Must be strictly smaller than the polynomial degree.}
  \item{\dots}{Not used by the method.}
}

\value{
x-coordinate values that solve the given inverse problem, either a vector if a "PP" is solved, or a list (an entry per replicate) of vectors if a "PPA" is solved. Note that an empty vector \code{numeric(0)} implies no solution.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}
