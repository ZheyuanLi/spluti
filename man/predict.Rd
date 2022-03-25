\name{predict.PP}

\alias{predict.PP}

\title{
"PP" and "PPA" methods for generic function predict()
}

\description{
Evaluate a "PP" or "PPA" object, or its derivatives, at a vector of x-coordinate values.
}

\usage{
\method{predict}{PP}(object, newx, deriv = 0, \dots)
}

\arguments{
  \item{object}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
  \item{newx}{A vector of x-coordinate values where the piecewise polynomial or its derivatives is evaluated. All values must be within the range of knots in \code{object}; out-of-boundary prediction is not allowed.}
  \item{deriv}{Number of derivatives to take. Not allowed to exceed polynomial degree, otherwise the differentiated piecewise polynomial becomes 0.}
  \item{\dots}{Not used by the method.}
}

\value{
Evaluated y-coordinate values, either a vector if a "PP" object is predicted, or a matrix (a column per replicate) if a "PPA" object is predicted.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}
