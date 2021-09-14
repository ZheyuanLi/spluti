\name{summary.PP}

\alias{summary.PP}
\alias{summary.PPA}

\title{
"PP" and "PPA" methods for generic function summary()
}

\description{
Print polynomial equations as formatted strings for all polynomial pieces in a "PP" or "PPA" object.
}

\usage{
\method{summary}{PP}(object, nhead = NULL, \dots)

\method{summary}{PPA}(object, nhead = NULL, \dots)
}

\arguments{
  \item{object}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
  \item{nhead}{Number of piecewise polynomial to export. Normally not used by users (hence left as \code{NULL} and later set as the number of piecewise polynomial). A small value of 6 is internally used by \code{\link{print.PP}} for a preview.}
  \item{\dots}{Not used by the method.}
}

\value{
Formatted polynomial equations, either a character vector if a "PP" object is summarized, or a list (an entry per replicate) of character vectors if a "PPA" object is summarized.
}

\note{
The polynomial equations always use "x" as variable.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

## Here I simply use Poly2PP() to create example "PP" and "PPA" objects

#################
## "PP" object ##
#################

PP <- Poly2PP(pc = runif(4), knots = 1:7)

## view piecewise polynomial coefficients (a matrix)
PP$coef

## summarize the object
summary(PP)

##################
## "PPA" object ##
##################

## 2 replicates
PPA <- Poly2PP(pc = matrix(runif(4 * 2), nrow = 4), knots = 1:7)

## view piecewise polynomial coefficients (now a 3-dimensional array)
PPA$coef

## summarize the object
summary(PPA)
}
