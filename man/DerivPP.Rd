\name{DerivPP}

\alias{DerivPP}

\title{
Differentiate a "PP" or "PPA" object
}

\description{
Differentiate all polynomial pieces simultaneously in a "PP" or "PPA" object to obtain derivatives of the object.
}

\usage{
DerivPP(PPObj, deriv = 0)
}

\arguments{
  \item{PPObj}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
  \item{deriv}{Number of derivatives to take. Not allowed to exceed the polynomial degree of \code{PPObj}, otherwise the differentiated piecewise polynomial becomes 0.}
}

\value{
A "PP" or "PPA" object.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

###############
## Example 1 ##
###############

## construct a random B-spline and represent it as a "PP" object
## by default, the polynomial degree is 3
spl <- BSpl2PP(x = rnorm(10), df = 7, intercept = TRUE, coef = runif(7))
spl

## take 1st derivative, resulting in piecewise quadratic polynomials
d1_spl <- DerivPP(spl, deriv = 1)

## take 2nd derivative, resulting in piecewise linear polynomials
d2_spl <- DerivPP(spl, deriv = 2)

## take 3rd derivative, resulting in piecewise constants
d3_spl <- DerivPP(spl, deriv = 3)

###############
## Example 2 ##
###############

## construct 10 random B-splines and represent them as a "PPA" object
## by default, the polynomial degree is 3
spls <- BSpl2PP(x = rnorm(10), df = 7, intercept = TRUE, coef = matrix(runif(70), 7))
spls

## take 1st derivative, resulting in an array of piecewise quadratic polynomials
d1_spls <- DerivPP(spls, deriv = 1)
d1_spls

## take 2nd derivative, resulting in an array of piecewise linear polynomials
d2_spls <- DerivPP(spls, deriv = 2)
d2_spls

## take 3rd derivative, resulting in an array of piecewise constants
d3_spls <- DerivPP(spls, deriv = 3)
d3_spls
}
