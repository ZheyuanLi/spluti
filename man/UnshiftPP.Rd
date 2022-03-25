\name{UnshiftPP}

\alias{UnshiftPP}

\title{
Change a "PP" or "PPA" object from shifted to unshifted form
}

\description{
By default, a polynomial piece in a "PP" or "PPA" object has shifted form that looks like \eqn{f_i(x) = a_{0,i} + a_{1,i}(x - x_i) + a_{2,i}(x - x_i) ^ 2 + \ldots + a_{d,i}(x - x_i) ^ d}, where \eqn{x_i} is the \eqn{i}-th knot. \code{UnshiftPP} can do some extra calculations and express the polynomial piece in unshifted form as \eqn{f_i(x) = b_{0,i} + b_{1,i}x + b_{2,i}x ^ 2 + \ldots + b_{d,i}x ^ d}.
}

\usage{
UnshiftPP(PPObj)
}

\arguments{
\item{PPObj}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
}

\value{
A "PP" or "PPA" object in unshifted form.
}

\note{
Unshifted form is numerically less stable to \code{predict}, \code{plot} and \code{solve}, thus \pkg{spluti} discourages users from using it. When given a "PP" or "PPA" object in unshifted form, most routines will stop and report that such form is not supported, except for \code{print}, \code{summary} and the auxiliary routine \code{PPAEntry}. The idea is that unshifted form should only be used for display purpose, not for any serious computations.}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

######################
## 1. a "PP" object ##
######################

## construct a random cubic B-spline and represent it as a "PP" object
spl <- BSpl2PP(x = rnorm(10), df = 7, intercept = TRUE, coef = runif(7))

## the object is in shifted form
spl$shift  ## TRUE

## print equations of the first 2 polynomial pieces
summary(spl, nhead = 2)

## now convert the object to unshifted form
spl_unshifted <- UnshiftPP(spl)
spl_unshifted$shift  ## FALSE

## print equations of the first 2 polynomial pieces
summary(spl_unshifted, nhead = 2)

#######################
## 2. a "PPA" object ##
#######################

## construct 2 random cubic B-splines and represent them as a "PPA" object
spls <- BSpl2PP(x = rnorm(10), df = 7, intercept = TRUE, coef = matrix(runif(14), 7))

## print equations of the first 2 polynomial pieces for each replicate
summary(spls, nhead = 2)

## now convert the object to unshifted form and print equations again
spls_unshifted <- UnshiftPP(spls)
summary(spls_unshifted, nhead = 2)
}
