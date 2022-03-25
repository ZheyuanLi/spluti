\name{PP}

\alias{PP}
\alias{PPA}
\alias{PPAEntry}

\title{
The "PP" and "PPA" classes
}

\description{
Define the "PP" and "PPA" classes, explain their stored contents, and how to extract a "PP" object from a "PPA" object.
}

\usage{
PPAEntry(PPAObj, i)
}

\arguments{
\item{PPAObj}{A "PPA" object.}
\item{i}{Position index informing which (single) replicate to extract.}
}

\details{

\subsection{One single spline}{

Consider a spline \eqn{f(x)} of degree \eqn{d} (or order \eqn{(d + 1)}) and \eqn{k} knots \eqn{x_1, x_2, \ldots, x_k}. In its piecewise polynomial representation, the \eqn{i}-th polynomial piece \eqn{f_i(x)} on interval \eqn{[x_i, x_{i + 1}]} is parametrized in a \emph{shifted} form: \deqn{f_i(x) = a_{0,i} + a_{1,i}(x - x_i) + a_{2,i}(x - x_i) ^ 2 + \ldots + a_{d,i}(x - x_i) ^ d.} Coefficents for all pieces can then be stored as a \eqn{(d + 1) \times (k - 1)} matrix:

\tabular{cccc}{
\eqn{a_{0,1}} \tab \eqn{a_{0,1}} \tab \eqn{\ldots} \tab \eqn{a_{0,k-1}}\cr
\eqn{a_{1,1}} \tab \eqn{a_{1,1}} \tab \eqn{\ldots} \tab \eqn{a_{1,k-1}}\cr
\eqn{\vdots} \tab \eqn{\vdots} \tab \tab \eqn{\vdots}\cr
\eqn{a_{d-1,1}} \tab \eqn{a_{d-1,1}} \tab \eqn{\ldots} \tab \eqn{a_{d-1,k-1}}
} % tabular

If function \code{\link{UnshiftPP}} is used to transform the "PP" object from shifted form to unshifted form: \deqn{f_i(x) = b_{0,i} + b_{1,i}x + b_{2,i}x ^ 2 + \ldots + b_{d,i}x ^ d,} then the coefficient matrix will instead be:

\tabular{cccc}{
\eqn{b_{0,1}} \tab \eqn{b_{0,1}} \tab \eqn{\ldots} \tab \eqn{b_{0,k-1}}\cr
\eqn{b_{1,1}} \tab \eqn{b_{1,1}} \tab \eqn{\ldots} \tab \eqn{b_{1,k-1}}\cr
\eqn{\vdots} \tab \eqn{\vdots} \tab \tab \eqn{\vdots}\cr
\eqn{b_{d-1,1}} \tab \eqn{b_{d-1,1}} \tab \eqn{\ldots} \tab \eqn{b_{d-1,k-1}}
} % tabular
} % subsection: One single spline

\subsection{Several splines}{

To motivate this, consider a dataset with variables \code{x}, \code{y1}, \code{y2} and \code{y3}. In practice, \code{y1}, \code{y2}, \code{y3} and \code{y4} can be:

\itemize{
\item different variables;
\item repeated measurements of the same variable on the same subject;
\item the same variable measured on different subjects;
\item simulated data in bootstrapping.
} % itemize

Now a multi-response linear regression via:

\preformatted{lm(cbind(y1, y2, y3, y4) ~ bs(x, df = 10), data = dat)}

will produce 4 spline curves with the same knots and degree but different coefficients. Here we call them 4 \emph{replicates} of a spline. Since piecewise polynomial coefficient matrix for each replicate has the same size, they can be combined into a 3-dimensional array \code{A}, where \code{A[, , i]} gives the coefficient matrix associated with the \eqn{i}-th replicate.
} % subsection: Several splines
} % details

\value{
Piecewise polynomials can be stored in a list with 3 components:

\item{coef}{A \eqn{(d + 1) \times (k - 1)} matrix for a single spline, or a \eqn{(d + 1) \times (k - 1) \times n} 3-dimensional array for \eqn{n} replicates of a spline.}
\item{knots}{A vector of knots.}
\item{shift}{A TRUE/FALSE logical value indicating whether the piecewise polynomial is in shifted form.}

If \code{coef} is a matrix, the list has a single "PP" class and is called a "PP" object.

If \code{coef} is a 3-dimensional array, the list has classes "PPA" (primary) and "PP" (secondary). Primary class identifies an object so it is called a "PPA" object, but it \emph{inherits} the "PP" class.

If users want to extract a replicate from a "PPA" object for a "PP" object, they can use the auxiliary function \code{PPAEntry}. See \bold{Examples}.
}

\note{
Unless otherwise stated, all piecewise polynomials in \pkg{spluti} are in shifted form. Transformation routines transform splines to "PP" or "PPA" objects in shifted form; all methods and tools for piecewise polynomials take and return "PP" or "PPA" objects in shifted form. The only exception is \code{UnshiftPP}, which transforms "PP" or "PPA" objects from shifted to unshifted form. But \pkg{spluti} has minimal support for unshifted form. See \code{\link{UnshiftPP}} for more.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

## Here I simply use Poly2PP() to create example "PP" and "PPA" objects

######################
## 1. a "PP" object ##
######################

## express a random cubic polynomial in piecewise fashion
PP <- Poly2PP(pc = runif(4), knots = 1:7)

## view piecewise polynomial coefficients (a matrix)
PP$coef

## convert the object from shifted to unshifted form and note the coefficient change
PP_unshifted <- UnshiftPP(PP)
PP_unshifted$coef

#######################
## 2. a "PPA" object ##
#######################

## express two replicates of a random cubic polynomial in piecewise fashion
PPA <- Poly2PP(pc = matrix(runif(4 * 2), nrow = 4), knots = 1:7)

## view piecewise polynomial coefficients (now a 3-dimensional array)
PPA$coef

## extract the 1st replicate as a "PP" object
PP1 <- PPAEntry(PPA, 1)
PP1$coef  ## a matrix

## convert the object from shifted to unshifted form and note the coefficient change
PPA_unshifted <- UnshiftPP(PPA)
PPA_unshifted$coef

## extract the 2nd replicate as a "PP" object
PP2_unshifted <- PPAEntry(PPA_unshifted, 2)
PP2_unshifted$coef  ## a matrix
}
