\name{Poly2PP}

\alias{Poly2PP}

\title{
Convert a polynomial or many of its replicates to a "PP" or "PPA" object
}

\description{% motivation
Sometimes we want to add a polynomial \eqn{p(x)} defined on \eqn{(-\infty, \infty)} to a spline. Since a spline is always handled as a piecewise polynomial in shifted form by \pkg{spluti}, the addition is straightforward once \eqn{p(x)} is expressed in the same piecewise fashion. Taking coefficients of \eqn{p(x)} and knots of the spline, \code{Poly2PP} can convert \eqn{p(x)} to such required form. The resulting "PP" or "PPA" object is then ready for addition with another "PP" or "PPA" object using \code{\link{AddPP}}.
}

\usage{
Poly2PP(pc, knots)
}

\arguments{
\item{pc}{Polynomial coefficients in ascending order of power. It can either be a vector if there is only a single polynomial to transform, or a matrix if there are multiple replicates of a polynomial to transform. In the latter case, each column of the matrix should give coefficients of one replicate.}
\item{knots}{Knots of the resulting piecewise polynomial.}
}

\value{
A "PP" or "PPA" object.
}

\note{
In regression models, polynomial coefficients may be estimated by a polynomial term constructed from \code{stats::poly}. Note that this function by default constructs an orthogonal polynomial rather than an ordinary (or raw) polynomial, and the corresponding regression coefficients can not be interpreted as polynomial coefficients.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
## list of examples:
## 1. express a single polynomial in piecewise fashion
## 2. express multiple replicates of a polynomial in piecewise fashion
## 3. a fixed-effect spline and + random linear polynomials
## 4. a fixed-effect spline and + random quadratic polynomials

require(spluti)

#########################################################
## 1. express a single polynomial in piecewise fashion ##
#########################################################

## a linear polynomial "-1 + x"
PP <- Poly2PP(pc = c(-1, 1), knots = (-4):4)

## print the object
PP

## display new piecewise polynomial coefficients
PP$coef

## by express the piecewise polynomial in unshift form
## we will see that all columns are identical to 'pc'
UnshiftPP(PP)$coef

####################################
## 2. convert a single polynomial ##
####################################

## 2 replicates of a linear polynomial: "-1 + x" and "2 - 3x"
PPA <- Poly2PP(pc = cbind(c(-1, 1), c(2, -3)), knots = (-4):4)

## print the object
PPA

## display new piecewise polynomial coefficients
PPA$coef

##############################################################
## 3. a fixed-effect spline and + random linear polynomials ##
##############################################################

require(nlme)
model <- lme(y ~ ns(age, df = 7), random = ~ age | id, data = BMC)

## express the fixed-effect spline as a "PP" object
fx_spl <- RegBSpl2PP(model, "ns(age, df = 7)")
fx_intercept <- model$coefficients$fixed[1]
fx_spl <- AddInterceptPP(fx_spl, fx_intercept)

## express random linear polynomials as a "PPA" object
re_coef <- model$coefficients$random$id
re_poly <- Poly2PP(t(re_coef), fx_spl$knots)

## add "PP" and "PPA" to get estimated spline for each subject
spl <- AddPP(fx_spl, re_poly)

## plot those splines and overlay observations
plot(spl, ylim.hint = range(BMC$y))
with(BMC, points(age, y, col = 8))

#################################################################
## 4. a fixed-effect spline and + random quadratic polynomials ##
#################################################################

}
