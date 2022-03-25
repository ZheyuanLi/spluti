\name{AddPP}

\alias{AddPP}

\title{
Add two objects that inherit "PP" class
}

\description{
If \eqn{f(x)} is a piecewise polynomial with knots set \eqn{A} and \eqn{g(x)} is a piecewise polynomial with knots set \eqn{B}, their sum \eqn{h(x) = f(x) + g(x)} is still a piecewise polynomial but has knots set \eqn{C = A \cup B}. \code{AddPP} can add two piecewise polynomial objects, with knots automatically relocated. The function is versatile, able to add two "PP" objects, a "PP" object with a "PPA" object, and two "PPA" objects (with the same number of replicates).
}

\usage{
AddPP(a, b)
}

\arguments{
\item{a, b}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
}

\details{
% When \eqn{A = B}, calculation of \eqn{h(x) = f(x) + g(x)} is straightforward, by doing a piece by piece calculation \eqn{h_i(x) = f_i(x) + g_i(x)}. When \eqn{A \neq B}, we can first transform \eqn{f(x)} and \eqn{g(x)} to \eqn{F(x)} and \eqn{G(x)} that are both defined on knots set \eqn{C}, then calculate \eqn{h(x) = F(x) + G(x)}.

% MapPP
% Consider the transformation from \eqn{f(x)} to \eqn{F(x)}. Since \eqn{A \subset C} and \eqn{C = A + C \textbackslash A}, so \eqn{C} is obtained by inserting extra knots to \eqn{A}. For example, if it turns out that an extra knot \eqn{c} is inserted between \eqn{[x_i, x_{i + 1}]}, then we need to transform the \eqn{f_i(x)} to a piecewise polynomial defined on knots \eqn{x_i, c, x_{i + 1}}. In summary, the transformation from \eqn{f(x)} to \eqn{F(x)} replaces a polynomial piece to more polynomial pieces, whenever additional knots are inserted within that piece.

When adding a "PP" object and another, the result is still a "PP" object.

When adding a "PP" object and a "PPA" object, the piecewise polynomial in the "PP" object is added to each replicate of the "PPA" object, resulting in a "PPA" object.

When adding two "PPA" objects, replicates from two objects are added in one-two-one matching and the result is still a "PPA" object. It is required that the two objects have the same number of replicates.
}

\value{
A "PP" or "PPA" object.
}

\note{
If \eqn{f(x)} is a spline of polynomial degree \eqn{d_1} and \eqn{g(x)} a spline of degree \eqn{d_2}, their sum \eqn{h(x)} is not a spline unless \eqn{d_1 = d_2}. To see this, let's say \eqn{d_1 > d_2}, then \eqn{f^{(d_2)}(x)} is continuous but \eqn{g^{(d_2)}(x)} is not, hence \eqn{h^{(d_2)}(x) = f^{(d_2)}(x) + g^{(d_2)}(x)} is not continuous. See \bold{Examples}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
## list of examples:
## 1. "PP" + "PP", where a spline + a spline is still a spline
## 2. "PP" + "PP", where a spline + a spline is NOT a spline
## 3. "PP" + "PPA"
## 4. "PPA" + "PPA"

require(spluti)

####################
## 1. "PP" + "PP" ##
####################

## a random cubic B-spline with 5 degree of freedom
spl1 <- BSpl2PP(1:10, df = 5, intercept = TRUE, coef = rnorm(5))

## another random cubic B-spline with 8 degree of freedom
spl2 <- BSpl2PP(1:10, df = 8, intercept = TRUE, coef = rnorm(8))

## they have different knots
spl1$knots
spl2$knots

## but they can be added with knots relocated
spl <- AddPP(spl1, spl2)
spl$knots

## the result is still a spline because both splines have the same polynomial degree
op <- par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(spl)  ## piecewise cubic
plot(spl, deriv = 1)  ## piecewise quadratic
plot(spl, deriv = 2)  ## piecewise linear
plot(spl, deriv = 3)  ## piecewise constant
par(op)

############################
## 2. another "PP" + "PP" ##
############################

## a random cubic B-spline with 5 degree of freedom
spl1 <- BSpl2PP(1:10, df = 5, intercept = TRUE, coef = rnorm(5))

## a random quadratic B-spline with 8 degree of freedom
spl2 <- BSpl2PP(1:10, df = 8, degree = 2, intercept = TRUE, coef = rnorm(8))

## they can be added
PP <- AddPP(spl1, spl2)

## but the resulting piecewise polynomial is not a spline
## because its 2nd derivative is not continuous
op <- par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(PP)  ## piecewise cubic
plot(PP, deriv = 1)  ## piecewise quadratic
plot(PP, deriv = 2)  ## not continuous!!
plot(PP, deriv = 3)  ## piecewise constant
par(op)

#####################
## 3. "PP" + "PPA" ##
#####################

## a random cubic B-spline with 5 degree of freedom
PP <- BSpl2PP(1:10, df = 5, intercept = TRUE, coef = rnorm(5))

## 4 replicates of a random cubic B-spline with 8 degree of freedom
PPA <- BSpl2PP(1:10, df = 8, intercept = TRUE, coef = matrix(rnorm(8 * 4), nrow = 8))

## the sum is a "PPA" object, and the order in which objects are given does not matter
PPA1 <- AddPP(PP, PPA)
PPA2 <- AddPP(PPA, PP)
identical(PPA1, PPA2)  ## TRUE

######################
## 4. "PPA" + "PPA" ##
######################

## 4 replicates of a random cubic B-spline with 5 degree of freedom
PPA1 <- BSpl2PP(1:10, df = 5, intercept = TRUE, coef = matrix(rnorm(5 * 4), nrow = 5))

## 4 replicates of a random cubic B-spline with 8 degree of freedom
PPA2 <- BSpl2PP(1:10, df = 8, intercept = TRUE, coef = matrix(rnorm(8 * 4), nrow = 8))

## the sum is still a "PPA" object
AddPP(PPA1, PPA2)
}
