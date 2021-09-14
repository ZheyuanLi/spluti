\name{BSpl2PP}

\alias{BSpl2PP}

\title{
Transform a B-spline or many of its replicates to a "PP" or "PPA" object
}

\description{
Construct B-spline basis, then with user-provided basis coefficients, transform the corresponding B-spline(s) to piecewise polynomial(s).
}

\usage{
BSpl2PP(x, fun = "bs", df = NULL, degree = 3, intercept = FALSE,
        knots = NULL, Boundary.knots = range(x), coef)
}

\arguments{
  \item{fun}{Character name of the function used to construct B-spline basis. Currently only "bs" and "ns" from package \pkg{splines} are supported. However, it is planned that additional functions "pbs" and "pns" be added in \pkg{spluti} to provide basis construction of periodic B-spline and periodic natural cubic spline.}
  \item{x}{The predictor variable. Missing values are allowed but omitted. This argument is passed to \code{fun}.}
  \item{df}{Desired degree of freedom (number of B-spline basis) of the spline. This argument is passed to \code{fun}.}
  \item{degree}{Polynomial degree of the spline. This argument is passed to \code{fun}. Note that the degree and the order of a B-spline are not the same thing. The order is (degree + 1).}
  \item{intercept}{Default is FALSE, implying that the spline is pinned to 0 at its left boundary knot. If TRUE, an intercept is included in the spline so that the spline can move away from 0 at this boundary. This argument is passed to \code{fun}.}
  \item{knots}{Interior knots for the spline. This argument is passed to \code{fun}. If provided, they must lie within range of non-NA x-values. If not provided, they will be automatically placed at quantiles of non-NA x-values by \code{fun}.}
  \item{Boundary.knots}{Clamped boundary knots for the spline, default to the range of the non-NA x-values. If set to other values, they must lie beyond range of non-NA x-values. This argument is passed to \code{fun}.}
  \item{coef}{Coefficients for the constructed B-spline basis. If a vector, a single spline is produced as a "PP" object; if a matrix (where each column gives a replicate of basis coefficients), multiple replicates of a spline are produced as a "PPA" object.}
}

\details{
\subsection{How degree of freedom is calculated in a B-spline}{

If a B-spline has \eqn{k} interior knots, then counting in the two boundary knots there are \eqn{(k + 2)} knots and hence \eqn{(k + 1)} polynomial pieces. If the polynomial is of degree \eqn{d}, then each piece has \eqn{(d + 1)} coefficients, resulting in a total of \eqn{(k + 1)(d + 1)} coefficients. However, the \eqn{0} to \eqn{(d - 1)}-th derivatives of a spline is continuous at interior knots, thus polynomial coefficients are subject to \eqn{kd} constraints, leaving \eqn{(k + 1)(d + 1) - kd = (k + d + 1)} effective number of coefficients, known as the \emph{degree of freedom}. It is also identical to the number of B-spline basis for this spline. Removing intercept from the spline further decreases degree of freedom by 1. In addition, restricting the B-spline to be periodic reduces degree of freedom by an extra 1. Generally, there is \deqn{df = (k + d + 1) - (1 - \textrm{intercept}) - \textrm{periodic} = k + d + \textrm{intercept} - \textrm{periodic},} where intercept and periodic are TRUE/FALSE or 1/0.

Natural cubic spline has \eqn{d = 3}, and there are two extra natural boundary conditions, thus there is \deqn{df = k + 3 + \textrm{intercept} - \textrm{periodic} - 2 = k + 1 + \textrm{intercept} - \textrm{periodic},} where intercept and periodic are TRUE/FALSE or 1/0.
}

\subsection{Two usual ways to specify B-spline construction}{

To construct B-spline basis with \code{fun}, we can either specify \code{df} or \code{knots}. With either one set, the other can be implied. We rarely given both at the same time, to avoid inconsistency due to misspecification.

Note that when only specifying \code{df}, there is a minimum degree of freedom at \eqn{k = 0}. As a result, it must be satisfied that \eqn{df \ge d + \textrm{intercept} - \textrm{periodic}} for B-spline and \eqn{df \ge 1 + \textrm{intercept} - \textrm{periodic}} for natural cubic spline.
}

\subsection{Difference between \emph{BSpl2PP} and \emph{RegBSpl2PP}}{

See \code{\link{RegBSpl2PP}} for details.
}
}

\value{
A "PP" or "PPA" object.
}

\note{
\code{BSpl2PP} forwards to \code{fun} all arguments it needs for B-spline basis construction, but \code{BSpl2PP} does not allow users to specify x-values beyond the range of boundary knots, which is valid in using \code{fun} itself.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

######################
## 1. a "PP" object ##
######################

## construct a random natural cubic spline and represent it as a "PP" object
spl <- BSpl2PP(x = rnorm(10), fun = "ns", df = 7, intercept = TRUE, coef = runif(7))
spl  ## print the "PP" object

#######################
## 2. a "PPA" object ##
#######################

## construct 10 random natural cubic splines and represent them as a "PPA" object
## note that all splines are pinned at 0 at the left boundary, when intercept = FALSE
spls1 <- BSpl2PP(x = rnorm(10), fun = "ns", df = 7, intercept = FALSE,
                 coef = matrix(runif(70), 7))
spls1  ## print the "PPA" object
plot(spls1)  ## plot the "PPA" object

## by contrast when intercept = TRUE, all splines have a "free" left boundary value
spls2 <- BSpl2PP(x = rnorm(10), fun = "ns", df = 7, intercept = TRUE,
                 coef = matrix(runif(70), 7))
spls2  ## print the "PPA" object
plot(spls2)  ## plot the "PPA" object
}
