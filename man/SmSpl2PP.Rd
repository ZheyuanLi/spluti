\name{SmSpl2PP}

\alias{SmSpl2PP}

\title{
Reparametrize a fitted smoothing spline to a "PP" object
}

\description{
\code{stats::smooth.spline} can smooth (x, y) data using a natural cubic spline, and evaluate the spline or its derivatives at new x-coordinate values. However, there is no built-in method to solve an inverse problem. \code{SmSpl2PP} transforms a smoothing spline to a "PP" object, so that an inverse problem can be solved using the "PP" method for \code{solve}.
}

\usage{
SmSpl2PP(SmSpl)
}

\arguments{
  \item{SmSpl}{A smoothing spline, typically a fitted model returned by \code{stats::smooth.spline}}.
}

\details{
The function deals with a fitted natural cubic smoothing spline model created by \code{stats::smooth.spline}. It extracts knots from the fitted spline, evaluates the spline at knots, then calls \code{\link{ICSpl2PP}} to construct a natural cubic interpolation spline.
}

\value{
A "PP" object.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

## fit a smoothing spline
sm <- smooth.spline(wtoyreg$x, wtoyreg$y1)

## coerce the "smooth.spline" object to a "PP" object
spl <- SmSpl2PP(sm)

## print the "PP"
spl

## plot the "PP"
plot(spl)

## find all the roots
xr <- solve(spl)
points(xr, rep.int(0, length(xr)))
abline(h = 0, lty = 2)

## find all stationary / saddle points
xs <- solve(spl, deriv = 1)

## predict the "PP" at stationary / saddle points
ys <- predict(spl, xs)
points(xs, ys, pch = 19)
}
