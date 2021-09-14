\name{ICSpl2PP}

\alias{ICSpl2PP}

\title{
Construct and transform an interpolation cubic spline to a "PP" object
}

\description{% motivation
\code{stats::spline} and \code{stats::splinefun} can interpolate (x, y) data using various types of cubic splines, and evaluate the spline or its derivatives at new x-coordinate values. However, there is no built-in method to solve an inverse problem. \code{ICSpl2PP} restructures an interpolation spline as a "PP" object, so that an inverse problem can be solved using the "PP" method for \code{solve}.
}

\usage{
ICSpl2PP(x, y, method)
}

\arguments{
  \item{x, y}{Vectors giving the coordinates of the points to be interpolated.}
  \item{method}{Methods for cubic spline interpolation, namely "fmm", "natural", "periodic", "monoH.FC" and "hyman" as used by \code{stats::spline} and \code{stats::splinefun}.}
}

\value{
A "PP" object.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
## list of interpolation cubic spline examples
## 1. a non-monotonic, non-periodic spline with "fmm" method
## 2. a non-monotonic, non-periodic spline with "natural" method
## 3. a periodic spline with "periodic" method
## 4. a monotonic spline with "monoH.FC" method
## 5. a monotonic spline with "hyman" method

require(spluti)

###############################################################
## 1. a non-monotonic, non-periodic spline with "fmm" method ##
###############################################################

## perform interpolation and transform the spline to a "PP" object
fmm <- ICSpl2PP(toyinterp$x, toyinterp$y1, "fmm")

## print this "PP" object
fmm

## print polynomial equations as formatted strings for all polynomial pieces
summary(fmm)

## plot the spline and its 1st ~ 3rd derivatives
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
plot(fmm, show.knots = TRUE)
points(toyinterp$x, toyinterp$y1, pch = 19)
title("interpolation spline")
plot(fmm, show.knots = TRUE, deriv = 1)
title("1st derivative")
plot(fmm, show.knots = TRUE, deriv = 2)
title("2nd derivative")
plot(fmm, show.knots = TRUE, deriv = 3)
title("3rd derivative")
par(op)

## backsolve the spline given y = 2.85
xr_fmm <- solve(fmm, b = 2.85)
plot(fmm)
abline(h = 2.85, lty = 2)
points(xr_fmm, rep.int(2.85, length(xr_fmm)), pch = 19)

## find all extrema
xs_fmm <- solve(fmm, deriv = 1)
ys_fmm <- predict(fmm, xs_fmm)
plot(fmm)
points(xs_fmm, ys_fmm, pch = 19)
title("local extrema")

###################################################################
## 2. a non-monotonic, non-periodic spline with "natural" method ##
###################################################################

## perform interpolation and transform the spline to a "PP" object
natural <- ICSpl2PP(toyinterp$x, toyinterp$y1, "natural")

## print this "PP" object
natural

## print polynomial equations as formatted strings for all polynomial pieces
summary(natural)

## plot the spline and its 1st ~ 3rd derivatives
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
plot(natural, show.knots = TRUE)
points(toyinterp$x, toyinterp$y1, pch = 19)
title("interpolation spline")
plot(natural, show.knots = TRUE, deriv = 1)
title("1st derivative")
plot(natural, show.knots = TRUE, deriv = 2)
title("2nd derivative")
plot(natural, show.knots = TRUE, deriv = 3)
title("3rd derivative")
par(op)

## backsolve the spline given y = 2.85
xr_natural <- solve(natural, b = 2.85)
plot(natural)
abline(h = 2.85, lty = 2)
points(xr_natural, rep.int(2.85, length(xr_natural)), pch = 19)

## find all extrema
xs_natural <- solve(natural, deriv = 1)
ys_natural <- predict(natural, xs_natural)
plot(natural)
points(xs_natural, ys_natural, pch = 19)
title("local extrema")

#################################################
## 3. a periodic spline with "periodic" method ##
#################################################

## perform interpolation and transform the spline to a "PP" object
periodic <- ICSpl2PP(toyinterp$x, toyinterp$y2, "periodic")

## print this "PP" object
periodic

## print polynomial equations as formatted strings for all polynomial pieces
summary(periodic)

## plot the spline and its 1st ~ 3rd derivatives
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
plot(periodic, show.knots = TRUE)
points(toyinterp$x, toyinterp$y2, pch = 19)
title("interpolation spline")
plot(periodic, show.knots = TRUE, deriv = 1)
title("1st derivative")
plot(periodic, show.knots = TRUE, deriv = 2)
title("2nd derivative")
plot(periodic, show.knots = TRUE, deriv = 3)
title("3rd derivative")
par(op)

## backsolve the spline given y = 2.85
xr_periodic <- solve(periodic, b = 2.85)
plot(periodic)
abline(h = 2.85, lty = 2)
points(xr_periodic, rep.int(2.85, length(xr_periodic)), pch = 19)

## find all extrema
xs_periodic <- solve(periodic, deriv = 1)
ys_periodic <- predict(periodic, xs_periodic)
plot(periodic)
points(xs_periodic, ys_periodic, pch = 19)
title("local extrema")

\dontrun{
##################################################
## 4. a monotonic spline with "monoH.FC" method ##
##################################################

## perform interpolation and transform the spline to a "PP" object
monoH.FC <- ICSpl2PP(toyinterp$x, toyinterp$y3, "monoH.FC")

## print this "PP" object
monoH.FC

## print polynomial equations as formatted strings for all polynomial pieces
summary(monoH.FC)

## plot the spline and its 1st ~ 3rd derivatives
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
plot(monoH.FC, show.knots = TRUE)
points(toyinterp$x, toyinterp$y3, pch = 19)
title("interpolation spline")
plot(monoH.FC, show.knots = TRUE, deriv = 1)
title("1st derivative")
plot(monoH.FC, show.knots = TRUE, deriv = 2)
title("2nd derivative")
plot(monoH.FC, show.knots = TRUE, deriv = 3)
title("3rd derivative")
par(op)

## backsolve the spline given y = 2.85
xr_monoH.FC <- solve(monoH.FC, b = 2.85)
plot(monoH.FC)
abline(h = 2.85, lty = 2)
points(xr_monoH.FC, rep.int(2.85, length(xr_monoH.FC)), pch = 19)
}
###############################################
## 5. a monotonic spline with "hyman" method ##
###############################################

## perform interpolation and transform the spline to a "PP" object
hyman <- ICSpl2PP(toyinterp$x, toyinterp$y3, "hyman")

## print this "PP" object
hyman

## print polynomial equations as formatted strings for all polynomial pieces
summary(hyman)

## plot the spline and its 1st ~ 3rd derivatives
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
plot(hyman, show.knots = TRUE)
points(toyinterp$x, toyinterp$y3, pch = 19)
title("interpolation spline")
plot(hyman, show.knots = TRUE, deriv = 1)
title("1st derivative")
plot(hyman, show.knots = TRUE, deriv = 2)
title("2nd derivative")
plot(hyman, show.knots = TRUE, deriv = 3)
title("3rd derivative")
par(op)

## backsolve the spline given y = 2.85
xr_hyman <- solve(hyman, b = 2.85)
plot(hyman)
abline(h = 2.85, lty = 2)
points(xr_hyman, rep.int(2.85, length(xr_hyman)), pch = 19)
}
