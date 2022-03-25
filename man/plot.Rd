\name{plot.PP}

\alias{plot.PP}

\title{
"PP" and "PPA" methods for generic function plot()
}

\description{
Plot a "PP" or "PPA" object, or its derivatives.
}

\usage{
\method{plot}{PP}(x, spread = 3, deriv = 0, show.knots = FALSE, ylim.hint = NULL, \dots)
}

\arguments{
  \item{x}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
  \item{spread}{A graphical parameter determining how dense the evaluation grid is set up for plotting. Note that this is a multiplier factor. For a spline of polynomial degree \eqn{d}, there will be \code{spread * (d + 1)} evenly spaced evaluation points between two adjacent knots.}
  \item{deriv}{Number of derivatives to take. Not allowed to exceed polynomial degree, otherwise the differentiated piecewise polynomial becomes 0.}
  \item{show.knots}{If \code{TRUE}, grey vertical dotted lines will be drawn at knots of the spline.}
  \item{ylim.hint}{A hint on the y-axis range of the plot. This does not enforce the y-axis range; rather, it suggests that y-axis of the plot should at least cover this range. The actual y-axis range will be jointly determined by the evaluated spline values used for plotting. The argument proves useful when users want to reserver enough space in the plotting region so that they can later draw extra data on the plot. See \bold{Examples}.}
  \item{\dots}{Not used by this method. See \bold{Note} below.}
}

\value{
The function invisibly returns data for plotting in a list of two variables:
\itemize{
\item{\code{x}, a vector of x-coordinate values where splines are evaluated;}
\item{\code{y}, evaluated spline values. either a vector if a "PP" object is plotted, or a matrix (a column per replicate) if a "PPA" object is plotted.}
}
}

\note{
While the method contains \code{\dots} in its arguments, it is unused so do not expect to customize the graphical display of the plot, by passing in general graphical parameters like \code{xlab}, \code{ylab}, etc. However, users can take the invisibly returned \code{x}, \code{y} values by the method, and produce a fancy plot with for example \pkg{ggplot2} package.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

## fit a smoothing spline
sm <- smooth.spline(wtoyreg$x, wtoyreg$y3)

## coerce the "smooth.spline" object to a "PP" object
spl <- SmSpl2PP(sm)

## plot the fitted smoothing spline against data
op <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
## some data are out of the plotting region hence not displayed
plot(spl)
points(wtoyreg$x, wtoyreg$y3, col = 8)
## reserve space so that all data can be displayed in the plotting region
plot(spl, ylim.hint = range(wtoyreg$y3))
points(wtoyreg$x, wtoyreg$y3, col = 8)
par(op)
}
