\name{ToyReg}
\docType{data}

\alias{ToyReg}
\alias{wtoyreg}
\alias{ltoyreg}

\title{
Two toy datasets for regression
}

\description{
Two toy datasets for making reproducible examples on how to use \pkg{spluti} for regression splines.
}

\usage{
wtoyreg
ltoyreg
}

\format{
\code{wtoyreg} (wide-form toy dataset) is a data frame of 101 observations, with:
\itemize{
\item{1 independent variable \code{x}, taking evenly spaced values on [0, 10];}
\item{4 response variables \code{y1} to \code{y4}, observed with i.i.d. Gaussian noise from 4 random replicates of a cubic B-spline on [0, 10] with 6 degree of freedom (including intercept), respectively, at above \code{x} values.}
}
This dataset is useful for demonstrating single-response or multi-response linear regression model.

\code{ltoyreg} (long-form toy dataset) is adapted from \code{wtoyreg}, by independently and randomly selecting 50 observations from each of (\code{x}, \code{y1}) to (\code{x}, \code{y4}) pair, and stacking those subsamples to long format. It is a data frame of 200 observations, with three variables:
\itemize{
\item{\code{x}, concatenating subsampled \code{x} values for each replicate;}
\item{\code{y}, concatenating subsampled \code{y1}, \code{y2}, \code{y3} and \code{y4} values;}
\item{\code{id}, a 4-level factor variable (with levels "1", "2", "3", "4") indicating which replicate the (\code{x}, \code{y}) pair comes from.}
}
This dataset is useful for demonstrating single-response linear regression model fitting a spline per factor level.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\note{
There are not many factor levels in \code{ltoyreg} to fit a reasonable linear mixed model with a random effect per factor level. For demonstration on linear mixed models, consider the \code{BMC} longitudinal dataset.
}

\examples{
require(spluti)

## scatter plot of 'wtoyreg'
with(wtoyreg, matplot(x, cbind(y1, y2, y3, y4), type = "p", pch = 1, ylab = "y"))

## scatter plot of 'ltoyreg'
with(ltoyreg, plot(x, y, pch = 1, col = as.integer(id)))
}
