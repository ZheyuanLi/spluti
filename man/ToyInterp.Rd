\name{ToyInterp}
\docType{data}

\alias{toyinterp}

\title{
A toy dataset for interpolation
}

\description{
A toy dataset for making reproducible examples on how to use \pkg{spluti} for interpolation splines.
}

\usage{
toyinterp
}

\format{
A data frame with 4 variables, where
\itemize{
\item{\code{x} contains 10 unevenly spaced data in ascending order of value;}
\item{\code{y1} contains 10 random numbers from N(3, 1) distribution;}
\item{\code{y2} modifies \code{y1} by replacing its last value by its first value;}
\item{\code{y3} sorts \code{y1} into ascending order of value.}
}
(x, y1), (x, y2) and (x, y3) data are intended for general, periodic and monotonic interpolation splines, respectively.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}
