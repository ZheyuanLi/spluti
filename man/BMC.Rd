\name{BMC}
\docType{data}

\alias{BMC}

\title{
Bone Mineral Content dataset for demonstration purpose
}

\description{
A small longitudinal dataset containing 20 subjects' Bone Mineral Content (BMC) measured at growing ages.
}

\usage{
%data(BMC, package = "spluti")  % LazyData: false
BMC  % LazyData: true
}

\format{
A dataframe of 271 observations with three variables:
\itemize{
\item{\code{age}, the age of the subject when the measurement is taken;}
\item{\code{y}, the measured BMC;}
\item{\code{id}, a 20-level factor (with levels "1", "2", \ldots, "20") for subjects.}
}
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(spluti)

## colour for each subject
col <- 5 - as.integer(BMC$id) \%\% 5

## pch for each subject
pch <- ceiling(as.integer(BMC$id) / 5)

## differnet subject has different colour AND pch
with(BMC, plot(age, y, col = col, pch = pch))
}
