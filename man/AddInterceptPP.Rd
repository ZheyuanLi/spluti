\name{AddInterceptPP}

\alias{AddInterceptPP}

\title{
Add regression intercept to a "PP" or "PPA" object
}

\description{% motivation
A regression B-spline constructed by \code{splines::bs} and \code{splines::ns} does not contain intercept by default, so that the spline is pinned to 0 at its left boundary knot. This is a reasonable setup as otherwise the intercept in the regression model will be unidentifiable. However, for plotting and prediction purpose or to solve an inverse problem, we may want to incorporate the estimated model intercept to the fitted spline, and \code{AddInterceptPP} is a handy function for doing this.
}

\usage{
AddInterceptPP(PPObj, intercept)
}

\arguments{
\item{PPObj}{An object that inherits "PP" class, typically a "PP" or "PPA" object.}
\item{intercept}{Regression intercept, must be provided by users. It can be either a single value or a vector of values, see \bold{Details}.}
}

\details{
Once a regression spline is transformed to a "PP" or "PPA" object, \code{AddInterceptPP} can be used to add intercept to the piecewise polynomial. This function spares users the need to manually update polynomial coefficients, but does not automatically extract intercept value(s) from a fitted regression model. Therefore, users are required to provide intercept value(s) themselves, which usually means that they need to examine carefully the returned coefficients of a model object, and correctly pick out intercept value(s). The extraction process depends on model class and model structure and thus varies a lot. The \bold{Examples} section below offers a good list of examples but is far from complete. Users are encouraged to contribute new examples by contacting the maintainer.

\code{AddInterceptPP} is clever about its input \code{PPObj} and \code{intercept}:
\itemize{
\item{"PP" + single intercept value = "PP";}
\item{"PP" + vector of intercept values = "PPA", i.e., each intercept value is added to the "PP" object in turn;}
\item{"PPA" + single intercept value = "PPA", i.e., a common intercept value is added to all replicates in the "PPA" object;}
\item{"PPA" + vector of intercept values = "PPA", i.e., in one-to-one matching, an intercept is added to a replicate. It is required that there are as many values in \code{intercept} as replicates in \code{PPObj}.}
}
}

\value{
A "PP" or "PPA" object.
}

\note{
In other context than regression, \code{AddInterceptPP} can be used, for example, to add a constant to a "PP" or "PPA" object.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
## list of regression examples
## 1. a single spline, via 'stats::lm'
## 2. 4 replicates of a spline, via 'stats::lm'
## 3. a spline for each factor level, via 'stats::lm'
## 4. a fixed-effect spline + random intercept, via 'nlme::lme'
## 5. a fixed-effect spline + a random-effect spline, via 'nlme::lme'

require(spluti)

#########################################
## 1. a single spline, via 'stats::lm' ##
#########################################

## fit 'y1' as a spline of 'x'
model <- lm(y1 ~ bs(x, df = 5), data = wtoyreg)

## convert the fitted bs() term to a "PP" and plot it
## the fitted bs() term is away from data by an evident vertical shift
spl <- RegBSpl2PP(model, "bs(x, df = 5)")
plot(spl, ylim.hint = range(wtoyreg$y1))
with(wtoyreg, points(x, y1, col = 8))

## after including intercept, the fitted spline now follows data
intercept <- model$coefficients[1]
spl <- AddInterceptPP(spl, intercept)
plot(spl, ylim.hint = range(wtoyreg$y1))
with(wtoyreg, points(x, y1, col = 8))

##################################################
## 2. 4 replicates of a spline, via 'stats::lm' ##
##################################################

## fit 'y1', 'y2', 'y3' and 'y4' as a spline of 'x'
model <- lm(cbind(y1, y2, y3, y4) ~ bs(x, df = 5), data = wtoyreg)

## convert the fitted bs() term to a "PPA"
spl <- RegBSpl2PP(model, "bs(x, df = 5)")

## add model intercept (a vector of 4 values) to "PPA" (4 replicates of a spline)
## now the fitted spline follows data
intercept <- model$coefficients[1, ]
spl <- AddInterceptPP(spl, intercept)
plot(spl, ylim.hint = with(wtoyreg, range(y1, y2, y3, y4)))
with(wtoyreg, matpoints(x, cbind(y1, y2, y3, y4), pch = 1, col = 8))

########################################################
## 3. a spline for each factor level, via 'stats::lm' ##
########################################################

## fit a spline for each factor level
model <- lm(y ~ bs(x, df = 5) + id + bs(x, df = 5):id, data = ltoyreg)

## baseline-level spline
spl_base <- RegBSpl2PP(model, "bs(x, df = 5)")
intercept_base <- model$coefficients[1]
spl_base <- AddInterceptPP(spl_base, intercept_base)

## factor-level spline
spl_fctr_coef <- model$coefficients[10:24]
spl_fctr_coef <- matrix(spl_fctr_coef, ncol = 3)
spl_fctr <- BSpl2PP(ltoyreg$x, fun = "bs", df = 5, coef = spl_fctr_coef)
intercept_fctr <- model$coefficients[7:9]
spl_fctr <- AddInterceptPP(spl_fctr, intercept_fctr)

##################################################################
## 4. a fixed-effect spline + random intercept, via 'nlme::lme' ##
##################################################################

require(nlme)
model <- lme(y ~ ns(age, df = 7), random = ~ 1 | id, data = BMC)

## fixed-effect spline
fx_spl <- RegBSpl2PP(model, "ns(age, df = 7)")  ## a "PP"
fx_intercept <- model$coefficients$fixed[1]  ## single value
fx_spl <- AddInterceptPP(fx_spl, fx_intercept)  ## still a "PP"

## random intercept
re_intercept <- model$coefficients$random$id[, 1]  ## a vector of 20 values

## add two parts together
spl <- AddInterceptPP(fx_spl, re_intercept)  ## now a "PPA"

## plot fitted spline and data
plot(spl, ylim.hint = range(BMC$y))
with(BMC, points(age, y, col = 8))

########################################################################
## 5. a fixed-effect spline + a random-effect spline, via 'nlme::lme' ##
########################################################################

## fitting a linear mixed model
require(nlme)
model <- lme(y ~ ns(age, df = 7), random = ~ ns(age, df = 3) | id, data = BMC)

## fixed-effect spline
fx_spl <- RegBSpl2PP(model, "ns(age, df = 7)")  ## a "PP"
fx_intercept <- model$coefficients$fixed[1]  ## a single value
fx_spl <- AddInterceptPP(fx_spl, fx_intercept)  ## still a "PP"

## random-effect spline
re_spl_coef <- model$coefficients$random$id[, -1, drop = FALSE]
re_spl <- BSpl2PP(BMC$age, fun = "ns", df = 3, coef = t(re_spl_coef))  ## a "PPA"
re_intercept <- model$coefficients$random$id[, 1]  ## a vector of 20 values
re_spl <- AddInterceptPP(re_spl, re_intercept)  ## still a "PPA"

## add two splines together
spl <- AddPP(fx_spl, re_spl)  ## "PP" + "PPA" = "PPA"

## plot fitted spline and data
plot(spl, ylim.hint = range(BMC$y))
with(BMC, points(age, y, col = 8))
}
