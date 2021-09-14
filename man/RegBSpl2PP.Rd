\name{RegBSpl2PP}

\alias{RegBSpl2PP}
\alias{RegBSpl2PP.lm}
\alias{RegBSpl2PP.mlm}
\alias{RegBSpl2PP.lme}

\title{
Transform a B-spline or many of its replicates in a fitted regression model to a "PP" or "PPA" object
}

\description{
When a regression model contains a B-spline term constructed by \code{splines::bs} or \code{splines::ns}, \code{RegBSpl2PP} may be used to transform this spline to a piecewise polynomial. There are conditions for this to work; see \bold{Details}.
}

\usage{
RegBSpl2PP(RegModel, BSplTerm)

\method{RegBSpl2PP}{lm}(RegModel, BSplTerm)

\method{RegBSpl2PP}{mlm}(RegModel, BSplTerm)

\method{RegBSpl2PP}{lme}(RegModel, BSplTerm)
}

\arguments{
  \item{RegModel}{A fitted regression model with one or more spline term specified by \code{splines::bs} or \code{splines:ns}.}
  \item{BSplTerm}{A character string giving the name of the spline term to be transformed. The name must be correctly provided, otherwise the function throws a "BSplTerm not found" error. However, the error message also lists names of all terms, making it easy for users to do copy-and-paste for a fix.}
}

\details{

\subsection{Conditions for \emph{RegBSpl2PP} to work successfully}{

\code{RegBSpl2PP} is a convenient generic function for users working with fitted regression models, but it only works successfully in two conditions:

\enumerate{
\item{the model class can be recognized by \code{RegBSpl2PP};}
\item{the basis construction parameters for the B-spline term requested by users can be found in the fitted model object.}
}

When these are met, \code{RegBSpl2PP} can extract those construction parameters and automatically find basis coefficients from regression coefficients to proceed to the transformation.

\bold{Condition 1} may be relaxed as more methods are written for \code{RegBSpl2PP}, but there is a bottom line: the fitted model must be produced by a regression routine that has a formula interface where users can specify their model structure using formula(e). For example, \code{stats::lm}, \code{stats::glm}, \code{nlme::lme} have such formula interface hence "lm", "mlm" and "lme" methods have been provided. (There is no need for a "glm" method, as a "glm" object inherits "lm" class and the "lm" method suffices.) "lme4" methods corresponding to \code{lme4::lmer} and \code{lme4::glmer} is also planned. But there is no way to write a method for \code{glmnet::glmnet} for example.

\bold{Condition 2} should generally be satisfied if the regression routine has been properly written, as for \emph{Safe Prediction} the predict call (containing all basis construction parameters) for \code{splines::bs} and \code{splines::ns} terms should be available in the fitted regression model. However, we have found that such a predict call does not exist in a "lme" object returned by \code{nlme::lme}, if a B-spline term is specified in the formula for the random effect.
}

\subsection{What to do when conditions are unmet}{

Consider using \code{\link{BSpl2PP}} instead. It is a more primitive way to transform B-spline(s) to piecewise polynomial(s). To use it in this settings, users need to

\enumerate{
\item{ensure that parameters, particularly \code{x}.} 
\item{basis coefficients are correctly extracted from regresssion coefficients.}
}

}
}

\value{
A "PP" or "PPA" object.
}

\note{
It is incorrect to think values of a spline term as fitted values of the regression model, even if this spline is the only term in the RHS of the model formula. See \code{\link{AddInterceptPP}}.
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

## get more about this example under ?AddInterceptPP

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
