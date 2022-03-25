\name{spluti}
\docType{package}

\alias{spluti}
\alias{spluti-package}

\title{
An introduction to spluti
}

\description{

\subsection{Piecewise polynomial definition of a spline}{

A spline \eqn{f(x)} of degree \eqn{d} (or order \eqn{(d + 1)}) and \eqn{k} knots \eqn{x_1, x_2, \ldots, x_k}, is by definition a set of \eqn{(k - 1)} piecewise polynomial of degree \eqn{d} that are \eqn{(d - 1)} times continuously differentiable at the knots. The \eqn{i}-th polynomial piece \eqn{f_i(x)} on interval \eqn{[x_i, x_{i + 1}]} is often parametrized in a \emph{shifted} form: \deqn{f_i(x) = a_{0,i} + a_{1,i}(x - x_i) + a_{2,i}(x - x_i) ^ 2 + \ldots + a_{d,i}(x - x_i) ^ d.} This form is appealing because polynomial coefficients are directly associated with derivatives at knots: \deqn{f^{(m)}(x_i) = f_i^{(m)}(x_i) = m! \times a_{mi}.} It also make it easy to integrate the spline.
} % subsection: Piecewise polynomial definition of a spline

\subsection{Representation of a spline in applications}{

Piecewise polynomial representation of a spline is conventionally used in interpolation problems, but not in regression problems and smoothing problems. For example, R routines \code{stats::spline} and \code{stats::splinefun} for interpolation cubic spline readily saves piecewise polynomial information in its returned environment, but routine \code{stats::smooth.spline} for fitting a smoothing spline uses B-spline representation for the resulting natural cubic spline. In addition, splines constructed by \code{splines::bs} and \code{splines::ns}, later passed to regression routines like \code{stats::lm}, \code{stats::glm}, \code{nlme::lme}, \code{lme4::lmer}, \code{lme4::glmer}, also represent splines with B-spline basis.
} % subsection: Representation of a spline in applications

\subsection{Why piecewise polynomial representation is useful}{

Coefficients for B-spline basis are difficult to interpret. It might thus be helpful to reparametrize a fitted regression spline or smoothing spline into piecewise polynomial for appreciation of broder audience.

In fact, piecewise polynomial representation makes it easier to post-process a fitted spline. For example, evaluation of a spline and its derivatives is straightforward. Moreover, an inverse problem: \deqn{f^{(m)}(x) = b} with \eqn{0 \le m \le {d - 1}} can be easily solved, by simply finding real roots of the polynomial \eqn{f^{(m)}(x) - b}, piece by piece, which can be done by for example, R routine \code{base::polyroot}. Many important practical questions can be stated as such inverse problem, for example, finding extrema and extreme values of a spline is identical to the \eqn{m = 0}, \eqn{b = 0} case.
} % subsection: Why piecewise polynomial representation is useful

\subsection{What spluti offers}{

Package \pkg{spluti} offers two sets of routines:

\enumerate{
\item routines that reparametrize splines constructed in applications to piecewise polynomial;
\item routines that process piecewise polynomial objects.
} % enumerate

The first set is best summarized by the table below, where links to relevant routines and abbreviations used in naming are listed.

% an empty first column helps aligning
\tabular{llll}{
\tab PP \tab Piecewise Polynomial Object Class \tab \link{PP} \cr
\tab PPA \tab Piecewise Polynomial Array Object class \tab \link{PPA} \cr
\tab Poly \tab Polynomial defined on \eqn{(-\infty, \infty)} \tab \code{\link{Poly2PP}} \cr
\tab ICSpl \tab Interpolation Cubic Spline \tab \code{\link{ICSpl2PP}} \cr
\tab SmSpl \tab Smoothing Spline \tab \code{\link{SmSpl2PP}} \cr
\tab BSpl \tab B-Spline \tab \code{\link{BSpl2PP}} \cr
\tab RegBSpl \tab Regression B-Spline \tab \code{\link{RegBSpl2PP}}
} % tabular

The second set contains on one hand "PP" method (and "PPA" method only if a different handling from "PP" method is necessary) for some standard generics, and on the other hand some extra utility routines proved to be practically useful. See summary tables below.

% an empty first column helps aligning
\tabular{llll}{
\tab \code{print} \tab print "PP" and "PPA" objects in a pretty way \tab \code{\link{print.PP}}, \code{\link{print.PPA}} \cr
\tab \code{summary} \tab format piecewise polynomial equations as strings \tab \code{\link{summary.PP}}, \code{\link{summary.PPA}} \cr
\tab \code{plot} \tab plot piecewise polynomial or its derivatives \tab \code{\link{plot.PP}} \cr
\tab \code{predict} \tab predict piecewise polynomial or its derivatives \tab \code{\link{predict.PP}} \cr
\tab \code{solve} \tab solve an inverse problem \tab \code{\link{solve.PP}}, \code{\link{solve.PPA}}
} % tabular

% an empty first column helps aligning
\tabular{lll}{
\tab add regression intercept to a regression spline \tab \code{\link{AddInterceptPP}} \cr
\tab differentiate piecewise polynomial \tab \code{\link{DerivPP}} \cr
\tab add two piecewise polynomial \tab \code{\link{AddPP}} \cr
\tab extract one "PP" object from a "PPA" object \tab \code{\link{PPAEntry}} \cr
} % tabular
} % subsection: What spluti offers

%\subsection{How spluti does reparametrization?}{

%The method can differ depending on how spline is constructed. One method that always work is to conducts reparametrization by exact least squares fitting. \eqn{d + 1} evenly spaced points are taken on \eqn{[x_i, x_{i + 1}]}, where spline values are predicted. Given those \eqn{(x, y)} data, a polynomial of degree \eqn{d} can be exactly fitted. The resulting piecewise polynomial coefficients \eqn{a_{0i}, a_{1i}, \ldots, a_{di}} are then stored in the \eqn{i}-th column of a \eqn{(d + 1) \times (k - 1)} coefficient matrix.
%}

%\subsection{Limitation of the piecewise polynomial representation}{

%For regression splines and smoothing splines, it is not easy to obtain confidence interval for either piecewise polynomial coefficients or 
%} % subsection: Limitation of the piecewise polynomial representation
} % description

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}
