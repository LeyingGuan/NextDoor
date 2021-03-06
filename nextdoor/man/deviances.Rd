% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deviance.R
\name{deviances}
\alias{deviances}
\title{Deviance function using the link function value for different families}
\usage{
deviances(ypred, y, family, model0)
}
\arguments{
\item{ypred}{predicted link function values}

\item{y}{observations}

\item{family}{model family}

\item{model0}{trained model}
}
\description{
Gaussian : prediction, binomial/multinomial:probability, poisson: expected value
}
