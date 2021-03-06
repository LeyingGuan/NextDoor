% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_model.R
\name{train_model}
\alias{train_model}
\title{Perform training}
\usage{
train_model(x, y, family = c("gaussian", "binomial", "poisson",
  "multinomial"), foldid = NULL, lambda = NULL, lambda_extra = NULL,
  lossfun = NULL, standardize = FALSE)
}
\arguments{
\item{x}{n by p training feature matrix}

\item{y}{length Y response vector}

\item{family}{model family}

\item{foldid}{fold id for across-validation, if it is NULL, the samples are randomly divided into nfolds folds.}

\item{lambda}{A user supplied lambda sequence. WARNING: use with care. Avoid supplying a single value for lambda.}

\item{lambda_extra}{pad lambda with a small value}

\item{lossfun}{The user provided loss function used cross-validation. It has the form lossfun(Y_pred, Y),
with Y_pred being the predicted values using the link function. If lossfun == NULL, by default, we use the deviance}

\item{standardize}{whether to standardize the data, standardize = T by default}
}
\value{
model: model sequence trained using cv.glmnet.

errors0: cross-validation error curves.
}
\description{
Perform training
}
