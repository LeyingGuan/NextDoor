#' Perform training
#'
#' @param x n by p training feature matrix
#' @param y length Y response vector
#' @param family model family
#' @param foldid fold id for across-validation, if it is NULL, the samples are randomly divided into nfolds folds.
#' @param lambda A user supplied lambda sequence. WARNING: use with care. Avoid supplying a single value for lambda.
#' @param lambda_extra pad lambda with a small value
#' @param standardize whether to standardize the data, standardize = T by default
#' @param lossfun The user provided loss function used cross-validation. It has the form lossfun(Y_pred, Y),
#'        with Y_pred being the predicted values using the link function. If lossfun == NULL, by default, we use the deviance
#' @return model: model sequence trained using cv.glmnet.
#' @return errors0: cross-validation error curves.
train_model<-function(x, y, family=c("gaussian","binomial","poisson","multinomial"), foldid =NULL, lambda = NULL, lambda_extra = NULL, 
                      lossfun = NULL, standardize = FALSE){
    cv_glm = cv.glmnet(x, y, family = family, lambda = lambda_extra, foldid = foldid, standardize = standardize,keep = TRUE)
    foldid = cv_glm$foldid
    ypred = cv_glm$fit.preval[,1:length(lambda)]
    if(is.null(lossfun)){
        errors0 = apply(ypred, 2, deviances, y = y, family = family, model0 = cv_glm)
    }else{
        errors0 = apply(ypred, 2, lossfun, y = y)
    }
    return(list(model = cv_glm, errors0 = errors0))
}

#' Perform training, model selection
#'
#' @param x n by p training feature matrix
#' @param y length Y response vector
#' @param family model family
#' @param foldid fold id for across-validation, if it is NULL, the samples are randomly divided into nfolds folds.
#' @param lambda A user supplied lambda sequence. WARNING: use with care. Avoid supplying a single value for lambda.
#' @param lambda_extra lambda padded with a small value.
#' @param standardize whether to standardize the data, standardize = T by default
#' @param lossfun The user provided loss function used cross-validation. It has the form lossfun(Y_pred, Y),
#'        with Y_pred being the predicted values using the link function. If lossfun == NULL, by default, we use the deviance.
#' @param alpha added error with covariance structure*alpha, by default, alpha = .1
#' @param epsilon added error with covariance structure being identity times min(covariance diagonal) times epsilon, by default, epsilon = 0.05^2
#' @param selectionType if selectionType == 0, pick the model with the smallest randomized error
#'        if selectionType == 1, use the 1se rule
#' @return S: selected features indexes.
#' @return index0: selected model index.
train_model2<-function(x, y, family=c("gaussian","binomial","poisson","multinomial"),foldid =NULL, lambda = NULL, lambda_extra=NULL,
                       lossfun = NULL,standardize = FALSE,epsilon, alpha , selectionType){
  p = ncol(x); n =nrow(x)
  cv_glm = cv.glmnet(x, y, family = family, lambda = lambda_extra, foldid = foldid, standardize = standardize,keep = TRUE)
  ypred = cv_glm$fit.preval[,1:length(lambda)]
  if(is.null(lossfun)){
    errors0 = apply(ypred, 2, deviances, y = y, family = family, model0 = cv_glm)
  }else{
    errors0 = apply(ypred, 2, lossfun, y = y)
  }
  m = ncol(errors0)
  epsilon1 = min(apply(errors0, 2, var))*epsilon
  sigmaHat = cov(errors0, errors0)/n
  diag(sigmaHat) = diag(sigmaHat) +epsilon1
  ADD1 = rnorm(m, sd =sqrt(epsilon1))
  ADD2 =mvrnorm(n =1,mu = rep(0, m), Sigma = sigmaHat)
  errors0_mean = apply(errors0, 2, mean)
  cv1 = ADD1+sqrt(alpha)*ADD2+errors0_mean
  index0 = which.min(cv1)
  one_sds = (cv_glm$cvsd)[1:length(lambda)]
  if(selectionType == 1){
    index0 = which(cv1 <= cv1[index0]+one_sds[index0])[1]
  }
  S = unlist(predict(cv_glm, s = lambda[index0], type = "n"))
  return(list(S = S, index0 = index0))
}





