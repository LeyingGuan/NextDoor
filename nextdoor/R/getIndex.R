#' Perform the model index selection
#' @description This function can be used to evaluate selection frequency for general training procedure with errors0 be the error matrix using original procedure with the Bootstrap samples with the same sequence of penalties
#' @param errors0 the n by m errors of the original model sequence
#' @param cv_glm Return from the function cv.glmnet keep = T.
#' @param alpha added error with covariance structure*alpha, by default, alpha = .1
#' @param epsilon added error with covariance structure being identity times min(covariance diagonal) times epsilon, by default, epsilon = 0.05^2
#' @param selectionType if selectionType == 0, pick the model with the smallest randomized error
#' @param one_sds if the selectionType is 1, the we choose the model with smallest index
#'        such that model error(randomized) <= one_sds[i] + min error(randomized)
#' @return model_index: de-biased estimate of the model error for the original procedure
#' @export
getIndex <-function(errors0,cv_glm,alpha = 0.1, epsilon = 0.05^2, selectionType = 0,
one_sds = rep(0, ncol(errors0))){
    m = ncol(errors0);n=nrow(errors0)
    epsilon1 = min(apply(errors0, 2, var))*epsilon
    sigmaHat = cov(errors0, errors0)/n
    diag(sigmaHat) = diag(sigmaHat) +epsilon1
    ADD1 = rnorm(m, sd =sqrt(epsilon1))
    ADD2 =mvrnorm(n =1,mu = rep(0, m), Sigma = sigmaHat)
    errors0_mean = apply(errors0, 2, mean)
    cv1 = ADD1+sqrt(alpha)*ADD2+errors0_mean
    index0 = which.min(cv1)
    one_sds = cv_glm$cvsd
    if(selectionType == 1){
        index0 = which(cv1 <= cv1[index0]+one_sds[index0])[1]
    } else if (selectionType != 0){
        stop("invalid criterion type")
    }
    return(model_index = index0)

}

