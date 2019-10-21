#' Perform the model selection, unbiased error estimation and the nextdoor test(p-value)
#'
#' @param errors0 the n by m errors of the original model sequence.
#' @param errors a list of n by m errors of the proximal model sequences.
#' @param S a vector corresponding of indexes of the proximal sequences.
#' @param nams a length p vector containing feature names. By default, nams = NULL, the features are named according to their column indexes in x.
#' @param K number of repetitions estimating the de-biased error
#' @param B number of bootstrap repetitions
#' @param alpha added error with covariance structure*alpha, by default, alpha = .1
#' @param epsilon added error with covariance structure being identity times min(covariance diagonal) times epsilon,
#'        by default, epsilon = 0.05^2
#' @param epsilon2 added error in the Bootstrap step being min(covariance diagonal) times epsilon2,
#'        by default, epsilon2 = 0.05^2
#' @param Bindex n by B index matrix for bootstrap. if Bindex == NULL,
#'        use the default Bootstrap, otherwise, use the provide matrix.
#' @param selectionType if selectionType == 0, pick the model with the smallest randomized error
#'        if selectionType == 1, use the 1se rule
#' @param pv if pv == True, estimate the p-values
#' @param rescale if rescale == True, perform the mean-rescaled Bootstrap
#' @param one_sds if the selectionType is 1, the we choose the model with smallest index
#'        such that model error(randomized) <= one_sds[i] + min error(randomized)
#' @param trace if trace == True, print the p-value process
#' @return debiased_errors0: de-biased estimate of the model error for the original procedure
#' @return debiased_errors: de-biased estimate of the model error for the process excluding a specific feature
#' @return worsen: estimated increase in prediction error
#' @return pv_value: p values for proximity analysis
#' @importFrom MASS mvrnorm
#' @examples
#' data(prostateCancerData);library(ranger); library(MASS)
#' data_train = prostateCancerData$train;nams=prostateCancerData$train$names
#' data = data.frame(cbind(data_train$feature, data_train$response)); colnames(data) = c(nams, "response")
#' ###get the predictions using the original models and using the nextdoor models with different hyperparameters.###########
#' ###here, for illustration purpose and simplicity, we used the errors from random forest without cross-validation.###########
#' set.seed(48)
#' alphas = c(1:10)/10; 
#' errors0 = array(NA, dim = c(length(data$response),length(alphas))); errors=list()
#' for(i in 1:length(alphas)){
#'     model0 = ranger(response~., data = data, alpha = alphas[i])
#'     errors0[,i] = (data$response-model0$predictions)^2
#'     for(j in 1:length(nams)){
#'        if(is.null(errors[j][[1]])){errors[[j]] = array(0, dim = c(length(data$response),length(alphas)))}
#'        data1 = data[,-j];model0 = ranger(response~., data = data1, alpha = alphas[i]) 
#'        errors[[j]][,i] = (data$response-model0$predictions)^2
#'     }
#' }
#'res = nextdoor(errors0 = errors0, errors = errors, S =c(1:length(nams)), nams=nams, B = 1000, alpha = 0.1, pv = TRUE, rescale = TRUE, selectionType = 0,trace = TRUE)
#' @export
nextdoor <- function(errors0, errors, S, nams=NULL, K = 100, B = 1000,
alpha = 0.1, epsilon = 0.05^2, epsilon2 = 0.05^2,
Bindex = NULL,pv = TRUE,  rescale = TRUE,
selectionType = 0, one_sds = rep(0, ncol(errors0)),
trace = T){
    sigmaLeft = list()
    sd1 = list()
    errorsMean = list()
    random_estimates = list()
    error0_random_estimates = 0
    error_random_estimates = list()
    errors0Mean  = apply(errors0,2,mean)
    m = ncol(errors0); n = nrow(errors0)
    sigmaHat0 = cov(errors0)/n
    epsilon1 = min(diag(sigmaHat0))*epsilon
    epsilon21 = min(diag(sigmaHat0))*epsilon2
    for(j in S){
        errors1 = errors[[j]]
        sigmaHat = cov(cbind(errors0, errors1))/n
        m1 = ncol(errors1)
        diag(sigmaHat) <-diag(sigmaHat)+epsilon1
        ADD1 = sapply(1:K, function(i) rnorm(m+m1, sd =sqrt(epsilon1)))
        ADD2 = t(mvrnorm(n =K,mu = rep(0, m+m1), Sigma = sigmaHat))
        errorsMean[[j]]  = apply(errors1,2,mean)
        errors0_random_add = sapply(1:K, function(i) errors0Mean+ADD1[1:m,i]+sqrt(alpha)*ADD2[1:m,i])
        errors0_random_minus = sapply(1:K, function(i) errors0Mean+ADD1[1:m,i]-sqrt(1/alpha)*ADD2[1:m,i])
        errors_random_minus = sapply(1:K, function(i) errorsMean[[j]]+ADD1[(m+1):(m+m1),i]-sqrt(1/alpha)*ADD2[(m+1):(m+m1),i])
        if(selectionType == 0){
            idxs = apply(errors0_random_add, 2, which.min)
        } else if (selectionType == 1){
            if(is.null(one_sds)){stop("NO 1se rule!")}
            idxs = apply(errors0_random_add, 2, function(x){
                ii = which.min(x)
                ii = which(x <= x[ii]+one_sds[ii])[1]
                ii
            })
        }else{
            stop("invalid criterion type")
        }
        error_random_estimates[[j]] = mean(sapply(1:K, function(i) errors_random_minus[idxs[i],i]))
        error0_random_estimates_temp =  mean(sapply(1:K, function(i) errors0_random_minus[idxs[i],i]))
        error0_random_estimates =   error0_random_estimates + error0_random_estimates_temp/length(S)
        random_estimates[[j]] = error_random_estimates[[j]] -error0_random_estimates_temp
        sigmaSVD = svd(sigmaHat)
        sigmaLeft[[j]] = sigmaSVD[[2]]%*%diag(sqrt(sigmaSVD[[1]]))
    }
    p_value = NULL
    if(pv){
        p_value = nextdoor_unconditional_test (errors0 = errors0, errors = errors, S = S, nams = nams,
        random_estimates =random_estimates , rescale = rescale,
        epsilon = epsilon1, epsilon2 = epsilon21, sigmaLeft = sigmaLeft,alpha = alpha,
        K = K, B = B, Bindex = Bindex,selectionType = 0,
        one_sds = one_sds, trace=TRUE)
        names(p_value)=nams[S]
    }
    return(list(debiased_errors0 = error0_random_estimates,
    debiased_errors = error_random_estimates,
    worsen = random_estimates,
    p_value = p_value))
}

