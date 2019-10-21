#' Provide Bootstrap p-values
#'
#' @param errors0 the n by m errors of the original model sequence
#' @param errors a list of n by m errors of the proximal model sequences
#' @param S a vector corresponding of feature indexes considered in the nextdoor models.
#' @param nams a length p vector containing feature names. 
#' @param random_estimates the randomized error estimate of the deterioration proximal - original
#' @param sigmaLeft Left eigenvector matrix for errors of each of the model sequence.
#' @param K number of repetitions estimating the de-biased error
#' @param B number of bootstrap repetitions
#' @param alpha added error with covariance structure*alpha, by default, alpha = .1
#' @param epsilon added error with covariance structure being identity times min(covariance diagonal) times epsilon,
#'        by default, epsilon = 0
#' @param epsilon2 added error in the Bootstrap step being min(covariance diagonal) times epsilon2,
#'        by default, epsilon2 = 0
#' @param Bindex n by B index matrix for bootstrap. if Bindex == NULL,
#'        use the default Bootstrap, otherwise, use the provide matrix
#' @param rescale if rescale == True, perform the mean-rescaled Bootstrap
#' @param selectionType if selectionType == 0, pick the model with the smallest randomized error
#'        if selectionType == 1, use the 1se rule
#' @param one_sds if the selectionType is 1, the we choose the model with smallest index
#'        such that model error(randomized) <= one_sds[i] + min error(randomized)
#' @param trace if trace == True, print the p-value process
#' @importFrom MASS mvrnorm
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib nextdoor
#' @return pvalues: p values for nextdoor analysis
nextdoor_unconditional_test <- function(errors0, errors, S, nams = NULL, random_estimates,
sigmaLeft, epsilon = NULL, epsilon2 = NULL,alpha = 0.1,
K = 100, B = 1000, Bindex = NULL, rescale = T,
selectionType = 0, one_sds = rep(0, ncol(errors0)),  trace=T){
    n = nrow(errors0)
    s = length(S)
    m = ncol(errors0)
    if(rescale){
        #the rescaling steps simply rescale the variance of the underlying mean to match the truth
        mean0 = apply(errors0, 2, mean)
        cov0 = cov(errors0, errors0)/n
        scaler = diag(cov0)/(m-1) - sum(cov0)/(m*(m-1))
        scaler = ifelse(scaler > 0, scaler, 0)
        scaler = (var(mean0) - scaler)/var(mean0)
        scaler = ifelse(scaler > 0, scaler, 0)
        mean1 = (mean0 - mean(mean0))*scaler
        errors0 = t(apply(errors0, 1, function(x) (x - mean0+mean1)))
    }
    if(is.null(Bindex)){
        Bindex = fastBoostrap(n = n, B = B)
    }else{
        B = ncol(Bindex)
        Bindex = Bindex - 1
    }
    try(if(sum(is.na(errors0)) > 0) stop("NaN in error matrix"))
    #bootstrap of Q_j - Q - (Ej - E) r - z > 0 -> Ej - E >0
    pvalues = c()
    for(j in S){
        if(trace){print(paste0("###nextdoor analysis to feature ", nams[j], "#####################"))}
        errors1 = errors[[j]]
        try(if(sum(is.na(errors1)) > 0) stop("NaN in error matrix"))
        try(if(nrow(errors1)!=nrow(errors0)) stop("wrong dimension"))
        try(if(ncol(errors1)!=ncol(errors0)) stop("wrong dimension"))
        errors0Mean  = apply(errors0,2,mean)
        errorsMean  = apply(errors1,2,mean)
        z = boostrap_proximity_random(errors0 = errors0, errors = errors1,
        errors0Mean = errors0Mean, errorsMean = errorsMean,
        Bindex = Bindex, epsilon = sqrt(epsilon), epsilon2 = sqrt(epsilon2),
        alpha = alpha, sigmaLeft = sigmaLeft[[j]], K = K,
        selectionType = selectionType, sds = one_sds)
        pvalues = c(pvalues, (sum(random_estimates[[j]] - z < 0)+.1)/(B+.1))
    }
    return(pvalues)
}
