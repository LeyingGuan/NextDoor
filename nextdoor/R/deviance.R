#' Deviance function using the link function value for different families
#'
#' @description Gaussian : prediction, binomial/multinomial:probability, poisson: expected value
#' @param ypred predicted link function values
#' @param y observations
#' @param family model family
#' @param model0 trained model
deviances <- function(ypred, y, family, model0){
    if(family == "gaussian"){
        errors = (ypred-y)^2
    }else if (family == "binomial"){
        y0 =model0$glmnet.fit$classnames
        y = ifelse(y == y0[1], 0, 1)
        ypred = ifelse(ypred < 10^(-4), 10^(-4), ypred)
        ypred = ifelse(ypred > 1-10^(-4), 1-10^(-4),  ypred)
        errors = -2*y*log(ypred)-2*(1-y)*log(1-ypred)
    }else if (family == "multinomial"){
        y0 = model0$glmnet.fit$classnames
        YY = matrix(0, ncol = ncol(ypred), nrow(ypred))
        for(i in 1:nrow(YY)){
            YY[i, which(y0 == y[i])] = 1
        }
        ypred = ifelse(ypred < 10^(-4), 10^(-4), ypred)
        ypred = ifelse(ypred > 1-10^(-4), 1-10^(-4),ypred)
        errors = -2*YY*log(ypred)-2*(1-YY)*log(1- ypred)
        errors = apply(errors,1,sum)
    }else if (family == "poisson"){
        errors = 2*(y*log(ypred)-(y-ypred))
    }else{
        stop("model familly not supported!")
    }
    return(errors)
}
