#' Perform model selection, unbiased error estimation and the nextdoor test(p-value/model score) for functions in glmnet
#'
#' @param x n by p training feature matrix.
#' @param y Length y response vector.
#' @param nams a length p vector containing feature names. By default, nams = NULL, the features are named according to their column indexes in x.
#' @param cv_glm Return from the function cv.glmnet keep = T.
#' @param family Response type, it can be one of "gaussian","binomial","poisson","multinomial". By default, family = "gaussian". It should be consistent with cv_glm.
#' @param lossfun A user-specific loss function for model evaluation. If loss = NULL, by default, we will use the deviance.  
#' @param standardize Whether to standardize the data, standardize = T by default. It should be consistent with cv_glm.
#' @param K number of repetitions estimating the de-biased error.
#' @param B number of bootstrap repetitions estimating the p-value.
#' @param alpha added errors' level with added errors being covariance structure*alpha. By default, alpha = .1.
#' @param epsilon added errors' level  with added errors being identity*min(covariance diagonal)\*epsilon. By default, epsilon = 0.05^2.
#' @param epsilon2 added errors' level in the Bootstrap with  added errors being min(covariance diagonal)*epsilon2.
#'        By default, epsilon2 = 0.05^2
#' @param Bindex n by B index matrix for bootstrap. If Bindex == NULL,
#'        use the randomly generated bootstrap samples, otherwise, use the provide matrix.
#' @param selectionType if selectionType == 0, pick the model with the smallest randomized error,
#'        if selectionType == 1, use the 1se rule.
#' @param pv if pv == True, estimate the p-values.
#' @param rescale if rescale == True, perform the mean-rescaled Bootstrap.
#' @param score if score == True, provide model scores.
#' @param B1 number of repetitions for the paired Bootstrap used to create model score.
#' @param Bindex1 n by B1 index matrix for paired bootstrap in the model score step. If Bindex1 == NULL,
#'        use the default Bootstrap, otherwise, use the provide matrix.
#' @param trace if trace == True, print the p-values as they are calculated.
#' @return model0: original model sequence.
#' @return models: list of nextdoor model sequences.
#' @return errors0: original error matrix.
#' @return errors: list of nextdoor model matrices.
#' @return debiased_errors0: de-biased estimate of the prediction error for the original process.
#' @return debiased_errors: de-biased estimates of the prediction error for the processes excluding a each of the selected feature.
#' @return worsen: estimated increase in prediction error.
#' @return p_value: p values for proximity analysis.
#' @return selection_frequency: frequency of features in S being selected.
#' @return model_score: model_score for nextdoor analysis.
#' @return result_table: organized result table.
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom MASS mvrnorm
#' @examples
#' require(glmnet)
#' data(prostateCancerData)
#' data_train = prostateCancerData$train
#' x = data_train$feature
#' y = data_train$response
#' nams = data_train$names; n=length(y)
#' set.seed(48)
#' cv_glm = cv.glmnet(x = x, y = y, keep = TRUE, family = "gaussian", standardize = FALSE, nfolds = 10, nlambda = 30)
#' R1 = nextdoor.glmnet(x = x, y = y, nams = nams, family = "gaussian", cv_glm = cv_glm, standardize = FALSE, alpha = .1,B = 1000, B1 = 20)
#' print(R1, digits =3)
#' @export
nextdoor.glmnet <- function(x, y, cv_glm, nams = NULL, family = "gaussian",  lossfun = NULL, standardize = T,
K = 100, B = 1000, alpha = 0.1, epsilon = 0.05^2, epsilon2 =0.05^2,
selectionType = 0, Bindex = NULL, pv = TRUE,  rescale = TRUE,
score = TRUE, B1 = 50, Bindex1 = NULL,trace = TRUE){
    n = length(y); p = ncol(x);lambda = cv_glm$lambda; foldid = cv_glm$foldid; ypred = cv_glm$fit.preval[,1:length(lambda)]
    if(is.null(lossfun)){
      errors0 = apply(ypred, 2, deviances, y = y, family = family, model0 = cv_glm)
    }else{
      errors0 = apply(ypred, 2, lossfun, y = y)
    }
    if(is.null(nams)){
      nams = as.character(rep(1,p))
    }
    index0 = getIndex(errors0 = errors0,cv_glm = cv_glm,alpha = alpha, epsilon =epsilon, selectionType = selectionType)
    S = unlist(predict(cv_glm, s = lambda[index0], type = "n"))
    one_sds = cv_glm$cvsd
    lambda_extra = c(lambda, min(lambda)/2)
    model0 = cv_glm
    models = list(); errors = list()
    deBias_error0 = NULL; deBias_errors = NULL
    p_value = NULL; model_score = NULL; selection_frequency = NULL; worsen = NULL; result_table = NULL
    if(is.null(S)){
        return(list(model0 = model0, models = models, index = index0, S = S,
        deBias_error0 = deBias_error0, deBias_errors = deBias_errors,
        p_value = p_value, selection_frequency = selection_frequency,
        model_score = model_score))
    }else{
        for(s in S){
            Rs = train_model(x = x[,-s], y = y, family=family, foldid = foldid, lambda = lambda,lambda_extra = lambda_extra, lossfun = lossfun)
            errors[[s]] = Rs$errors0
            if(ncol(errors[[s]]) > ncol(errors0)){
                errors[[s]] = errors[[s]][,1:ncol(errors0)]
            }else if(ncol(errors[[s]]) < ncol(errors0)){
                for(kk in 1:(ncol(errors0) - ncol(errors[[s]]))){
                    errors[[s]]  = cbind(errors[[s]] , errors[[s]][,ncol( errors[[s]])])
                }
            }
            models[[s]] = Rs$model
        }
        R1 = nextdoor(errors0 = errors0, errors = errors, S = S, nams = nams, K = K, B =B,
        alpha = alpha, epsilon = epsilon, epsilon2 =epsilon2,
        Bindex = Bindex ,pv = pv,  rescale = rescale,
        selectionType = selectionType, one_sds = one_sds,
        trace = trace)
        p_value = R1$p_value
        debiased_errors0 = R1$debiased_errors0
        debiased_errors = R1$debiased_errors
        worsen = R1$worsen
        if(trace){
            print(paste(nams[S],":", p_value))
        }
        if(score){
            if(is.null(Bindex1)){
                Bindex1 = fastBoostrap(n = n, B = B1)+1
            }else{
                B1 = ncol(Bindex1)
            }
            counts = rep(0,length(S))
            print("start estimating the selection frequency.")
            print("########Bootstrap Iteration#############")
            for(b in 1:B1){
                if(trace){print(b)}
                xb = x[Bindex1[,b],]
                yb = y[Bindex1[,b]]
                Rb = train_model2(x = xb, y = yb, family=family, foldid = foldid, lambda = lambda, lambda_extra = lambda_extra, 
                                  standardize = standardize,lossfun = lossfun,epsilon = epsilon,alpha = alpha, selectionType = selectionType)

                counts = counts+ S%in%Rb$S
            }
            print("########Bootstrap Iteration#############")
            counts = counts/B1
            selection_frequency = counts
            model_score = p_value/selection_frequency
        }
        result_table <- create_table(model0 = model0, models = models, nams = nams,
        errors0 = errors0, errors = errors,
        debiased_errors0 = debiased_errors0,
        debiased_errors = debiased_errors,
        index = index0, S = S, p = ncol(x),
        pvalue = p_value, selection_frequency = selection_frequency,
        model_score = model_score)
        return_obj <-list(model0 = model0,  models = models,
                          index = index0, S = S,
                          errors0 = errors0, errors = errors,
                          debiased_errors0 = debiased_errors0, debiased_errors = debiased_errors,
                          worsen = worsen, p_value = p_value,
                          selection_frequency = selection_frequency, model_score = model_score,
                          result_table = result_table)
        class(return_obj) <- "nextdoor"
        
        return(return_obj)
    }
}

