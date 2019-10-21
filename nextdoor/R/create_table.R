#' Create the result table for glmnet models
#' @param model0 models for the original procedure
#' @param models models for the proximal procedures
#' @param nams a length p vector containing feature names.
#' @param errors0 the n by m errors of the original model sequence
#' @param errors a list of n by m errors of the proximal model sequences
#' @param debiased_errors0 debiased cross-validation/test error estimate for selected model in the original model sequence
#' @param debiased_errors debiased cross-validation/test error estimates for selected model in the alternative model sequences
#' @param index penalty index for those representative models
#' @param S selected features
#' @param p size of full feature dimension
#' @param pvalue model pvalues, default is NULL
#' @param selection_frequency selection frequencies, default is NULL
#' @param model_score model scores, default is NULL
#' @return resultTable: organized result table for proximity analysis
create_table <- function(model0, models, nams,
errors0, errors,
debiased_errors0, debiased_errors,
index, S, p,
pvalue, selection_frequency, model_score){
    ##get all features
    features0 <- unlist(predict(model0, s = model0$lambda[index], type = 'n'))
    features <- c()
    for(j in S){
        a = c(1:p)[-j]
        addded_features = unlist(predict(models[[j]], s = model0$lambda[index], type = 'n'))
        features = union(features,setdiff(addded_features, features0))
    }
    features = c(features0, features)
    results = data.frame(matrix(NA, ncol = length(S)+1, nrow = length(features)+5))
    colnames(results) = c("original", nams[S])
    results[1+length(features),1] = mean(errors0[,index])
    results[2+length(features),1] = debiased_errors0
    jj = 1
    for(j in S){
        jj = jj+1
        results[1+length(features), jj] = mean(errors[[j]][,index])
        results[2+length(features),jj] = debiased_errors[[j]]
    }
    if(!is.null(selection_frequency)){
        results[3+length(features),-1] = selection_frequency
    }
    if(!is.null(pvalue)){
        results[4+length(features),-1] = pvalue
    }
    if(!is.null(model_score)){
        results[5+length(features),-1] = model_score
    }
    rownames(results)= c(nams[features], "cv_error", "debiased_error",
    "selection_frequency","model_pvalue",
    "model_score")
    colnames(results) = c("original", nams[S])
    coefj =  coef(model0, s = model0$lambda[index])[-1,1]
    for(j in 1:length(features)){
        if(coefj[features[j]] != 0){
            results[j,1] = coefj[features[j]]
        }
    }
    col_id = 1
    for(j in S){
        coefj = rep(0, p)
        coefj[-j] = coef(models[[j]], s = model0$lambda[index])[-1,1]
        for(jj in 1:length(features)){
            if(coefj[features[jj]] != 0){
                results[jj,col_id+1] = coefj[features[jj]]
            }
        }
        col_id = col_id+1
    }
    ll= order(-results[length(features)+2,-1])
    results[,-1] = results[,-1][,ll]
    colnames(results)[-1] = colnames(results)[-1][ll]
    results[1:length(S),] = results[ll,]
    rownames(results)[1:length(S)] = rownames(results)[1:length(S)][ll]
    return(results)
}

