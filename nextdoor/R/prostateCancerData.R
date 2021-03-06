#' Proximity analysis example with prostate cancer data
#'
#' The data is stored as a list. The training data with name data_train and the test data
#' with name data_test. It consists of $n=67$ training observations and 30 test observations. There are eight predictors.
#' The response is the log PSA for men who had  prostate cancer surgery.
#'
#' @docType data
#'
#' @usage data(prostateCancerData)
#'
#' @format An list object with element "train" and "test", each of them is again a list object with three elements "feature", "response", "names"
#'
#' @keywords datasets
#'
#' @references Friedman, J., Hastie, T. & Tibshirani, R. (2001), The elements of statistical learning, Vol. 1, Springer series in statistics New York.
#'
#' @examples
#' require(glmnet)
#' data(prostateCancerData)
#' data_train = prostateCancerData$train
#' x = data_train$feature
#' y = data_train$response
#' nams = data_train$names; n=length(y)
#' set.seed(48)
#' cv_glm = cv.glmnet(x = x, y = y, keep = TRUE, family = "gaussian", standardize = FALSE, nfolds = 3, nlambda = 30)
#' R1 = nextdoor.glmnet(x = x, y = y, nams = nams, family = "gaussian", cv_glm = cv_glm, standardize = FALSE, alpha = .1,B = 1000, B1 = 20)
#' print(round(R1$result_table,3))
"prostateCancerData"



