#' Print result summary of a nextdoor object
#'
#' @param x return from nextdoor.glmnet or nextdoor.
#' @param ... parameters for print.
#' @export
print.nextdoor <- function(x,...){
  if(!is.null(x['result_table'])){
    print(x[['result_table']],...)
  }else{
    stop("not an nextdoor object.")
  }
}
