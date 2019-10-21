#' Default Bootstrap
#' @param n number of samples
#' @param B number of Bootstrap repetitions
#' @return A: 0 indexes n by B Bootstrap index matrix
fastBoostrap <- function(n,B){
  aa0 = sample(0:(n-1), n*B, replace = T)
  A = matrix(aa0, nrow = n)
  A
}
