#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
arma::uword minIndex2(arma::vec vec, int index0, arma::vec sds){
    bool cont = true;
    int j = -1;
    double t =vec(index0)+sds(index0);
    while(cont){
      j += 1;
      if(vec(j) < t){
        cont = false;
      }
    }
    return j;
}

arma::mat mvrnormArma(int m, int K, arma::mat sigmaLeft) {
    arma::mat Y = arma::randn(m, K);
    return sigmaLeft*Y;
}
arma::mat normArma(int m, int K, double scale){
    arma::mat Y = arma::randn(m, K);
    return Y *scale;
}

// [[Rcpp::export]]
arma::vec boostrap_proximity_random(arma::mat errors0, arma::mat errors,
                         arma::vec errors0Mean, arma::vec errorsMean,
                         arma::mat Bindex, double epsilon, double epsilon2, double alpha, arma::mat sigmaLeft,
                         int K, int selectionType, arma::vec sds){
  arma::uword sample_size = errors.n_rows;
  arma::uword lambda_size = errors.n_cols;
  arma::uword B = Bindex.n_cols;
  arma::vec deltas = arma::zeros<arma::vec>(B);
  double trueDiff = 0;
  for(int b = 0; b < B; ++b){
    arma::vec errors0MeanB = arma::zeros<arma::vec>(lambda_size);
    arma::vec errorsMeanB = arma::zeros<arma::vec>(lambda_size);
    arma::mat errors0B = arma::zeros<arma::mat>(sample_size, lambda_size);
    arma::mat errorsB = arma::zeros<arma::mat>(sample_size, lambda_size);
    for(arma::uword i = 0; i < sample_size; i++){
      for(arma::uword j = 0; j < lambda_size; j ++){
        errors0B(i, j) = errors0(Bindex(i,b),j);
        errorsB(i, j) = errors(Bindex(i,b),j);
        errors0MeanB(j) = errors0MeanB(j)+errors0B(i, j)/sample_size;
        errorsMeanB(j) = errorsMeanB(j)+errorsB(i, j)/sample_size;
      }
    }
    arma::mat ADD1 = normArma(lambda_size*2, K, epsilon);
    arma::mat ADD2 = mvrnormArma(lambda_size*2, K, sigmaLeft);
    arma::mat cv_random0 = arma::repmat(errors0MeanB, 1, K) +ADD1.rows(0, lambda_size-1);
    arma::mat cv_random_addB0 = cv_random0+ sqrt(alpha)*ADD2.rows(0, lambda_size-1);
    arma::mat cv_random_minusB0 = cv_random0 - sqrt(1/alpha)*ADD2.rows(0, lambda_size-1);
    arma::uword idx0 = 0;
    for(arma::uword k = 0; k < K; k++){
      idx0 = cv_random_addB0.col(k).index_min();
      if(selectionType == 1){idx0 = minIndex2(cv_random_addB0.col(k), idx0, sds);}
      deltas[b] += errorsMeanB(idx0)+ADD1(lambda_size+idx0,k)- sqrt(1/alpha)*ADD2(lambda_size+idx0, k) -
      cv_random_minusB0(idx0,k);
    }
    trueDiff += ((errorsMean(idx0)-errors0Mean(idx0)));
    deltas[b] = deltas[b]/K+ADD1(0,0)*epsilon2/epsilon;
  }
  trueDiff  = trueDiff/B;
  for(int b = 0; b < B; b++){
    deltas[b] = deltas[b]-trueDiff;
  }
  return deltas;
}


