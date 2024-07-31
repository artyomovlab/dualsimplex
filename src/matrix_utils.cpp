#include "matrix_utils.h"

arma::mat correctByNorm(arma::mat& X) {
    arma::vec norm_(X.n_rows);
    for (int k = 0; k < X.n_rows; k++) {
        norm_[k] = arma::norm(X.row(k), 2);
    }
    arma::mat B = diagmat(1 / norm_) * X;
    B.elem(arma::find_nonfinite(B)).zeros();
    return B;
}

arma::rowvec find_cosine(const arma::mat& X) {
    arma::mat Y = arma::trans(X) * X;
    arma::mat res = Y / (arma::sqrt(arma::diagvec(Y)) * arma::trans(arma::sqrt(arma::diagvec(Y))));
    res.diag(0).fill(0);
    arma::rowvec X_dist = arma::max(res, 0);
    return X_dist;
}

arma::mat getDoubleProjection(const arma::mat& V, const arma::mat& R, const arma::mat& S) {
    arma::mat T = S.t() * S * V * R.t() * R;
    return T;
}

arma::uword getNegative(arma::mat X) {
    arma::uvec q1 = arma::find(X < 0);
    arma::vec B = arma::conv_to<arma::vec>::from(q1);
    return B.n_elem;
}

double getSum(arma::mat X, arma::mat M) { return arma::accu(X) / M.n_rows; }


Rcpp::List getNonnegativeLowRankApproximationWithSVD(const arma::mat& X,  
                                                     const int rank,
                                                     const int iterations,
                                                     const double left) {
  // truncated svd
  arma::mat Ur;
  arma::vec Sr;
  arma::mat Vr;
  arma::mat Yi;
  arma::mat errors_statistics(iterations, 2, arma::fill::zeros);
  svd(Ur,Sr,Vr,X);
  Ur = Ur.head_cols(rank);
  Vr = Vr.head_cols(rank);
  Sr = Sr.head(rank);
  Yi = Ur * arma::diagmat(Sr) * Vr.t();
  for (int i = 0; i < iterations; i++) {
    Yi.elem(arma::find(X < left)).fill(left);
    svd(Ur,Sr,Vr,Yi);
    Ur = Ur.head_cols(rank);
    Vr = Vr.head_cols(rank);
    Sr = Sr.head(rank);
    Yi = Ur * arma::diagmat(Sr) * Vr.t();
    // get statistics values
    // frobenius norm of negative elements
    double fro_norm = arma::norm( Yi.elem(arma::find(X < 0)), "fro" );
    // number of negatives
    arma::uword neg_count = static_cast<double>(getNegative(Yi));
    errors_statistics.row(i) = arma::rowvec{fro_norm, neg_count};
  }
   return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                             Rcpp::Named("errors") = errors_statistics);
}

//
//Rcpp::List getNonnegativeLowRankApproximationWithTangent(const arma::mat& X,
//                                                         const int rank,
//                                                         const int iterations,
//                                                         const double left){
//
//}
//
//
//Rcpp::List getNonnegativeLowRankApproximationWithHMT(const arma::mat& X,
//                                                     const int rank,
//                                                     const int p,
//                                                     const int k,
//                                                     const double rho,
//                                                     const int iterations,
//                                                     const double left){
//
//}
//
//Rcpp::List getNonnegativeLowRankApproximationWithTropp(const arma::mat& X,
//                                                       const int rank,
//                                                       const int k,
//                                                       const int l,
//                                                       const double rho,
//                                                       const int iterations,
//                                                       const double left){
//
//}
//
//Rcpp::List getNonnegativeLowRankApproximationWithGN(const arma::mat& X,
//                                                    const int rank,
//                                                    const int l,
//                                                    const double rho,
//                                                    const int iterations,
//                                                    const double left) {
//
//}
