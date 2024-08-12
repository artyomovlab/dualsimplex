#include "matrix_utils.h"

arma::mat correctByNorm(arma::mat& X) {
    arma::vec norm_(X.n_rows);
    for (unsigned int k = 0; k < X.n_rows; k++) {
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
                                                     const double left,
                                                     const double right) {
  // truncated svd
  arma::mat Ur;
  arma::vec Sr;
  arma::mat Vr;
  arma::mat Yi;
  arma::rowvec frobenius_statistics(iterations, arma::fill::zeros);
  arma::urowvec neg_elements_statistics(iterations, arma::fill::zeros);    
  svd(Ur,Sr,Vr,X);
  Ur = Ur.head_cols(rank);
  Vr = Vr.head_cols(rank);
  Sr = Sr.head(rank);
  Yi = Ur * arma::diagmat(Sr) * Vr.t();
  for (int i = 0; i < iterations; i++) {
    Yi.elem(arma::find(Yi < left)).fill(left);
    if (right > 0) {
        Yi.elem(arma::find(Yi > right)).fill(right);
    }
    svd(Ur,Sr,Vr,Yi);
    Ur = Ur.head_cols(rank);
    Vr = Vr.head_cols(rank);
    Sr = Sr.head(rank);
    Yi = Ur * arma::diagmat(Sr) * Vr.t();
    // get statistics values
    // frobenius norm of negative elements
    double fro_norm = arma::norm( Yi.elem(arma::find(Yi < 0)), "fro" );
    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics);
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
Rcpp::List getNonnegativeLowRankApproximationWithHMT(const arma::mat& X,
                                                     const int rank,
                                                     const int p,
                                                     const int k,
                                                     const int iterations,
                                                     const double left,
                                                     const double right) {

    arma::rowvec frobenius_statistics(iterations, arma::fill::zeros);
    arma::urowvec neg_elements_statistics(iterations, arma::fill::zeros);                                                  
    arma::mat Psi;
    arma::mat Z1, Z2;
    arma::mat Q, R;
    // initial truncated decomposition
    arma::mat Ur;
    arma::vec Sr;
    arma::mat Vr;
    arma::mat Yi;
    arma::svd(Ur,Sr,Vr,X);
    Ur = Ur.head_cols(rank);
    Vr = Vr.head_cols(rank);
    Sr = Sr.head(rank);
    Yi = Ur * arma::diagmat(Sr) * Vr.t();
    int n = X.n_cols;
    for (int i = 0; i < iterations; i++) {
    Yi.elem(arma::find(Yi < left)).fill(left);
    if (right > 0) {
        Yi.elem(arma::find(Yi > right)).fill(right);
    }
    //generate psi matrix (n, k, norm_dist, rho)
    Psi =  arma::randn(n, k);
    Z1 = Yi * Psi;
    // do QR decomposition
    arma::qr_econ(Q, R, Z1);
    for (int j = 0; j < p; j++) {
        Z2 = Q.t() * Yi;
        arma::qr_econ(Q, R, Z2.t());
        Z1 = Yi * Q;
        arma::qr_econ(Q, R, Z1);
    }
    Z2 = Q.t() * Yi;
    arma::svd(Ur,Sr,Vr,Z2);
    Ur = Ur.head_cols(rank);
    Vr = Vr.head_cols(rank);
    Sr = Sr.head(rank);
    Ur = Q * Ur;
    Yi = Ur * arma::diagmat(Sr) * Vr.t();
    // get statistics values
    // frobenius norm of negative elements
    double fro_norm = arma::norm( Yi.elem(arma::find(Yi < 0)), "fro" );
    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics);


}
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
Rcpp::List getNonnegativeLowRankApproximationWithGN(const arma::mat& X,
                                                   const int rank,
                                                   const int l,
                                                   const int iterations,
                                                   const double left,
                                                   const double right) {
    arma::rowvec frobenius_statistics(iterations, arma::fill::zeros);
    arma::urowvec neg_elements_statistics(iterations, arma::fill::zeros);    
    arma::mat Psi, Phi;
    arma::mat Z, W, V, U;
    arma::mat Q, R;
    // initial truncated decomposition
    arma::mat Ur;
    arma::vec Sr;
    arma::mat Vr;
    arma::mat Yi;
    Yi = X;
    int m = X.n_rows;
    int n = X.n_cols;

    for (int i = 0; i < iterations; i++) {
    Yi.elem(arma::find(Yi < left)).fill(left);
    if (right > 0) {
        Yi.elem(arma::find(Yi > right)).fill(right);
    }
    //generate psi matrix (n, rank)
    Psi =  arma::randn(n, rank);
    //generate phi matrix (l, m)
    Phi = arma::randn(l, m);
    Z = Yi * Psi;
    W = Phi * Z;
    // do truncated QR decomposition
    arma::qr_econ(Q, R, W);
    // R = R.head_rows(std::min(l, rank));
    // Q = Q.head_cols(std::min(l, rank));
    V =  (Phi * Yi).t() * Q;
    U = Z * arma::inv(R);
    Yi = U * V.t();

    // get statistics values
    // frobenius norm of negative elements
    double fro_norm = arma::norm( Yi.elem(arma::find(Yi < 0)), "fro" );
    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics);

}


Rcpp::List getNonnegativeLowRankApproximationWithTangentMethod(const arma::mat& X,  
                                                     const int rank,
                                                     const int iterations,
                                                     const double left,
                                                     const double right) {
  // truncated svd
  arma::mat Ur, U2r;
  arma::vec Sr, S2r;
  arma::mat Vr, V2r;
  arma::mat Yi;
  arma::mat Im, In;
  arma::mat G1, G2, Q1, Q2, R1, R2, Z;
  arma::rowvec frobenius_statistics(iterations, arma::fill::zeros);
  arma::urowvec neg_elements_statistics(iterations, arma::fill::zeros);    
  svd(Ur,Sr,Vr,X);
  Ur = Ur.head_cols(rank);
  Vr = Vr.head_cols(rank);
  Sr = Sr.head(rank);
  Yi = Ur * arma::diagmat(Sr) * Vr.t();
  int m = X.n_rows;
  int n = X.n_cols;
  Im.eye(m, m);
  In.eye(n, n);
  for (int i = 0; i < iterations; i++) {
    Yi.elem(arma::find(Yi < left)).fill(left);
    if (right > 0) {
        Yi.elem(arma::find(Yi > right)).fill(right);
    }
    G1 = Ur.t() * Yi;
    G2 = (Im - (Ur * Ur.t())) * Yi * Vr;

    arma::qr_econ(Q1, R1, (In-(Vr*Vr.t())* G1.t()));
    arma::qr_econ(Q2, R2, G2);
    Z = arma::join_rows(G1*Vr, R1.t());
    Z = arma::join_cols(Z, arma::join_cols(R2, arma::zeros(size(R2))));
    svd(U2r,S2r,V2r,Z);
    U2r = U2r.head_cols(rank);
    V2r = V2r.head_cols(rank);
    S2r = S2r.head(rank);
    Ur = arma::join_rows(Ur, Q2) * U2r;
    Vr = arma::join_rows(Vr, Q1) * V2r;
    Yi = Ur * arma::diagmat(Sr) * Vr.t();
    // get statistics values
    // frobenius norm of negative elements
    double fro_norm = arma::norm( Yi.elem(arma::find(Yi < 0)), "fro" );
    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics);
}