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
