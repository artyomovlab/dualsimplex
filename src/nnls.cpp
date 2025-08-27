#include "nnls.h"

arma::vec nnls_col(const arma::mat &A,
                   const arma::subview_col<double> &b,
                   int max_iter = 500,
                   double tol = 1e-6,
                   bool verbose = false) {
    arma::vec mu = -A.t() * b;
    arma::mat H = A.t() * A;
    arma::vec x(A.n_cols), x0(A.n_cols);
    x.fill(0);
    x0.fill(-9999);

    int i = 0;
    double tmp;
    while (i < max_iter && max(abs(x - x0)) > tol) {
        x0 = x;
        for (unsigned int k = 0; k < A.n_cols; k++) {
            tmp = x(k) - mu(k) / H.at(k, k);
            if (tmp < 0) tmp = 0;
            if (tmp != x(k)) mu += (tmp - x(k)) * H.col(k);
            x(k) = tmp;
        }
        ++i;
    }

    return x;
}


arma::vec nnls_nonzero_col(const arma::mat &A,
                   const arma::subview_col<double> &b,
                   int max_iter = 500,
                   double tol = 1e-6,
                   bool verbose = false) {
    arma::vec mu = -A.t() * b;
    arma::mat H = A.t() * A;
    arma::vec x(A.n_cols), x0(A.n_cols);
    x.fill(1e-5);
    x0.fill(-9999);

    int i = 0;
    double tmp;
    while (i < max_iter && max(abs(x - x0)) > tol) {
        x0 = x;
        for (unsigned int k = 0; k < A.n_cols; k++) {
            tmp = x(k) - mu(k) / H.at(k, k);
            if (tmp < 0) tmp = 1e-5;
            if (tmp != x(k)) mu += (tmp - x(k)) * H.col(k);
            x(k) = tmp;
        }
        ++i;
    }

    return x;
}
    
arma::mat nnls_C__(arma::mat A, arma::mat b, int max_iter, double tol) {
    // solving Ax = b, where x and b are both matrices
    if (A.n_rows != b.n_rows)
        throw std::invalid_argument("A and b must have the same number of rows.");
    arma::mat x(A.n_cols, b.n_cols);
    for (unsigned int i = 0; i < b.n_cols; i++) x.col(i) = nnls_col(A, b.col(i), max_iter, tol);

    return x;
}


arma::mat nnls_nonzero_C__(arma::mat A, arma::mat b, int max_iter, double tol) {
    // solving Ax = b, where x and b are both matrices
    if (A.n_rows != b.n_rows)
        throw std::invalid_argument("A and b must have the same number of rows.");
    arma::mat x(A.n_cols, b.n_cols);
    for (unsigned int i = 0; i < b.n_cols; i++) x.col(i) = nnls_nonzero_col(A, b.col(i), max_iter, tol);

    return x;
}
