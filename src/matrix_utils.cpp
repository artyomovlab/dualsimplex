#include "matrix_utils.h"

arma::mat correctByNorm(arma::mat& X) {
    arma::vec norm_(X.n_rows);
    for (unsigned int k = 0; k < X.n_rows; k++) {
        norm_(k) = arma::norm(X.row(k), 2);
    }
    arma::mat B = diagmat(1 / norm_) * X;
    B.elem(arma::find_nonfinite(B)).zeros();
    return B;
}

double  cosine_similarity(const arma::rowvec& a,const arma::rowvec& b) {
    double Y = arma::dot(a,  b);
    double res_dist = Y / ( arma::norm(a,2) *  arma::norm(b,2));
    return res_dist;
}

arma::rowvec cosine_between_rows(const arma::mat& X) {
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


arma::mat get_relative_coordinates(const arma::mat& projected_points,const arma::mat& solution_points) {
    arma::mat coefficients(projected_points.n_rows, solution_points.n_rows, arma::fill::zeros);
    double main_determinant = std::abs(arma::det(solution_points));
    double target_determinant = 1;
    arma::mat target_matrix;
    for (unsigned int i = 0; i < projected_points.n_rows; i++) {
        for (unsigned int vertex = 0; vertex < solution_points.n_rows; vertex++ ) {
            target_matrix = solution_points;
            target_matrix.row(vertex) = projected_points.row(i);
            target_determinant = std::abs(arma::det(target_matrix));
            coefficients(i, vertex) = target_determinant/main_determinant;
        }

    }
    // now normalize row
    coefficients.each_col() %=  (1 / arma::sum(coefficients, 1));
    return coefficients;
}


arma::mat get_relative_coordinates_closest(const arma::mat& projected_points,const arma::mat& solution_points) {
    arma::mat coefficients(projected_points.n_rows, solution_points.n_rows, arma::fill::zeros);
    double main_determinant = arma::det(solution_points);
    double target_determinant = 1.0;
    arma::mat target_matrix;
    for (unsigned int i = 0; i < projected_points.n_rows; i++) {
        for (unsigned int vertex = 0; vertex < solution_points.n_rows; vertex++ ) {
            target_matrix = solution_points;
            target_matrix.row(vertex) = projected_points.row(i);
            target_determinant = arma::det(target_matrix);
            coefficients(i, vertex) = target_determinant/main_determinant;
        }
    }
    // negative elements set to 0
    //Rcpp::Rcout << coefficients << "\n";

    coefficients.elem(arma::find(coefficients < 0)).fill(0);
    // now normalize row
    coefficients.each_col() %=  (1 / arma::sum(coefficients, 1));
    return coefficients;
}




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
  arma::urowvec approximation_fro_norm(iterations, arma::fill::zeros);
  arma::urowvec normalized_feature_movements(iterations, arma::fill::zeros);


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
    double fro_distance = arma::norm(X-Yi, "fro");
    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    approximation_fro_norm(i) = fro_distance;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics,
                              Rcpp::Named("approx_error") = approximation_fro_norm);
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
    arma::urowvec approximation_fro_norm(iterations, arma::fill::zeros);
    arma::urowvec normalized_feature_movements(iterations, arma::fill::zeros);


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
    double fro_distance = arma::norm(X-Yi, "fro");

    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    approximation_fro_norm(i) = fro_distance;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics,
                              Rcpp::Named("approx_error") = approximation_fro_norm);


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
    arma::urowvec approximation_fro_norm(iterations, arma::fill::zeros);

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
    double fro_distance = arma::norm(X-Yi, "fro");
    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    approximation_fro_norm(i) = fro_distance;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics,
                              Rcpp::Named("approx_error") = approximation_fro_norm);

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
  arma::rowvec approximation_fro_norm(iterations, arma::fill::zeros);
  arma::rowvec normalized_feature_movements(iterations, arma::fill::zeros);

  svd(Ur,Sr,Vr,X);
  Ur = Ur.head_cols(rank); // m*r
  Vr = Vr.head_cols(rank); // n*r
  Sr = Sr.head(rank); // r
  Yi = Ur * arma::diagmat(Sr) * Vr.t(); // m*n
  arma::mat attention_matrix = arma::zeros(size(Yi));    
  int m = X.n_rows; 
  int n = X.n_cols;
  Im.eye(m, m); // m*m
  In.eye(n, n); // n*n
  for (int i = 0; i < iterations; i++) {
    attention_matrix.elem(arma::find(Yi < left)) += 1;
    if (left >= 0) {
        Yi.elem(arma::find(Yi < left)).fill(left);
    }
    else if (left == -1) {
        // mean for the whole matrix
        double mean_value = arma::mean(Yi.elem(arma::find(Yi >= 0)));
        Yi.elem(arma::find(Yi < 0)).fill(mean_value);
    } else if (left == -2) {
        // median for the whole matrix
        double mean_value = arma::median(Yi.elem(arma::find(Yi >= 0)));
        Yi.elem(arma::find(Yi < 0)).fill(mean_value);
    } else if (left == -3) {
        //mean for specific gene
        arma::uvec indices = arma::find(Yi < 0);
        arma::umat mutrix_indices = arma::ind2sub( size(Yi), indices );
        arma::urowvec row_indicies = mutrix_indices.row(0);
        arma::mat res_mat = Yi.rows(row_indicies);
        res_mat.elem(arma::find(res_mat < 0)).fill(0);
        arma::vec row_means = arma::mean(res_mat, 1);
        Yi.elem(arma::find(Yi < 0)) = row_means;
    }

    if (right > 0) {
        Yi.elem(arma::find(Yi > right)).fill(right);
    }
    G1 = Ur.t() * Yi; // r*m x m*n -> r*n
    G2 = (Im - (Ur * Ur.t())) * Yi * Vr; // (m*m - (m*r x r*m)) * m*n * n*r -> m*r

    arma::qr_econ(Q1, R1, ((In-(Vr * Vr.t()))* G1.t())); //  QR((n*n - (n*n))* n*r ) -> QR(n*r) ->  (n*r) (r*r)
    arma::qr_econ(Q2, R2, G2); // (m*r) (r*r)
    Z = arma::join_rows(G1 * Vr, R1.t()); // join_rows((r*n x n*r), r*r) ->  r*2r
    Z = arma::join_cols(Z, arma::join_rows(R2, arma::zeros(size(R2))));  // 2r*2r
    svd(U2r,S2r,V2r,Z); //   
    U2r = U2r.head_cols(rank); // 2r*r
    V2r = V2r.head_cols(rank); // 2r*r
    S2r = S2r.head(rank); //r
    Ur = arma::join_rows(Ur, Q2) * U2r; // jc(m*r, m*r)x 2r*r -> m*r
    Vr = arma::join_rows(Vr, Q1) * V2r; // jc(n*r, n*r) x 2r*r -> n*r
    Yi = Ur * arma::diagmat(S2r) * Vr.t();
    // get statistics values
    // frobenius norm of negative elements
    double fro_norm = arma::norm( Yi.elem(arma::find(Yi < 0)), "fro" );
    double fro_distance = arma::norm(X-Yi, "fro");
    double mean_gene_movement = arma::mean(arma::vecnorm((arma::normalise(X, 1, 1) - arma::normalise(Yi, 1, 1)), 1, 1));
    // number of negatives
    arma::uword neg_count = getNegative(Yi);
    frobenius_statistics(i) = fro_norm;
    neg_elements_statistics(i) = neg_count;
    approximation_fro_norm(i) = fro_distance;
    normalized_feature_movements(i) = mean_gene_movement;
    }
    return Rcpp::List::create(Rcpp::Named("newX") = Yi,
                              Rcpp::Named("attention") = attention_matrix,
                              Rcpp::Named("frobenius_neg_norm") = frobenius_statistics,
                              Rcpp::Named("neg_count") = neg_elements_statistics,
                              Rcpp::Named("approx_error") = approximation_fro_norm,
                              Rcpp::Named("gene_movements") = normalized_feature_movements
                              );
}