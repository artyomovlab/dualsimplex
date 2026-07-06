#' Sinkhorn transform matrix
#'
#' @param V input matrix
#' @param max_iter Maximum iterations of the Sinkhorn scaling. Default is 20 iterations.
#' @param iter_start_check From which iteration should the function checks the convergence. By default, check started from iteration 5.
#' @param check_every_iter How offeten should we check the convergence. The default is check every 3 iterations.
#' @param epsilon The tolerance for convergece. Default value is same as R's built in `all.equal` function.
#' @param return_scaled_matrix wether to return actual scaled matrix. By default FALSE for efficiency.
#' @return scaling object (find it in dso$st$scaling)
#' @export
sinkhorn_scale <- function(
  V,
  max_iter = 20L,
  iter_start_check = 5L,
  check_every_iter = 3L,
  epsilon = sqrt(.Machine$double.eps),
  return_scaled_matrix = FALSE
) {
  scaling <- efficient_sinkhorn(V, max_iter, iter_start_check, check_every_iter, epsilon)
  if (return_scaled_matrix) {
  scaling$V_row  <- sinkhorn_sweep_c(V = V,
                                     D_vs_row = scaling$D_vs_row,
                                     D_vs_col = scaling$D_vs_col,
                                     iter = scaling$iterations,
                                     do_last_step = 0)
  scaling$V_column  <- sinkhorn_sweep_c(V = V,
                                 D_vs_row = scaling$D_vs_row,
                                 D_vs_col = scaling$D_vs_col,
                                 iter = scaling$iterations,
                                 do_last_step = 1)
  }
  return(scaling)
}

#' reverse Sinkhorn transform scaled matrix
#'
#' Calculate all normalizing matrices and restore "unscaled" W and H
#' columns of D_ws_col and D_hs_row will contain history of D_w and D_h matrices
#' W_column, H_column will contain "unscaled" W and H
#'
#' @param H_row solution matrix H_ss produced by optimization algorithm (X*R)
#' @param W_col solution matrix W_gs produced by optimization algorithm (t(S)*Omega)
#' @param scaling dso$scaling object containing normalization matrices
#' @return unscaled object
#' @import Rcpp
#' @import RcppArmadillo
#' @export
reverse_sinkhorn <- function(H_row, W_col, scaling) {
  unscaled <- reverse_sinkhorn_c(H_row,
                                 W_col,
                                 scaling$D_vs_row,
                                 scaling$D_vs_col,
                                 scaling$iterations)
  
  dimnames(unscaled$H) <- dimnames(H_row)
  dimnames(unscaled$W) <- dimnames(W_col)
  dimnames(unscaled$H_row) <- dimnames(H_row)
  dimnames(unscaled$Dv_inv_W_row) <- dimnames(W_col)
  
  return(unscaled)
}

#' reverse Sinkhorn transform scaled matrix
#'
#' Calculate all normalizing matrices and restore "unscaled" W and H
#' columns of D_ws_col and D_hs_row will contain history of D_w and D_h matrices
#' W_column, H_column will contain "unscaled" W and H
#'
#' @param H_col solution matrix H_ss produced by optimization algorithm (X*R).
#' @param W_row solution matrix W_gs produced by optimization algorithm (t(S)*Omega).
#' @param scaling dso$scaling object containing normalization matrices.
#' @param enforce_sum_to_one_H algorithm will use NNLS to find H as close to column normalized as possible.
#' @param enforce_sum_to_one_V algorithm will ensure result V is column normalized (leading to H column normalizing).
#' @return unscaled object
#' @import Rcpp
#' @import RcppArmadillo
#' @export
clean_reverse_sinkhorn <- function(H_col, W_row, scaling, enforce_sum_to_one_H = 0, enforce_sum_to_one_V = 0) {
  unscaled <- clean_reverse_sinkhorn_c(
    H_col,
    W_row,
    scaling$D_vs_row,
    scaling$D_vs_col,
    scaling$iterations,
    enforce_sum_to_one_H = enforce_sum_to_one_H,
    enforce_sum_to_one_V = enforce_sum_to_one_V
    )

  dimnames(unscaled$H) <- dimnames(H_col)
  dimnames(unscaled$W) <- dimnames(W_row)
  # Alternative solution - W could be col-norm or H col-norm
  dimnames(unscaled$H_2) <- dimnames(H_col)
  dimnames(unscaled$W_2) <- dimnames(W_row)
  return(unscaled)
}


#' reverse Sinkhorn transform scaled matrix
#'
#' Calculate all normalizing matrices and restore "unscaled" W and H
#' columns of D_ws_col and D_hs_row will contain history of D_w and D_h matrices
#' W_column, H_column will contain "unscaled" W and H
#'
#' @param h_col solution matrix H_ss produced by optimization algorithm (X*R)
#' @param w_row solution matrix W_gs produced by optimization algorithm (t(S)*Omega)
#' @param scaling dso$scaling object containing normalization matrices
#' @param v_inf_row result V_ss matrix (you can get it with `get_V_row` function)
#' @return unscaled object
#' @import Rcpp
#' @import RcppArmadillo
#' @export
geometrical_reverse_sinkhorn <- function(h_col, w_row, scaling, v_inf_row) {
  unscaled <- geometrical_reverse_sinkhorn_c(h_col,
                                 w_row,
                                 scaling$D_vs_row,
                                 scaling$D_vs_col,
                                 v_inf_row,
                                 scaling$iterations)


  dimnames(unscaled$H) <- dimnames(h_col)
  dimnames(unscaled$W) <- dimnames(w_row)
  dimnames(unscaled$H_col) <- dimnames(h_col)
  dimnames(unscaled$W_corrected) <- dimnames(w_row)
  return(unscaled)
}

#' Wrapper to call reverse sinkhorn on dso object
#'
#' @param solution_scaled contain solution matrices $H_row and $W_col
#' @param scaling dso$scaling object containing normalization matrices
#' @return unscaled object
#' @import Rcpp
#' @import RcppArmadillo
#' @export
reverse_solution_sinkhorn <- function(solution_scaled, scaling) {
  return(reverse_sinkhorn(
    solution_scaled$H_row,
    solution_scaled$W_col,
    scaling
  ))
}

#' Wrapper to call reverse sinkhorn on dso object
#'
#' @param solution_scaled contain solution matrices $H_row and $W_col.
#' @param scaling dso$scaling object containing normalization matrices.
#' @param ... everything which can be passed to `clean_reverse_sinkhorn_c` function.
#' @return unscaled object
#' @import Rcpp
#' @import RcppArmadillo
#' @export
clean_reverse_solution_sinkhorn<- function(solution_scaled, scaling, ...) {
  return(clean_reverse_sinkhorn(
    solution_scaled$H_col,
    solution_scaled$W_row,
    scaling, ...
  ))
}

#' Wrapper to call geometrical sinkhorn on dso object
#'
#' @param solution_scaled contain solution matrices $H_row and $W_col
#' @param scaling dso$scaling object containing normalization matrices
#' @param V_inf_row V_ss matrix obtained during forward scaling (can get with `get_V_row()` function)
#' @return unscaled object
#' @import Rcpp
#' @import RcppArmadillo
#' @export
geometrical_reverse_solution_sinkhorn<- function(solution_scaled, scaling, V_inf_row) {
  return(geometrical_reverse_sinkhorn(
    solution_scaled$H_col,
    solution_scaled$W_row,
    scaling,
    V_inf_row
  ))
}


#' Extended sinkhorn scaling to track all the matrices produced.
#' This is extremely usefull if you want to track back the results and find "true" coordinates of W and H in sinkhorned space.
#'
#' @param V input matrix V  (V = WH)
#' @param W input matrix W
#' @param H input matrix H
#' @param n_iter exact number of iterations performed

#' @export
extended_sinkhorn_scale <- function(V, W, H, n_iter) {
  extended_scaling_result <- extended_sinkhorn(V, W, H, n_iter)
  rownames(extended_scaling_result$V_row) <- rownames(V)
  colnames(extended_scaling_result$V_row) <- colnames(V)
  rownames(extended_scaling_result$V_col) <- rownames(V)
  colnames(extended_scaling_result$V_col) <- colnames(V)
  rownames(extended_scaling_result$H_row) <- rownames(H)
  colnames(extended_scaling_result$H_row) <- colnames(H)
  rownames(extended_scaling_result$H_col) <- rownames(H)
  colnames(extended_scaling_result$H_col) <- colnames(H)
  rownames(extended_scaling_result$H_row) <- rownames(H)
  colnames(extended_scaling_result$H_row) <- colnames(H)
  rownames(extended_scaling_result$W_col) <- rownames(W)
  colnames(extended_scaling_result$W_col) <- colnames(W)

  return(list(
    V_row = extended_scaling_result$V_row,
    V_col = extended_scaling_result$V_col,
    W_row = extended_scaling_result$W_row,
    W_col = extended_scaling_result$W_col,
    H_row = extended_scaling_result$H_row,
    H_col = extended_scaling_result$H_col,
    D_vs_row =  extended_scaling_result$D_vs_row,
    D_vs_col =  extended_scaling_result$D_vs_col,
    D_hs_row =  extended_scaling_result$D_hs_row,
    D_ws_col =  extended_scaling_result$D_ws_col
  ))
}