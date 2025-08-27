#' Sinkhorn transform matrix
#'
#' @param V input matrix
#' @param max_iter Maximum iterations of the Sinkhorn scaling. Default is 20 iterations.
#' @param iter_start_check From which iteration should the function checks the convergence. By default, check started from iteration 5.
#' @param check_every_iter How offeten should we check the convergence. The default is check every 3 iterations
#' @param epsilon The tolerance for convergece. Default value is same as R's built in `all.equal` function.
#' @return scaling object (find it in dso$st$scaling)
#' @export
sinkhorn_scale <- function(
  V,
  max_iter = 20L,
  iter_start_check = 5L,
  check_every_iter = 3L,
  epsilon = sqrt(.Machine$double.eps)
) {
  efficient_sinkhorn(V, max_iter, iter_start_check, check_every_iter, epsilon)
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
  rownames(extended_sinkhorn_scale$V_row) <- rownames(V)
  colnames(extended_sinkhorn_scale$V_row) <- colnames(V)
  rownames(extended_sinkhorn_scale$V_col) <- rownames(V)
  colnames(extended_sinkhorn_scale$V_col) <- colnames(V)
  rownames(extended_sinkhorn_scale$H_row) <- rownames(H)
  colnames(extended_sinkhorn_scale$H_row) <- colnames(H)
  rownames(extended_sinkhorn_scale$H_col) <- rownames(H)
  colnames(extended_sinkhorn_scale$H_col) <- colnames(H)
  rownames(extended_sinkhorn_scale$H_row) <- rownames(H)
  colnames(extended_sinkhorn_scale$H_row) <- colnames(H)
  rownames(extended_sinkhorn_scale$W_col) <- rownames(W)
  colnames(extended_sinkhorn_scale$W_col) <- colnames(W)

  return(list(
    V_row = extended_sinkhorn_scale$V_row,
    V_col = extended_sinkhorn_scale$V_col,
    W_row = extended_sinkhorn_scale$W_row,
    W_col = extended_sinkhorn_scale$W_col,
    H_row = extended_sinkhorn_scale$H_row,
    H_col = extended_sinkhorn_scale$H_col,
    D_vs_row =  extended_scaling_result$D_vs_row,
    D_vs_col =  extended_scaling_result$D_vs_col,
    D_hs_row =  extended_scaling_result$D_hs_row,
    D_ws_col =  extended_scaling_result$D_hs_col
  ))
}