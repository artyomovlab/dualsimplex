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
