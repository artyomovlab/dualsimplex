#' Sinkhorn transform matrix
#'
#' @param V input matrix
#' @param iterations number of iterations
#' @return scaling object (find it in dso$st$scaling)
#' @export
sinkhorn_scale <- function(V, iterations = 20) {
  scaling <- sinkhorn_scale_c(V, iterations)
  
  dimnames(scaling$V_row) <- dimnames(V)
  dimnames(scaling$V_column) <- dimnames(V)
  
  scaling$iterations <- iterations
  
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
