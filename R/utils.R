calc_partial_dist <- function(data_proj_single, with_dims) {
  data_part <- data_proj_single
  data_part[, -with_dims] <- 0
  return(euclidean_dist(data_part))
}

euclidean_dist <- function(data) {
  return(sqrt(rowSums(data^2)))
}

#' Linearize log transformed gene expression data
#'
#' @param ge expression set matrix
#' @return linearized data
#' @export
linearize_dataset <- function(ge) {
  if (is_logscale(ge))
    return(2^ge - 1)
  return(ge)
}

#' Remove duplicated genes
#'
#' @param ge expression set matrix
#' @return clean data
#' @export
replace_duplicate_genes <- function(ge) {
  ge <- ge[order(rownames(ge), -abs(rowSums(ge))), ]
  ge <- ge[!duplicated(rownames(ge)), ]
  return(ge)
}

#' Log scale data
#'
#' @param ge expression set matrix
#' @return log scaled data
#' @export
log_dataset <- function(ge) {
  if (is_logscale(ge))
    return(ge)
  return(log2(ge + 1))
}

rownorm <- function(ge) {
  return(ge / rowSums(ge))
}

#' Check if log scale
#'
#' @param x expression set matrix
#' @return TRUE or FALSE
#' @export
is_logscale <- function(x) {
  qx <- quantile(as.numeric(x), na.rm = T)
  if (qx[5] - qx[1] > 100 || qx[5] > 100) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Remove zero rows
#'
#' @param ge expression set matrix
#' @return clean data
#' @export
remove_zero_rows <- function(ge) {
  zero_rows <- apply(ge, 1, function(x) all(x == 0))
  return(ge[!zero_rows, ])
}

#' Remove zero cols
#'
#' @param ge expression set matrix
#' @return clean data
#' @export
remove_zero_cols <- function(ge) {
  zero_cols <- apply(ge, 2, function(x) all(x == 0))
  return(ge[,!zero_cols ])
}


#' Remove zero rows
#'
#' @param data expression set matrix
#' @param expr predicate to filter the data
#' @param ... any other parameters for dplyr::filter
#' @return clean data
#' @importFrom rlang parse_expr
#' @export
filter_str <- function(data, expr, ...) {
  return(dplyr::filter(data, eval(rlang::parse_expr(expr)), ...))
}

#' Reade genes form a list
#'
#' @param path path to file
#' @return R list
#' @import stringr
#' @export
read_gene_list <- function(path) {
  genes_string <- readChar(path, file.info(path)$size)
  genes_list <- stringr::str_split(genes_string, " ")[[1]]
  return(genes_list)
}

rasterize_if_needed <- function(plot) {
  plot <- if(getOption("dualsimplex-rasterize", default=FALSE)) ggrastr::rasterise(plot, dpi=600) else plot
  return(plot)
}

#' Calculate zero_distance and plane_distance with truncated SVD
#'
#' @param approximated a low-rank approximation by SVD
#' @param original the original data matrix
#' @param residual the residual matrix
#' @param margin direction for distance calculation. 1 is for rows, 2 is for columns.
#' 
#' @return R list containing zero_distance and plane_distance
calc_dist_from_truncated_svd <- function(approximated, original = NULL, residual = NULL, margin) {
  if (is.null(original) & is.null(residual)) stop("One of the 'original' or 'residual' matrix should be provided")
  if (is.null(residual)) residual <- original - approximated

  dist_fns  <- list(
    # margin = 1: by row
    row_dist_from_mean,

    # margin = 2: by column
    col_dist_from_mean
  )

  list(
    zero_distance = dist_fns[[margin]](approximated),
    plane_distance =   dist_fns[[margin]](residual)
  )
}

#' Calculate Euclidean distance around the mean for each row
#'
#'@param mtx a dense matrix
#'@importFrom matrixStats rowSds
#'@return  a numeric vector
row_dist_from_mean <- function(mtx) {
  # Notice that standard deviation estimates is devided by sqrt(N-1)
  matrixStats::rowSds(mtx) * sqrt(ncol(mtx) - 1)
}

#' Calculate Euclidean distance around the mean for each column
#'
#'@param mtx a dense matrix
#'@importFrom matrixStats colSds
#'@return  a numeric vector
col_dist_from_mean <- function(mtx) {
  # Notice that standard deviation estimates is devided by sqrt(N-1)
  matrixStats::colSds(mtx) * sqrt(nrow(mtx) - 1)
}