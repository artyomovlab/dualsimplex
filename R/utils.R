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
  qx <- stats::quantile(as.numeric(x), na.rm = T)
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
  if (requireNamespace('ggrastr', quietly = TRUE)) {
  plot <- if(getOption("dualsimplex-rasterize", default=FALSE)) ggrastr::rasterise(plot, dpi=600) else plot
  }
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

##### Metric functions #####

##### Metric functions #####

#' Calculate pearson correlation. From -1 to 1. Further away from zero is better.
#'
#'@param estimated_vector estimated row or column vector
#'@param true_vector row or column vector to compare with
#'@return  a value
#'@export
pearson_correlation_function <- function(estimated_vector, true_vector) {
    return(stats::cor(estimated_vector, true_vector, method ="pearson"))
}

#' Calculate 1 - squared(pearson) metric value.  From 0 to 1. Lower values are better.
#'
#'@param estimated_vector estimated row or column vector
#'@param true_vector row or column vector to compare with
#'@return  a value
#'@export
pearson_loss_function <- function(estimated_vector, true_vector) {
  return(1-pearson_correlation_function(estimated_vector, true_vector)^2)
}

#' Calculate spearman correlation. From -1 to 1. Further away from zero is better.
#'
#'@param estimated_vector estimated row or column vector
#'@param true_vector row or column vector to compare with
#'@return  a value
#'@export
spearman_correlation_function <- function(estimated_vector, true_vector) {
    return(stats::cor(estimated_vector, true_vector, method ="spearman"))
}

#' Calculate 1 - squared(spearman) metric value.  From 0 to 1. Lower values are better.
#'
#'@param estimated_vector estimated row or column vector
#'@param true_vector row or column vector to compare with
#'@return  a value
#'@export
spearman_loss_function <- function(estimated_vector, true_vector) {
  return(1-spearman_correlation_function(estimated_vector, true_vector)^2)
}

#' Calculate cosine similarity function.  From 0 to 1. Higher is better
#'
#'@param estimated_vector estimated row or column vector
#'@param true_vector row or column vector to compare with
#'@return  a value
#'@export
cosine_similarity_function <- function(estimated_vector, true_vector) {
  return(cosine_similarity(estimated_vector, true_vector))
}

#' Calculate 1 - cosine metric value.  From 0 to 1. Lower values are better.
#'
#'@param estimated_vector estimated row or column vector
#'@param true_vector row or column vector to compare with
#'@return  a value
#'@export
cosine_loss_function <- function(estimated_vector, true_vector) {
  return(1 - cosine_similarity(estimated_vector, true_vector))
}

#' Calculate RMSE metric value.  From 0 to 1. Lower values are better.
#'
#'@param estimated_vector estimated row or column vector
#'@param true_vector row or column vector to compare with
#'@return  a value
#'@export
rmse_loss_function <- function(estimated_vector, true_vector) {
  return(sqrt((1/(length(estimated_vector)) * sum( (true_vector - estimated_vector)^2 ))))
}

#' Calculate normalized RMSE metric value.  From 0 to 1. Lower values are better.
#'
#'@param estimated_vector estimated row or column vector. Should be all positive.
#'@param true_vector row or column vector to compare with.  Should be all positive.
#'@return  a value
#'@export
normalized_rmse_loss_function <- function(estimated_vector, true_vector) {
    true_vector <- (true_vector - min(true_vector)) / (max(true_vector) - min(true_vector))
    estimated_vector <- (estimated_vector - min(estimated_vector)) / (max(estimated_vector) - min(estimated_vector))
    return(rmse_loss_function(estimated_vector, true_vector))
}















