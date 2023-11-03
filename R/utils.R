calc_partial_dist <- function(data_proj_single, with_dims) {
  data_part <- data_proj_single
  data_part[, -with_dims] <- 0
  return(euclidean_dist(data_part))
}

euclidean_dist <- function(data) {
  return(sqrt(rowSums(data^2)))
}

linearize_dataset = function(ge) {
  if (is_logscale(ge))
    return(2^ge - 1)
  return(ge)
}

log_dataset = function(ge) {
  if (is_logscale(ge))
    return(ge)
  return(log2(ge + 1))
}

rownorm <- function(ge) {
  return(ge / rowSums(ge))
}

is_logscale = function(x) {
  qx <- quantile(as.numeric(x), na.rm = T)
  if (qx[5] - qx[1] > 100 || qx[5] > 100) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

remove_zero_rows <- function(ge) {
  zero_rows <- apply(data_raw, 1, function(x) all(x == 0))
  return(ge[!zero_rows, ])
}

filter_str <- function(data, expr, ...) {
  return(dplyr::filter(data, eval(rlang::parse_expr(expr)), ...))
}

read_gene_list <- function(path) {
  genes_string <- readChar(path, file.info(path)$size)
  genes_list <- stringr::str_split(genes_string, " ")[[1]]
  return(genes_list)
}

rasterize_if_needed <- function(plot) {
  plot <- if(getOption("linseed-rasterize", default=FALSE)) ggrastr::rasterise(plot, dpi=600) else plot
  return(plot)
}
  
  