#' Simple default filters for genes and mad
#'
#' This is not used anymore but can be used as external method
#'
#' @param eset Expression set (annotated matrix)
#' @param keep_n_genes how many top mad genes to leave
#' @return modified set
#' @export
apply_default_filters <- function(eset, keep_n_genes) {
  eset <- predicate_filter(eset, "!RPLS & !LOC & !ORF", genes = T)
  eset <- top_filter(eset, feature = "mad", top_n = keep_n_genes, genes = T)
  return(eset)
}

#' Simple filter by row name
#'
#' will remove annotation row name not present in `keep_names`
#'
#' @param eset Expression set (annotated matrix)
#' @param keep_names which genes to preserve
#' @param genes if TRUE filter rows otherwise columns of eset
#' @return modified set
#' @export
names_filter <- function(eset, keep_names, genes = T) {
  anno <- get_anno(eset, genes)
  anno_flt <- anno[rownames(anno) %in% keep_names, ]
  return(set_anno(anno_flt, eset, genes))
}

#' Simple filter by predicate applied to annotation
#'
#' will remove genes rows/columns not passing the predicate function
#'
#' @param eset Expression set (annotated matrix)
#' @param predicate predicate to test for annotation. (should be True/False function)
#' @param genes if TRUE filter rows otherwise columns of eset
#' @return modified set
#' @export
predicate_filter <- function(eset, predicate, genes = T) {
  anno <- get_anno(eset, genes)
  anno_flt <- filter_str(anno, predicate)
  return(set_anno(anno_flt, eset, genes))
}

#' Simple filter based on annotation values
#'
#' will remove or keep only existing values having TRUE
#'
#' @param eset Expression set (annotated matrix)
#' @param by_cols which columns of annotation to test
#' @param genes if TRUE filter rows otherwise columns of eset
#' @param remove_true should we remove or keep these columns
#' @return modified set
#' @export
bool_filter <- function(eset, by_cols, genes = T, remove_true = T) {
  sign <- if (remove_true) "!" else ""
  if (!is.null(by_cols)) {
      eset <- predicate_filter(
        eset,
        paste0(
          sign,
          paste(by_cols, collapse = paste0(" & ", sign))
        ),
        genes
      )
  }
  return(eset)
}

#' Simple threshold filter based on annotation values
#'
#' will remove or keep only values below the threshold
#'
#' @param eset Expression set (annotated matrix)
#' @param feature numerical feature to test
#' @param threshold numerical threshold value
#' @param genes if TRUE filter rows else columns
#' @param keep_lower should we remove or keep values below threshold
#' @return modified set
#' @export
threshold_filter <- function(
  eset,
  feature,
  threshold,
  genes = T,
  keep_lower = T
) {
  anno <- get_anno(eset, genes)
  if (keep_lower) {
    anno_flt <- dplyr::filter(anno, get(feature) < threshold)
  } else {
    anno_flt <- dplyr::filter(anno, get(feature) > threshold)
  }
  return(set_anno(anno_flt, eset, genes))
}


#' Simple top_filter filter based on annotation values
#'
#' will remove or keep only top N max/min values based on params
#'
#' @param eset Expression set (annotated matrix)
#' @param feature numerical feature to test
#' @param top_n how many top values to keep
#' @param genes if TRUE filter rows else columns
#' @param max if TRUE will return top highest, else top lowest
#' @return modified set
#' @export
top_filter <- function(eset, feature, top_n, genes = T, max = T) {
  anno <- get_anno(eset, genes)
  slice_fun <- if (max) dplyr::slice_max else dplyr::slice_min
  anno_flt <- slice_fun(anno, order_by = .data[[feature]], n = top_n)
  anno_flt <- anno[rownames(anno) %in% rownames(anno_flt), ]
  return(set_anno(anno_flt, eset, genes))
}

#' Simple quantile filter based on annotation values
#'
#' will remove or keep only top N falling in specific quantile based on params
#'
#' @param eset Expression set (annotated matrix)
#' @param feature numerical feature to test
#' @param quant which quantile
#' @param genes if TRUE filter rows else columns
#' @param keep_lower if TRUE will return top highest, else top lowest
#' @return modified set
#' @export
quantile_filter <- function(eset, feature, quant, genes = T, keep_lower = T) {
  feature_col <- get_anno(eset, genes, feature)
  threshold <- stats::quantile(feature_col, quant)
  return(threshold_filter(eset, feature, threshold, genes, keep_lower))
}


#' Simple n_sigma filter based on annotation values
#'
#' will keep only values lying in a n_sigma interval for the feature
#'
#' @param eset Expression set (annotated matrix)
#' @param feature numerical feature to test
#' @param n_sigma how many standard deviations to keep
#' @param genes if TRUE filter rows else columns
#' @return modified set
#' @export
n_sigma_filter <- function(eset, feature, n_sigma = 3, genes = T) {
  feature_col <- get_anno(eset, genes, feature)
  feature_col <-  abs(feature_col - mean(feature_col))
  threshold <-  n_sigma * stats::sd(feature_col)
  return(threshold_filter(eset, feature, threshold, genes, keep_lower = T))
}
