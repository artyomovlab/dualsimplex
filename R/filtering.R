#' Simple default filters for genes and mad
#'
#' This is not used anymore but can be used as external method
#'
#' @param eset Expression set (annotated matrix)
#' @param keep_n_genes how many top mad genes to leave
#' @return modified set
#' @export
apply_default_filters <- function(eset, keep_n_genes) {
  eset <- predicate_filter(eset, "!RPLS & !LOC & !ORF", for_features = T)
  eset <- top_filter(eset, annotation_feature = "mad", top_n = keep_n_genes, for_features = T)
  return(eset)
}

#' Simple filter by row name
#'
#' will remove annotation row name not present in `keep_names`
#'
#' @param eset Expression set (annotated matrix)
#' @param keep_names which genes to preserve
#' @param for_features if TRUE filter rows otherwise columns of eset
#' @return modified set
#' @export
names_filter <- function(eset, keep_names, for_features = T) {
  anno <- get_anno(eset, for_features)
  anno_flt <- anno[rownames(anno) %in% keep_names, ]
  return(set_anno(anno_flt, eset, for_features))
}

#' Simple filter by predicate applied to annotation
#'
#' will remove rows/columns not passing the predicate function
#'
#' @param eset Expression set (annotated matrix)
#' @param predicate predicate to test for annotation. (should be True/False function)
#' @param for_features if TRUE filter rows otherwise columns of eset
#' @return modified set
#' @export
predicate_filter <- function(eset, predicate, for_features = T) {
  anno <- get_anno(eset, for_features)
  anno_flt <- filter_str(anno, predicate)
  return(set_anno(anno_flt, eset, for_features))
}

#' Simple filter based on annotation values
#'
#' will remove or keep only existing values having TRUE
#'
#' @param eset Expression set (annotated matrix)
#' @param by_cols which columns of annotation to test
#' @param for_features if TRUE filter rows otherwise columns of eset
#' @param remove_true should we remove or keep these columns
#' @return modified set
#' @export
bool_filter <- function(eset, by_cols, for_features = T, remove_true = T) {
  sign <- if (remove_true) "!" else ""
  if (!is.null(by_cols)) {
      eset <- predicate_filter(
        eset,
        paste0(
          sign,
          paste(by_cols, collapse = paste0(" & ", sign))
        ),
        for_features
      )
  }
  return(eset)
}

#' Simple threshold filter based on annotation values
#'
#' will remove or keep only values below the threshold
#'
#' @param eset Expression set (annotated matrix)
#' @param annotation_feature numerical feature to test
#' @param threshold numerical threshold value
#' @param for_features if TRUE filter rows else columns
#' @param keep_lower should we remove or keep values below threshold
#' @return modified set
#' @export
threshold_filter <- function(
  eset,
  annotation_feature,
  threshold,
  for_features = TRUE,
  keep_lower = TRUE
) {
  anno <- get_anno(eset, for_features)
  if (keep_lower) {
    anno_flt <- dplyr::filter(anno, get(annotation_feature) < threshold)
  } else {
    anno_flt <- dplyr::filter(anno, get(annotation_feature) >= threshold)
  }
  return(set_anno(anno_flt, eset, for_features))
}

#' Simple range filter based on annotation values
#'
#' will remove or keep only values in specified range
#'
#' @param eset Expression set (annotated matrix)
#' @param annotation_feature numerical feature to test
#' @param threshold_lower numerical threshold value
#' @param threshold_upper numerical threshold value
#' @param for_features if TRUE filter rows else columns
#' @param keep_within should we keep or remove values within the range
#' @return modified set
#' @export
range_filter <- function(
  eset,
  annotation_feature,
  threshold_lower,
  threshold_upper,
  for_features = T,
  keep_within = T
) {
  anno <- get_anno(eset, for_features)
  if (keep_within) {
    anno_flt <- dplyr::filter(anno, between(get(annotation_feature),threshold_lower, threshold_upper))
  } else {
    anno_flt <- dplyr::filter(anno, !between(get(annotation_feature),threshold_lower, threshold_upper))
  }
  return(set_anno(anno_flt, eset, for_features))
}




#' Simple top_filter filter based on annotation values
#'
#' will remove or keep only top N max/min values based on params
#'
#' @param eset Expression set (annotated matrix)
#' @param annotation_feature numerical feature to test
#' @param top_n how many top values to keep
#' @param for_features if TRUE filter rows else columns
#' @param max if TRUE will return top highest, else top lowest
#' @return modified set
#' @export
top_filter <- function(eset, annotation_feature, top_n, for_features = T, max = T) {
  anno <- get_anno(eset, for_features)
  slice_fun <- if (max) dplyr::slice_max else dplyr::slice_min
  anno_flt <- slice_fun(anno, order_by = .data[[annotation_feature]], n = top_n)
  anno_flt <- anno[rownames(anno) %in% rownames(anno_flt), ]
  return(set_anno(anno_flt, eset, for_features))
}


#' Simple quantile filter based on annotation values
#'
#' will remove or keep only top N falling in specific quantile based on params
#'
#' @param eset Expression set (annotated matrix)
#' @param annotation_feature numerical feature to test
#' @param quant which quantile
#' @param for_features if TRUE filter rows else columns
#' @param keep_lower if TRUE will return top highest, else top lowest
#' @return modified set
#' @export
quantile_filter <- function(eset, annotation_feature, quant, for_features = T, keep_lower = T) {
  feature_col <- get_anno(eset, for_features, annotation_feature)
  threshold <- stats::quantile(feature_col, quant)
  return(threshold_filter(eset, annotation_feature, threshold, for_features, keep_lower))
}


#' Simple n_sigma filter based on annotation values
#'
#' will keep only values lying in a n_sigma interval for the feature
#'
#' @param eset Expression set (annotated matrix)
#' @param annotation_feature numerical feature to test
#' @param n_sigma how many standard deviations to keep
#' @param for_features if TRUE filter rows else columns
#' @return modified set
#' @export
n_sigma_filter <- function(eset, annotation_feature, n_sigma = 3, for_features = T) {
  feature_col <- get_anno(eset, for_features, annotation_feature)
  sigma <-  stats::sd(feature_col)
  lower_bound <- mean(feature_col) - n_sigma * sigma
  upper_bound <- mean(feature_col) + n_sigma * sigma
  return(range_filter(eset, annotation_feature, threshold_lower = lower_bound,  threshold_upper = upper_bound, for_features, keep_within = T))
}


#' Simple n_sigma filter based on annotation values
#'
#' will keep only values lying in a n_sigma interval for the feature
#'
#' @param eset Expression set (annotated matrix)
#' @param annotation_features numerical features to test
#' @param n_sigma how many standard deviations to keep
#' @param for_features if TRUE filter rows else columns
#' @param keep_lower wether to take values inside interval
#' @return modified set
#' @export
mahalanobis_n_sigma_filter <- function(eset, annotation_features, n_sigma = 3, for_features = T, keep_lower=T) {
  feature_cols <- get_anno(eset, for_features, annotation_features)
  Sx <- stats::cov(feature_cols)
  distance_values <- sqrt(stats::mahalanobis(feature_cols, colMeans(feature_cols), Sx))
  anno <- get_anno(eset, for_features)
  if (keep_lower) {
    anno_flt <- anno[distance_values < n_sigma,]
  } else {
    anno_flt <- anno[distance_values >= n_sigma,]
  }
  return(set_anno(anno_flt, eset, for_features))
}
