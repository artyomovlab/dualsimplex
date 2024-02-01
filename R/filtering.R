apply_default_filters <- function(eset, keep_n_genes) {
  eset <- predicate_filter(eset, "!RPLS & !LOC & !ORF", genes = T)
  eset <- top_filter(eset, feature = "mad", top_n = keep_n_genes, genes = T)
  return(eset)
}

names_filter <- function(eset, keep_names, genes = T) {
  anno <- get_anno(eset, genes)
  anno_flt <- anno[rownames(anno) %in% keep_names, ]
  return(set_anno(anno_flt, eset, genes))
}

predicate_filter <- function(eset, predicate, genes = T) {
  anno <- get_anno(eset, genes)
  anno_flt <- filter_str(anno, predicate)
  return(set_anno(anno_flt, eset, genes))
}

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

top_filter <- function(eset, feature, top_n, genes = T, max = T) {
  anno <- get_anno(eset, genes)
  slice_fun <- if (max) dplyr::slice_max else dplyr::slice_min
  anno_flt <- slice_fun(anno, order_by = .data[[feature]], n = top_n)
  anno_flt <- anno[rownames(anno) %in% rownames(anno_flt), ]
  return(set_anno(anno_flt, eset, genes))
}

quantile_filter <- function(eset, feature, quant, genes = T, keep_lower = T) {
  feature_col <- get_anno(eset, genes, feature)
  threshold <- stats::quantile(feature_col, quant)
  return(threshold_filter(eset, feature, threshold, genes, keep_lower))
}

n_sigma_filter <- function(eset, feature, n_sigma = 3, genes = T) {
  feature_col <- get_anno(eset, genes, feature)
  threshold <- mean(feature_col) + n_sigma * sd(feature_col)
  return(threshold_filter(eset, feature, threshold, genes, keep_lower = T))
}
