#' Add dannotations for rows (genes) and columns (samples) of the matrix
#'
#' Will add (log_mean, log_median, log_sd, log_mad) for both rows and columns.
#' For rows additionaly will produce columns for TRUE/FALSE regex filters specified in `add_gene_names_regex_anno`
#' and also add TRUE/FALSE columns specified in gene_name_lists
#' For columns additionaly will produce TRUE/FALSE columns pecified in sample_name_lists
#'
#' @param eset Expression set
#' @param gene_name_lists named list of lists. Each sublist contains names of rows which should have TRUE value in annotaiton column.
#' @param sample_name_lists named list of lists. Each sublist contains names of columns which should have TRUE value in annotation column.
#' @return annotated expression set. (fData, pData) now contain annotations
#' @export
add_default_anno <- function(
  eset,
  gene_name_lists = NULL,
  sample_name_lists = NULL
) {
  eset <- add_default_gene_anno(eset, gene_name_lists)
  eset <- add_default_sample_anno(eset, sample_name_lists)
  return(eset)
}

add_default_gene_anno <- function(eset, name_lists = NULL) {
  eset <- add_data_stats_anno(eset, genes = T)
  eset <- add_gene_names_regex_anno(eset)
  if (!is.null(name_lists)) {
    eset <- add_name_lists_anno(eset, name_lists, genes = T)
  }
  return(eset)
}

add_default_sample_anno <- function(eset, name_lists = NULL) {
  eset <- add_data_stats_anno(eset, genes = F)
  if (!is.null(name_lists)) {
    eset <- add_name_lists_anno(eset, name_lists, genes = F)
  }
  return(eset)
}

# eset object initialization:
create_eset <- function(data) {
  ExpressionSet(
    assayData = data,
    featureData = create_gene_anno(data),
    phenoData = create_sample_anno(data)
  )
}

create_sample_anno <- function(data) {
  sample_anno <- data.frame(matrix(nrow=ncol(data), ncol=0))
  rownames(sample_anno) <- colnames(data)
  return(AnnotatedDataFrame(sample_anno))
}

create_gene_anno <- function(data) {
  gene_anno <- data.frame(matrix(nrow=nrow(data), ncol=0))
  rownames(gene_anno) <- rownames(data)
  return(AnnotatedDataFrame(gene_anno))
}


#' Annotate gene names in rows according to selected regex filters
#'
#' Create columns with TRUE/FALSE values wether gene names pass basic regex filters
#'@param eset Expression set
#'@return  annotated Expression set
add_gene_names_regex_anno <- function(eset) {
  fData(eset)$RPLS <- grepl("^(RPL|RPS).+", rownames(exprs(eset)), ignore.case = T)
  fData(eset)$LOC <- grepl("^LOC\\d+", rownames(exprs(eset)), ignore.case = T)
  fData(eset)$ORF <- grepl("^C\\w+orf\\d+", rownames(exprs(eset)), ignore.case = T)
  fData(eset)$SNOR <- grepl("^SNOR.+", rownames(exprs(eset)), ignore.case = T)
  fData(eset)$MT <- grepl("^MT-.+", rownames(exprs(eset)), ignore.case = T)
  fData(eset)$RRNA <- grepl(".+rRNA$", rownames(exprs(eset)), ignore.case = T)
  return(eset)
}

add_name_lists_anno <- function(eset, name_lists, genes = T) {
  anno <- get_anno(eset, genes)
  for (anno_name in names(name_lists)) {
    anno[, anno_name] <- rownames(anno) %in% name_lists[[anno_name]]
  }
  eset <- set_anno(anno, eset, genes)
  return(eset)
}

#' Annotate the data with statistics for rows and columns
#'
#' Will add (log_mean, log_median, log_sd, log_mad) for both rows and columns.
#'
#'@param eset Expression set
#'@param genes if TRUE apply to rows, otherwise columns
#'@return  annotated Expression set
add_data_stats_anno <- function(eset, genes = T) {
  anno <- get_anno(eset, genes)
  margin <- if (genes) 1 else 2
  anno$log_mean <- log(apply(exprs(eset), margin, mean) + 1)
  anno$log_median <- log(apply(exprs(eset), margin, median) + 1)
  anno$log_sd <- log(apply(exprs(eset), margin, sd) + 1)
  anno$log_mad <- log(apply(exprs(eset), margin, mad) + 1)
  eset <- set_anno(anno, eset, genes)
  return(eset)
}


#' Distances annotation
#' 
#' Add distance information in the expression set object
#' 
#' @param eset expression set to be modified
#' @param V_row Sinkhorn scaled matrix
#' @param proj a projectino object
#' 
#' @return an expression set object with distance information
#' 
#' TODO (cjlee): make it more efficient. by Rcpp maybe?
add_distances_anno <- function(eset, V_row, proj) {
  approx <- with(proj@meta, t(S) %*% Sigma %*% R)
  residual <- V_row - approx

  # For features (using V_row)
  feature_dists <- calc_dist_from_truncated_svd(approx, residual = residual, margin = 1)
  fData(eset)$plane_distance <- feature_dists$plane_distance
  fData(eset)$zero_distance <- feature_dists$zero_distance


  # For samples (using V_column)
  d <- ncol(V_row) / nrow(V_row)
  sample_dists <- calc_dist_from_truncated_svd(approx * d, residual = residual * d, margin = 2)
  pData(eset)$plane_distance <- sample_dists$plane_distance
  pData(eset)$zero_distance <- sample_dists$zero_distance

  return(eset)
}

add_knn_distances_anno <- function(eset, proj_full, annotation_columns,  k_neighbors, genes=T) {
  anno <- get_anno(eset, genes)
  for (anno_name in annotation_columns) {
    knns <-  FNN::get.knnx(proj_full$X[anno[[anno_name]], ], proj_full$X, k = k_neighbors)
    distances <- apply(knns$nn.dist, 1, min, na.rm=T)
    anno[, paste0(anno_name, "_subset_distance")] <- distances
  }
  eset <- set_anno(anno, eset, genes)
  return(eset)
}

#' Get annotations for the data
#'
#' Returns fData and pData of the Expresison set.
#'
#' @param eset Expression set
#' @param genes if TRUE return row annotations (genes) if FALSE return column annotaitons (samples)
#' @param feature name of specific annotation you want to extract
#' @return annotation object of the data
#' @export
get_anno <- function(eset, genes = T, feature = NULL) {
  if (genes) {
    anno <- fData(eset)
  } else {
    anno <- pData(eset)
  }
  if (!is.null(feature)) {
    anno <- anno[, feature]
  }
  return(anno)
}

#' Set custom annotations for the data
#'
#' sets fData and pData of the Expresison set using provided data
#'
#' @param anno annotation object to be set for the data
#' @param eset Expression set
#' @param genes if TRUE set row annotations (genes) if FALSE set column annotaitons (samples)
#' @export
set_anno <- function(anno, eset, genes = T) {
  if (genes) {
    eset <- eset[rownames(anno), ]
    fData(eset) <- anno
  } else {
    eset <- eset[, rownames(anno)]
    pData(eset) <- anno
  }
  return(eset)
}


############ PLOTTING ############

#' Plot distribution of specific feature from annotation
#'
#' just geom_histogram for the feature selected
#'
#' @param eset Expression set
#' @param feature annotation variable name (e.g. log_mad)
#' @param genes if TRUE plot for plot for rows, otherwise columns
#' @param col_by name of the feature to color by
#' @param bins bin number for a historgramm
#' @return ggplot object
#' @export
plot_feature <- function(eset, feature, genes = T, col_by = NULL, bins = 100) {
  anno <- get_anno(eset, genes)
  fill <- if (is.null(col_by)) "grey40" else "white"
  plt <- if (is.numeric(anno[, feature])) {
    ggplot(anno, aes_string(x = feature, color = col_by)) +
      geom_histogram(bins = bins, fill = fill) +
      theme_minimal() +
      ggtitle(paste0(nrow(anno), if (genes) " genes" else " samples"))
  } else {
    stop("Non-numeric features are not supported")
  }
  return(plt)
}

#' Plot scatterplot for a pair of annotation features
#'
#' Uses ggplot geom_point() you can add any valid geom_point(...) parameters after all arguments
#'
#' @param eset Expression set
#' @param feature_1 annotation variable name (e.g. log_mad)
#' @param feature_2 annotation variable name (e.g. log_mad)
#' @param genes if TRUE plot for plot for rows, otherwise columns
#' @param col_by name of the feature to color by
#' @param ... any valid geom_point(...) parameters
#' @return ggplot object
#' @export
plot_feature_pair <- function(
  eset,
  feature_1,
  feature_2,
  genes = T,
  col_by = NULL,
  ...
) {
  anno <- get_anno(eset, genes)
  fill <- if (is.null(col_by)) "grey40" else "white"
  plt <- if (is.numeric(anno[, feature_1]) && is.numeric(anno[, feature_2])) {
    ggplot(anno, aes_string(x = feature_1, y = feature_2, color = col_by)) +
      geom_point(...) +
      theme_minimal() +
      ggtitle(paste0(nrow(anno), if (genes) " genes" else " samples"))
  } else {
    stop("Non-numeric features are not supported")
  }
  return(plt)
}

#' Plot distribution for multiple features from annotation
#'
#' geom_histogramm for selected features. Uses cowplot::plot_grid()
#'
#' @param eset Expression set.
#' @param genes if TRUE plot for plot for rows, otherwise columns.
#' @param col_by name of the feature to color by.
#' @param bins bin number for a historgramm.
#' @param features annotation variable name (e.g. log_mad).
#' @param ncol how many columns to make in cowplot.
#' @param labels additional label for points.
#' @return cowplot object.
#' @export
plot_numeric_features <- function(
  eset,
  genes = T,
  col_by = NULL,
  bins = 100,
  features = NULL,
  ncol = NULL,
  labels = "Numeric features"
) {
  anno <- get_anno(eset, genes)
  if (is.null(features)) {
    numeric_cols <- colnames(anno)[sapply(colnames(anno), function(x) is.numeric(anno[, x]))]
    features <- numeric_cols
  }
  plotlist <- lapply(features, function(feature) {
    plot_feature(eset, feature, genes, col_by, bins) + ggtitle("")
  })
  if (is.null(ncol)) {
    return(plotlist)
  } else {
    return(
      cowplot::plot_grid(
        plotlist = plotlist,
        ncol = ncol,
        labels = paste(if (genes) "Genes," else "Samples,", labels)
      )
    )
  }
}

#' Plot description of categorical features
#'
#' You need Hmisc package to run this.
#'
#' @param eset Expression set
#' @param genes if TRUE plot for plot for rows, otherwise columns
#' @return Hmisc result plot
#' @importFrom Hmisc describe
#' @export
describe_cat_features <- function(eset, genes = T) {
  anno <- get_anno(eset, genes)
  cat_cols <- colnames(anno)[sapply(colnames(anno), function(x) is.factor(anno[, x]) | is.logical(anno[, x]))]
  cat_features_anno <- anno[, cat_cols, drop = F]
  return(Hmisc::describe(cat_features_anno))
}
