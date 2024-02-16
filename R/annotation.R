# Simple function for manual use
annotate_data <- function(data, distances_svd_num) {
  eset <- create_eset(data)
  eset <- add_default_anno(eset)
  scaling <- sinkhorn_scale(exprs(eset))
  proj_full <- svd_project(scaling, dims = NULL)
  eset <- add_distances_anno(eset, proj_full, distances_svd_num)
  return(eset)
}

############ MAIN LOGIC ############
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


# Basic annotation logic
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


# Distances annotation

# TODO: this function shamelessly peeks into further stages than
# annotation, but its required for filtering
add_distances_anno <- function(eset, proj_full, n_cell_types) {
  fData(eset)$plane_distance <- calc_partial_dist(
    proj_full$X,
    with_dims = -c(1, 2:n_cell_types)
  )

  pData(eset)$plane_distance <- calc_partial_dist(
    proj_full$Omega,
    with_dims = -c(1, 2:n_cell_types)
  )

  fData(eset)$zero_distance <- calc_partial_dist(
    proj_full$X,
    with_dims = 2:n_cell_types
  )

  pData(eset)$zero_distance <- calc_partial_dist(
    proj_full$Omega,
    with_dims = 2:n_cell_types
  )

  return(eset)
}

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

describe_cat_features <- function(eset, genes = T) {
  anno <- get_anno(eset, genes)
  cat_cols <- colnames(anno)[sapply(colnames(anno), function(x) is.factor(anno[, x]) | is.logical(anno[, x]))]
  cat_features_anno <- anno[, cat_cols, drop = F]
  return(Hmisc::describe(cat_features_anno))
}
