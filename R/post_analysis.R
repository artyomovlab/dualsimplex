############ SOLUTION STATS ############
plot_proportions_distribution <- function(H) {
  to_plot <- as.data.frame(t(H))
  to_plot <- tidyr::pivot_longer(
    to_plot,
    colnames(to_plot),
    names_to = "cell_type",
    values_to = "proportion"
  )
  plt <- ggplot(
    to_plot,
    aes(x = proportion, fill = cell_type)
  ) + geom_histogram(
    stat = "bin",
    binwidth = 0.05,
    boundary = 0
  ) + xlim(0, 1) + ylab("samples")
  return(plt)
}

plot_basis_distribution <- function(W, max_expr = 25, logp1 = T) {
  if (logp1) {
    to_plot <- log(W + 1)
  } else {
    to_plot <- W
  }
  to_plot <- as.data.frame(to_plot)
  to_plot <- tidyr::pivot_longer(
    to_plot,
    colnames(to_plot),
    names_to = "cell_type",
    values_to = "expression"
  )
  plt <- ggplot(
    to_plot,
    aes(x = expression, fill = cell_type)
  ) + geom_histogram(
    stat = "bin",
    binwidth = 0.5,
    boundary = 0
  ) + xlim(0, max(max_expr, max(to_plot$expression))) + ylab("genes")
  return(plt)
}


############ MARKERS ############
get_signature_markers <- function(signature, n_marker_genes = 100, stat = "mean_fc") {
  stat <- get_fold_change(signature, stat)
  markers <- lapply(colnames(stat), function(ct) {
    as.data.frame(stat) %>% slice_max(n = n_marker_genes, order_by = !!sym(ct), with_ties = F) %>% rownames()
  })
  names(markers) <- colnames(stat)
  return(markers)
}

revert_marker_list <- function(marker_list) {
  rev_mlist <- list()
  for (ct in names(marker_list)) {
    for (gene in marker_list[[ct]])
      rev_mlist[[gene]] <- ct
  }
  return(rev_mlist)
}

which_marker <- function(gene_names, marker_list) {
  rev_mlist <- revert_marker_list(marker_list)
  which_marker_anno <- rep(NA, length(gene_names))
  names(which_marker_anno) <- gene_names
  which_marker_anno[names(rev_mlist)] <- rev_mlist
  which_marker_anno <- factor(unlist(which_marker_anno[gene_names]), levels = c(names(marker_list), NA))
  return(which_marker_anno)
}

get_fold_change <- function(signature, stat = "mean_fc") {
  cell_types_stats <- matrix(0, ncol = ncol(signature), nrow = nrow(signature))
  rownames(cell_types_stats) <- rownames(signature)
  colnames(cell_types_stats) <- colnames(signature)
  for (ct_col_i in 1:ncol(cell_types_stats)) {
    stat_fun <- if (stat == "mean_fc") 
      function(gene_row) {
        gene_row[ct_col_i] / mean(gene_row[-ct_col_i])
      }
    else if (stat == "min_fc")
      function(gene_row) {
        gene_row[ct_col_i] / min(gene_row[-ct_col_i] + 1)
      }
    else
      stop(paste0("Unknown stat: ", stat))
    cell_types_stats[, ct_col_i] <- apply(signature, 1, stat_fun)
  }
  cell_types_stats
}

get_limma_fold_change <- function(signature){
  limma_fc <- matrix(0,nrow=nrow(signature),ncol=self$cell_types)
  colnames(limma_fc) <- paste0("LimmaFC_", 1:self$cell_types)
  for (ct in 1:self$cell_types) {
    expr <- signature[,grepl("^Cell.*",colnames(signature))]
    condition <- as.integer(paste0("Cell_type_",ct) != colnames(expr))
    design <- model.matrix(~ condition)
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    stats <- topTable(fit, number = nrow(signature))
    limma_fc[,ct] <- -stats[rownames(signature),"logFC"]
  }
  return(cbind(limma_fc, signature))
}

cat_markers <- function(marker_list) {
  for (ct in names(marker_list)) {
    cat(ct)
    cat("\n")
    cat(marker_list[[ct]])
    cat("\n\n")
  }
}

############ SINGLE CELL ############
add_list_markers <- function(so, markers, assay = "RNA") {
  for (ct in names(markers)) {
    if (any(markers[[ct]] %in% rownames(GetAssayData(so, assay = assay)))) {
      message(paste("Adding", ct))
      so <- AddModuleScore(so, features = list(markers[[ct]]), name = ct, assay = assay)
    } else {
      message(paste("Skipping", ct))
    }
  }
  so
}

plot_marker_enrichment <- function(
  so, marker_names, ncol = 4, limits = NULL, ggadd = function(plt) plt, wrap = T
) {
  p <- Seurat::FeaturePlot(
    so,
    features = paste0(marker_names, 1),
    #  label = T,
    #  label.size = 3,
    order = T,
    combine = FALSE,
    reduction = "umap"
  )
  
  for(i in 1:length(p)) {
    ct <- marker_names[[i]]
    p[[i]] <- ggadd(
      p[[i]] +
        scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")), limits = limits, oob = scales::squish) +
        ggtitle(ct)
    )
    #NoLegend()
  }
  
  if (wrap) {
    cowplot::plot_grid(plotlist = p, ncol = ncol)
  } else {
    p
  }
}
