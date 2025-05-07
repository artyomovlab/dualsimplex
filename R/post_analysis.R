############ SOLUTION STATS ############

#' Plot distribution of values for given matix H
#'
#' Plot distribution of values for given matix H. Color by rows of the matrix
#'
#' @param H coefficients matrix to be plotted
#' @return ggplot object
#' @export
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
    aes(x = .data$proportion, fill = .data$cell_type)
  ) + geom_histogram(
    stat = "bin",
    binwidth = 0.05,
    boundary = 0
  ) + xlim(0, 1) + ylab("samples")
  return(plt)
}

#' Plot distribution of values for given matix W
#'
#' Plot distribution of values for given matix W. Color by columns of the matrix
#'
#' @param W coefficients matrix to be plotted
#' @param max_expr maximum value of exprassion to be set (higher values will set to this)
#' @param logp1 TRUE if want to log scale basis matrix
#' @return ggplot object
#' @export
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
    aes(x = expression, fill = .data$cell_type)
  ) + geom_histogram(
    stat = "bin",
    binwidth = 0.5,
    boundary = 0
  ) + xlim(0, max(max_expr, max(to_plot$expression))) + ylab("genes")
  return(plt)
}


############ MARKERS ############

#' Get signature markers based on basis matrix.
#'
#' Will calculate fold change for rows and return top n rownames
#'
#' @param signature basis matrix (W) obtained by method or any other
#' @param n_marker_genes number of genes to return for each column
#' @param stat which statistic to calculate for fold change
#' @return list of lists of genes (for each column of the matrix `signature`)
#' @export
get_signature_markers <- function(signature, n_marker_genes = 100, stat = "mean_fc") {
  stat <- get_fold_change(signature, stat)
  markers <- lapply(colnames(stat), function(ct) {
    as.data.frame(stat) %>% slice_max(n = n_marker_genes, order_by = !!sym(ct), with_ties = F) %>% rownames()
  })
  names(markers) <- colnames(stat)
  return(markers)
}

#' Convert cell_type:gene list to gene:cell_type list
#'
#' Utility method
#'
#' @param marker_list list of lists of genes.
#' @return same list but now gene is key, cell type is value
#' @export
revert_marker_list <- function(marker_list) {
  rev_mlist <- list()
  for (ct in names(marker_list)) {
    for (gene in marker_list[[ct]])
      rev_mlist[[gene]] <- ct
  }
  return(rev_mlist)
}

#' Assign cell types to marker genes
#'
#' Utility method
#'
#' @param gene_names gene list to assing cell types
#' @param marker_list list of lists of marker genes for each cell type
#' @return return marker label for gene if it is marker
#' @export
which_marker <- function(gene_names, marker_list) {
  rev_mlist <- revert_marker_list(marker_list)
  which_marker_anno <- rep(NA, length(gene_names))
  names(which_marker_anno) <- gene_names
  which_marker_anno[names(rev_mlist)] <- rev_mlist
  which_marker_anno <- factor(unlist(which_marker_anno[gene_names]), levels = c(names(marker_list), NA))
  return(which_marker_anno)
}

#' Get fold change values for genes
#'
#' Utility method
#'
#' @param signature basis matrix (W)
#' @param stat statistic to calculate
#' @param colnorm normalize result value or not
#' @return cell_types_stats (fold change values for rows)
#' @export
get_fold_change <- function(signature, stat = "mean_fc", colnorm = T) {
  cell_types_stats <- matrix(0, ncol = ncol(signature), nrow = nrow(signature))
  rownames(cell_types_stats) <- rownames(signature)
  colnames(cell_types_stats) <- colnames(signature)
  if (colnorm) {
    signature <- apply(signature, 2, function(x) x / sum(x))
  }
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

#' Convert marker list to string
#'
#' Nice string of markers for matrix W
#'
#' @param marker_list markers to convert (list of lists)
#' @return string with all markers listed
#' @export
cat_markers <- function(marker_list) {
  for (ct in names(marker_list)) {
    cat(ct)
    cat("\n")
    cat(marker_list[[ct]])
    cat("\n\n")
  }
}

############ SINGLE CELL ############

#' Add markers to Seurat object
#'
#' Test if marker is in object using (GetAssayData), add markers using (AddModuleScore)
#'
#' @param so Seurat object
#' @param markers markers list of list
#' @param assay assay name
#' @return Seurat object
#' @export
add_list_markers <- function(so, markers, assay = "RNA") {
  for (ct in names(markers)) {
    if (any(markers[[ct]] %in% rownames(Seurat::GetAssayData(so, assay = assay)))) {
      message(paste("Adding", ct))
      so <- Seurat::AddModuleScore(so, features = list(markers[[ct]]), name = ct, assay = assay)
    } else {
      message(paste("Skipping", ct))
    }
  }
  so
}

#' Transform seurat marker genes to our genes
#'
#' Basically just split df by cell type
#'
#' @param markers_seurat marker genes extracted from seurat
#' @param allowed_genes which genes to preserve
#' @param n_markers how many markers we need
#' @return list of lists
#' @export
convert_sc_markers <- function(markers_seurat, allowed_genes, n_markers = 100) {
  gb <- markers_seurat %>% filter(.data$p_val_adj < 0.001) %>%
    group_by(.data$cluster) %>%
    filter(.data$gene %in% allowed_genes) %>% slice_max(order_by=.data$avg_log2FC, n=n_markers)
  sc_marker_list <- sapply(gb %>% group_split(), function(x) list(x$gene))
  names(sc_marker_list) <- tolower(gsub(" |-", "_", (gb %>% group_keys())$cluster))
  return(sc_marker_list)
}

#' Plot enrichement for markers
#'
#' Using Seurat::FeaturePlot to plot
#'
#' @param so Seurat object
#' @param marker_names marker genes for each cell type (list of lists)
#' @param ncol number of columns in result plot (cowplot::plot_grid)
#' @param limits color limits (to adjust color)
#' @param ggadd add text to each plot (function taking plot and index of plot)
#' @param wrap whould we do grid or not
#' @param ... any valid Seurat::FeaturePlot params
#' @return single plot or multiple plots
#' @import RColorBrewer
#' @import cowplot
#' @import scales
#' @export
plot_marker_enrichment <- function(
  so, marker_names, ncol = 4, limits = NULL, ggadd = function(plt, i) plt, wrap = T, ...
) {
  p <- Seurat::FeaturePlot(
    so,
    features = paste0(marker_names, 1),
    #  label = T,
    #  label.size = 3,
    order = T,
    combine = FALSE,
    ...
  )

  for(i in 1:length(p)) {
    ct <- marker_names[[i]]
    p[[i]] <- ggadd(
      p[[i]] +
        scale_colour_gradientn(
          colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")),
          limits = if (is.list(limits)) limits[[i]] else limits,
          oob = scales::squish
        ) +
        ggtitle(ct),
      i
    )
    #NoLegend()
  }

  if (wrap) {
    cowplot::plot_grid(plotlist = p, ncol = ncol)
  } else {
    p
  }
}
