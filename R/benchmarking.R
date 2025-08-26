############ MAIN LOGIC ############

#' Get similarity matrix for two matrices using specified metric and comparing rows.
#' If row has sd 0 it will add some random value to the row avoiding NaNs
#'
#'@param matrix_a first matrix
#'@param matrix_b second matrix
#'@param metric which metric to use
#'@param sd_fix_value standard deviation of the random values to add to the row in case it has sd of 0
#'@return result similarity matrix. Rows from the first matrix, columns from the second matrix.
#'@export
get_similarity_matrix_for_rows <- function(matrix_a, matrix_b, metric="rmse", sd_fix_value = 1e-4) {
  num_rows_a <-  dim(matrix_a)[[1]]
  num_rows_b <-  dim(matrix_b)[[1]]
  similarities <-  lapply(c(1:num_rows_a), function(index_a){
    vector_of_a <- matrix_a[index_a,]
    candidates <- lapply(c(1:num_rows_b), function(index_b){
      vector_of_b <- matrix_b[index_b,]
      if (metric %in% c('pearson', 'spearman')) {
        random_vector <-  abs(stats::rnorm(matrix_b[index_b,],0, sd=sd_fix_value))
        if (stats::sd(vector_of_a) == 0) {
          vector_of_a <- vector_of_a + random_vector
        }
        if (stats::sd(vector_of_b) == 0) {
          vector_of_b <- vector_of_b + random_vector
        }
      }
      result <- switch(metric,
                       "rmse" = 1 - rmse_loss_function(vector_of_a, vector_of_b),
                       "pearson" = pearson_correlation_function(vector_of_a, vector_of_b),
                       "spearman" = spearman_correlation_function(vector_of_a, vector_of_b),
                       "cosine" = cosine_similarity_function(vector_of_a, vector_of_b)
                       )
      return(result)
    })
    current_a_component_sim_row <- do.call(cbind, candidates)
    colnames(current_a_component_sim_row) <-  rownames(matrix_b)
    return(current_a_component_sim_row)
  })
  similarity_matrix <- do.call(rbind, similarities)
  rownames(similarity_matrix) <- rownames(matrix_a)
  colnames(similarity_matrix) <- rownames(matrix_b)
  return(abs(similarity_matrix))
}



#' Function to gess order of rows for the predicted markers to have the highest sum correlation with true matrix.
#' If the number of true matrix rows is smaller, will mark some rows as "extra rows" to keep predicted rows.
#' Implemented using lp.assign method.
#' Ensure your columns are in the matching order!
#' @param predicted_matrix predicted matrix
#' @param true_matrix true matrix
#' @return list containing row numbers for selected predicted cell type for given true cell type.
#' @importFrom lpSolve lp.assign
#' @export
guess_order <- function(predicted_matrix, true_matrix) {
  # similarity matrix (rows are predicted rows, columns are estimated columns)
  similarity_matrix <- get_similarity_matrix_for_rows(matrix_a = predicted_matrix, matrix_b = true_matrix, metric = "pearson")

# Add extra rows if we don't have enough true cell types
if (dim(similarity_matrix)[[1]] > dim(similarity_matrix)[[2]]) {
  need_to_add <-dim(similarity_matrix)[[1]] - dim(similarity_matrix)[[2]]
  separate_columns <- lapply(c(1:need_to_add), function(x) matrix(rep(0, dim(similarity_matrix)[[1]])))
  new_columns <-  do.call(cbind, separate_columns)
  colnames(new_columns) <-  paste0("extra_", c(1:need_to_add))
  similarity_matrix <- cbind(similarity_matrix, new_columns)
}
  asignment_result <- lpSolve::lp.assign(t(similarity_matrix), direction="max")
  new_order <-  apply(asignment_result$solution, MARGIN = 1, which.max)
  new_order <- as.list(new_order)
  names(new_order) <-  colnames(similarity_matrix)
  return(new_order)
}



#' Generate (prediction,true) pairs for each basis vector
#'
#' Will try to gess the right order of columns (since algorithm can shuffle it) for matrix W
#'
#' @param pred_basis  MxK predicted matrix W
#' @param true_basis MxK true matrix W
#' @return list of 2 matrices with ordered columns (1-1, 2-2, ..)
#' @export
coerce_pred_true_basis <- function(pred_basis, true_basis) {
  new_col_order <- guess_order(t(pred_basis), t(true_basis))[colnames(true_basis)] # Return only matching cell types
  pred_basis <- pred_basis[, unlist(new_col_order)]
  colnames(pred_basis) <- colnames(true_basis)
  rownames(pred_basis) <- rownames(true_basis)
  list(pred_basis, true_basis)
}

#' Generate (prediction,true) pairs for each coefficients vector
#'
#' Will try to gess the right order of rows (since algorithm can shuffle it) for matrix H
#'
#' @param pred_props  KxN predicted matrix H
#' @param true_props KxN true matrix H
#' @return list of 2 matrices with ordered rows (1-1, 2-2, ..)
#' @export
coerce_pred_true_props <- function(pred_props, true_props) {
  new_row_order <- guess_order(pred_props, true_props)[rownames(true_props)] # Return only matching cell types
  pred_props <- pred_props[unlist(new_row_order), ]
  colnames(pred_props) <- colnames(true_props)
  rownames(pred_props) <- rownames(true_props)
  list(pred_props, true_props)
}





toMatrix <- function(x) {
    if (is.data.frame(x)) {
        # Convert data frame (or tibble) to a plain matrix
        return(as.matrix(x))
    }
    if (is.matrix(x)) {
        # Return if already a matrix
        return(x)
    }
    stop("Invalid type for plotting: ", paste(class(x), collapse = ", "))
}


############ PLOTTING ############

#' Draw a plot of estimated proportions
#'
#' Draws a plot of estimated proprotions
#' If ggplot2 and reshape2 are installed will use them and return ggplot object
#' Otherwise will use standart R functions
#'
#' @param ... matricies, data frames, NMF objects of estimated proportions or paths to file
#' @param point_size point size for plot
#' @param line_size line size for plot
#' @param pnames experiment titles
#'
#' @return ggplot object
#'
#'
#' @import ggplot2
#' @import reshape2
#' @export
plotProportions <- function(..., pnames = NULL, point_size=2, line_size=1) {
    proportions <- list(...)
    proportions <- lapply(proportions, toMatrix)

    newRowNames <- do.call(function(...) {
        mapply(function(...) {
            dots <- list(...)
            rn <- unlist(dots)
            paste0(rn, collapse = "\n")
        }, ...)
    }, lapply(proportions, rownames))

    proportions <- lapply(proportions, function(p) {
        rownames(p) <- newRowNames
        p
    })

    names(proportions) <- pnames


    cellTypes <- nrow(proportions[[1]])
    results.m <- melt(proportions)
    results.m[, 4] <- as.factor(results.m[, 4])

    results.m <- results.m[sample(nrow(results.m)), ]

    gplot <- ggplot(results.m,
                             aes(x = as.numeric(.data$Var2),
                                          y = .data$value,
                                          fill = .data$Var1,
                                          color = .data$L1)) +
        geom_line(size=line_size) +
        geom_point(size=point_size) +
        scale_x_discrete(labels = colnames(proportions[[1]])) +
        facet_grid(Var1 ~ .) +
        ylab("proportions") +
        ylim(0, 1.1) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
                       axis.text.x = element_text(angle = 45,
                                                           hjust = 1)) +
        guides(fill = "none")
    if (length(proportions) > 1) {
        gplot <- gplot + theme(legend.title = element_blank(),
            legend.position = "top")

    } else {
        gplot <- gplot + theme(legend.position = "none")
    }
    gplot
}


#' Plot predicted/true lines plot for proportions/coefficients (matrix H)
#'
#' Using here method exported from linseed package.
#'
#' @param ptp result of coerce_pred_true_props or list of two matrices with the same order of rows
#' @return lineplot
#' @export
plot_ptp_lines <- function(ptp) {
  plotProportions(as.data.frame(ptp[[1]]), as.data.frame(ptp[[2]]), pnames = c("predicted", "true"))
}

#' Plot predicted/true scatter plot for proportions/coefficients (matrix H)
#'
#' Uses ggplot with rasterization to plot and  gridExtra::grid to arrange plots
#'
#' @param ptp result of coerce_pred_true_props or list of two matrices with the same order of rows
#' @return  gridExtra::grid.arrange() plot
#' @importFrom Metrics rmse
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob
#' @export
plot_ptp_scatter <- function(ptp) {
  plot_list <- list()
  cts <- rownames(ptp[[1]])
  for (ct in cts) {
    to_plot <- as.data.frame(list(
      DualSimplex = ptp[[1]][ct, ],
      other = ptp[[2]][ct, ]
    ))

    rmse_value <- round(Metrics::rmse(ptp[[1]][ct, ], ptp[[2]][ct, ]), 2)

    mx <- 1

    plot_list[[ct]] <- ggplot(to_plot, aes(x=.data$other, y=.data$DualSimplex)) +
      rasterize_if_needed(geom_point(size = 5, shape = 1, colour = "black")) +
      ggtitle(paste0(ct)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1.5) +
      xlab(NULL) +
      ylab(NULL) +
      xlim(0, mx) +
      ylim(0, mx) +
      theme_classic(base_size = 25)

    plot_list[[ct]] <- plot_list[[ct]] + annotate(
      geom = "text",
      x = 0.8, y = 0.1,
      label = paste0("RMSE = ", rmse_value),
      size = 6, family = "sans", parse = FALSE
    )
  }

  common_x_label <- "True cell type proportion in sample"
  common_y_label <- "Estimated cell type proportion in sample"
  common_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = length(cts))

  y.grob <- grid::textGrob(common_y_label, gp=grid::gpar(fontsize=20), rot=90)
  x.grob <- grid::textGrob(common_x_label, gp=grid::gpar(fontsize=20))

  gridExtra::grid.arrange(gridExtra::arrangeGrob(common_plot, left = y.grob, bottom = x.grob))
}

#' Plot predicted/true scatter plot for basis (matrix W)
#'
#' Uses ggplot with rasterization to plot and  gridExtra::grid to arrange plots
#'
#' @param ptb result of coerce_pred_true_props or list of two matrices with the same order of rows
#' @param max_expr limit for the scale. set manually if you want to adjust how much is shown
#' @param pt_alpha alphavalue for points
#' @return list of plots grid extra
#' @importFrom Metrics rmse
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob
#' @importFrom grDevices adjustcolor
#' @export
plot_ptb_scatter <- function(ptb, max_expr = 1, pt_alpha = 0.2) {
  plot_list <- list()
  cts <- colnames(ptb[[1]])
  for (ct in cts) {
    to_plot <- as.data.frame(list(
      DualSimplex = log2(ptb[[1]][, ct] + 1),
      other = log2(ptb[[2]][, ct] + 1)
    ))

    # Calculate RMSE value
    rmse_value <- round(Metrics::rmse(to_plot$DualSimplex, to_plot$other), 2)

    mx <- max(max_expr, max(max(to_plot$other, to_plot$DualSimplex))) + 1
    plot_list[[ct]] <- ggplot(to_plot, aes(x = .data$other, y = .data$DualSimplex)) +
      rasterize_if_needed(geom_point(colour = adjustcolor("black", alpha.f = pt_alpha))) +
      annotate("text", x = mx, y = mx/10, label = paste0("RMSE = ", rmse_value),
               hjust = 1, vjust = 0, size = 6, family = "sans") +
      ggtitle(ct) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1.5) +
      xlab(NULL) +
      ylab(NULL) +
      xlim(0, mx) +
      ylim(0, mx) +
      theme_classic(base_size = 25)
  }


  common_x_label <- "Actual log2 gene expression"
  common_y_label <- "Estimated log2 gene expression"

  common_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = length(cts))

  y.grob <- grid::textGrob(common_y_label, gp = grid::gpar(fontsize = 20), rot = 90)

  x.grob <- grid::textGrob(common_x_label, gp = grid::gpar(fontsize = 20))

  gridExtra::grid.arrange(gridExtra::arrangeGrob(common_plot, left = y.grob, bottom = x.grob))
}

#' Plot correlation heatmap for two matrices.
#' Plots correlation of rows.
#'
#' @param estimated_matrix Matrix we produced.
#' @param true_matrix Matrix we compare to.
#' @export
plot_correlation_matrix <- function(estimated_matrix, true_matrix) {
  gradient_colors <-  c("#798234FF", "#A3AD62FF", "#D0D3A2FF", "#FDFBE4FF", "#F0C6C3FF", "#DF91A3FF", "#D46780FF")
  invis <- element_blank()
  estimated_matrix <- t(estimated_matrix)
  true_matrix <- t(true_matrix)

  plot.data <- stats::cor(estimated_matrix, true_matrix, method = "p")

  plot.data <- as.data.frame(plot.data)
  plot.data[['estimated']] <- colnames(estimated_matrix)
  plot.data[['estimated']] <-  factor( plot.data[['estimated']], levels =  plot.data[['estimated']])
  plot.data <- reshape2::melt(plot.data, value.name = "Correlation",id="estimated")

  result_plot <-   ggplot(plot.data,
                          aes(x = .data[['variable']], y = .data[['estimated']],
                              fill =  .data[['Correlation']],
                              label = round( .data[['Correlation']], 2))) +
    geom_tile(colour = "black") +
    geom_text(size = 2.5) +
    theme_bw() +
    scale_fill_gradientn(colors=gradient_colors, limits =c(-1, 1)) +
    labs(x = "True Component", y = "Predicted Component") +
    theme(panel.border = invis,
          axis.ticks = invis,
          legend.position = "none",
          panel.grid = invis,
          axis.text.x = element_text(angle = 45, hjust = 1))
    return(result_plot)
}

#' Reorder rows  to match each other and plot correlation matrix
#'
#' @param estimated_matrix Matrix we produced.
#' @param true_matrix Matrix we compare to.
#' @param normalize Flag wether we should colum normalize matrices before plotting. Default TRUE.
#' @return list containing plot and cell type mapping.
#' @export
reorder_and_plot_correlation <- function(estimated_matrix, true_matrix, normalize=TRUE) {
  estimated_matrix <-  estimated_matrix[, colnames(true_matrix)]
  new_row_order <-  guess_order(estimated_matrix, true_matrix)
  reordered_matrix <-  estimated_matrix[unlist(new_row_order), ]
  rownames(reordered_matrix) <-  names(new_row_order)

  if (normalize) {
      reordered_matrix <- t(t(reordered_matrix)/colSums(reordered_matrix))
      true_matrix <- t(t(true_matrix)/colSums(true_matrix))
  }

  result_plot <- plot_correlation_matrix(reordered_matrix, true_matrix)
  cell_type_mapping <- as.list(rownames(estimated_matrix)[unlist(new_row_order)])
  names(cell_type_mapping) <-  names(new_row_order)
  return(list(correlation_plot = result_plot, cell_type_mapping=cell_type_mapping))
}


