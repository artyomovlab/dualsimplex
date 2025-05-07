############ MAIN LOGIC ############

#' Generate (prediction,true) pairs for each basis vector
#'
#' Will try to gess the right order of columns (since algorithm can shuffle it) for matrix W
#'
#' @param pred_basis  MxK predicted matrix W
#' @param true_basis MxK true matrix W
#' @return list of 2 matrices with ordered columns (1-1, 2-2, ..)
#' @export
coerce_pred_true_basis <- function(pred_basis, true_basis) {
  pred_basis <- pred_basis[, guess_order(t(pred_basis), t(true_basis))]
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
  pred_props <- pred_props[guess_order(pred_props, true_props), ]
  colnames(pred_props) <- colnames(true_props)
  rownames(pred_props) <- rownames(true_props)
  list(pred_props, true_props)
}


#' Function to gess order of rows
#'
#' Just calculates pairwise distance and takes argmin for each column.
#'
#' @param predicted predicted matrix
#' @param actual true matrix
#' @return reordered row names
#' @import combinat
#' @export
guess_order <- function(predicted, actual) {
  ctn <- nrow(predicted)
  allPerms <- combinat::permn(ctn)

  vals <- sapply(allPerms, function(perm) {
    sum(abs(predicted[perm, ] - actual))
  })
  perm <- allPerms[[which.min(vals)]]
  perm
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
#' @param max_expr all bigger values will be replaced to this value
#' @param pt_alpha alphavalue for points
#' @return list of plots grid extra
#' @importFrom Metrics rmse
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob
#' @importFrom grDevices adjustcolor
#' @export
plot_ptb_scatter <- function(ptb, max_expr = 25, pt_alpha = 0.2) {
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
      annotate("text", x = 15, y = 2, label = paste0("RMSE = ", rmse_value),
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
