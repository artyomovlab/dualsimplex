############ MAIN LOGIC ############

coerce_pred_true_basis <- function(pred_basis, true_basis) {
  pred_basis <- pred_basis[, guess_order(t(pred_basis), t(true_basis))]
  colnames(pred_basis) <- colnames(true_basis)
  rownames(pred_basis) <- rownames(true_basis)
  list(pred_basis, true_basis)
}

coerce_pred_true_props <- function(pred_props, true_props) {
  pred_props <- pred_props[guess_order(pred_props, true_props), ]
  colnames(pred_props) <- colnames(true_props)
  rownames(pred_props) <- rownames(true_props)
  list(pred_props, true_props)
}

guess_order <- function(predicted, actual) {
  ctn <- nrow(predicted)
  allPerms <- combinat::permn(ctn)

  vals <- sapply(allPerms, function(perm) {
    sum(abs(predicted[perm, ] - actual))
  })
  perm <- allPerms[[which.min(vals)]]
  perm
}

r2 <- function(a, b) {
  cor(as.vector(a), as.vector(b))^2
}

rmse <- Metrics::rmse


############ PLOTTING ############

plot_ptp_lines <- function(ptp) {
  linseed::plotProportions(ptp[[1]], ptp[[2]], pnames = c("predicted", "true"))
}

plot_ptp_scatter <- function(ptp) {
  plot_list <- list()
  cts <- rownames(ptp[[1]])
  for (ct in cts) {
    to_plot <- as.data.frame(list(
      linseed = ptp[[1]][ct, ],
      other = ptp[[2]][ct, ]
    ))

    rmse_value <- round(rmse(ptp[[1]][ct, ], ptp[[2]][ct, ]), 2)

    mx <- 1

    plot_list[[ct]] <- ggplot(to_plot, aes(x=other, y=linseed)) +
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

plot_ptb_scatter <- function(ptb, max_expr = 25, pt_alpha = 0.2) {
  plot_list <- list()
  cts <- colnames(ptb[[1]])
  for (ct in cts) {
    to_plot <- as.data.frame(list(
      linseed = log2(ptb[[1]][, ct] + 1),
      other = log2(ptb[[2]][, ct] + 1)
    ))

    # Calculate RMSE value
    rmse_value <- round(rmse(to_plot$linseed, to_plot$other), 2)

    mx <- max(max_expr, max(max(to_plot$other, to_plot$linseed))) + 1
    plot_list[[ct]] <- ggplot(to_plot, aes(x = other, y = linseed)) +
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