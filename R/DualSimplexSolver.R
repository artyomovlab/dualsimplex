#' @title DualSimplex Object.
#'
#' @description
#' The DualSimplex Object Class.
#' Provides an interface to perform:
#' sinkhorn transformation,
#' svd projection,
#' optimization of simplex corners (NMF solution).
#'
#' @docType class
#' @importFrom R6 R6Class
#' @useDynLib DualSimplex
#' @export
#' @return Object of \code{\link{R6Class}} -- an interface to work with data.
#' @format \code{\link{R6Class}} object.
#' @examples
#' M <-  8000 # number of genes (rows)
#' N <-  200 # number of samples (columns)
#' K <-  3 # number of main components
#' sim <- create_simulation(n_genes = M, n_samples = N, n_cell_types = K, with_marker_genes = FALSE)
#' dso <- DualSimplexSolver$new()
#' dso$set_data(sim$data) # run Sinkhorn procedure
#' dso$project(K) # project to SVD space
#' dso$plot_projected("zero_distance",
#'                    "zero_distance",
#'                    with_solution = TRUE,
#'                    use_dims = list(2:3)) # visualize the projection
#'
#' @name DualSimplexObject
#' @import Rcpp
#' @import RcppArmadillo
#' @import dplyr
#' @import ggplot2
#' @import Matrix
#' @import Biobase
#' @import irlba
#' @import knitr
#'
DualSimplexSolver <- R6Class(
  classname = "DualSimplexSolver",
  private = list(
    display_dims = NULL,
    save_dir = NULL,
    reset_since = function(step) {
      private$check_step_name(step)
      reset_sel <- seq_along(self$st) >= which(names(self$st) == step)
      self$st[reset_sel] <- list(NULL)
    },
    check_step_name = function(step) {
      if (!step %in% names(self$st)) {
        stop(paste("Internal error, steps misconfigured:", step))
      }
    },
    set_data_first = function() {
      if (is.null(self$st$data)) {
        stop("Call set_data first")
      }
    },
    project_first = function() {
      if (is.null(self$st$n_cell_types)) {
        stop("Call project first (see plot_svd for n_cell_types argument)")
      }
    },
    initialize_first = function() {
      if (is.null(self$st$solution_proj)) {
        stop("Call init_solution first")
      }
    },
    optimize_first = function() {
      if (is.null(self$st$solution_proj$optim_history)) {
        stop("Call optim_solution first")
      }
    },
    finalize_first = function() {
      if (is.null(self$st$solution)) {
        stop("Call finalize_solution first")
      }
    },
    add_filtering_log_step = function(step_name, filtering_params = NA) {
      if (is.null(self$st$filtering_log)) {
        self$st$filtering_log <- list(
          stats_df = data.frame(matrix(ncol = 4, nrow = 0)),
          object_log = list()
        )
        colnames(self$st$filtering_log$stats_df) <- c(
          "step_name",
          "n_genes",
          "n_samples",
          "filtering_params"
        )
      }

      step_num <- length(self$st$filtering_log$object_log) + 1

      self$st$filtering_log$stats_df[step_num, ] <- c(
        paste(step_num, step_name, sep = "_"),
        nrow(self$st$data),
        ncol(self$st$data),
        filtering_params
      )
      self$st$filtering_log$object_log[[step_num]] <- list(
        Sigma = self$st$proj$meta$Sigma
      )
    },
    resolve_color_col = function(color, genes) {
      private$set_data_first()
      anno <- get_anno(self$st$data, genes)
      if (is.character(color) && length(color) == 1) {
        if (!is.null(self$st$solution) && (color %in% colnames(self$st$solution$W))) {
          name <- color
          if (genes) {
            color <- log(self$st$solution$W[, color] + 1)
          } else {
            color <- self$st$solution$H[color, ]
          }
        } else if (genes && !is.null(self$st$marker_genes) && (color == "markers")) {
          name <- color
          color <- which_marker(rownames(self$st$data), self$st$marker_genes)
        } else if (color %in% colnames(anno)) {
          name <- color
          color <- anno[, color]
        } else {
          name <- NULL
        }
      } else {
        name <- NULL
      }
      return(list(color = color, name = name))
    },

    check_max_dim = function(dims, max_dim = self$st$proj_ops$max_dim) {
      private$set_data_first()
      if (max(dims) > max_dim) stop("Not enough dimension in ops. Run `calc_svd_ops` with larger max_dim parameter")
    }
  ),
  public = list(
    #' @field st contain the "state" of the current object. (data, solution, projections etc..).
    st = list(
      filtering_log = NULL,   # Auto calculated
      data = NULL,            # Set by user
      scaling = NULL,         # Auto calculated
      proj_ops = NULL,        # Auto calculated
      n_cell_types = NULL,    # Set by user
      dims = NULL,            # Auto calculated
      proj = NULL,            # Auto calculated, proj$umap is triggered by user
      solution_proj = NULL,   # Triggered by user
      solution = NULL,        # Triggered by user
      solution_orig = NULL,   # Auto calculated
      marker_genes = NULL     # Auto calculated
    ),

    #' @description
    #' Set data to the object.
    #' In general it could be any matrix with names on columns and rows. Expression set will be created.
    #'
    #' @param data input data matrix
    #' @param gene_anno_lists named list of lists. Each sublist contains names of rows which should have TRUE value in annotaiton column.
    #' @param sample_anno_lists named list of lists. Each sublist contains names of columns which should have TRUE value in annotation column.
    #' @param sinkhorn_iterations number of sinkhorn iterations to perform
    set_data = function(data, gene_anno_lists = NULL, sample_anno_lists = NULL, sinkhorn_iterations=20, max_dim = 50L, tol = 1e-05) {
      if (any(sapply(dimnames(data), is.null)))
        stop("Genes and samples should be named")
      if (any(sapply(dimnames(data), anyDuplicated)))
        stop("Gene and sample names should not contain duplicates")
      if (any(rowSums(as.matrix(data)) == 0))
        stop("The data matrix should not contain all zero rows. Use remove_zero_rows() method")
      if (any(colSums(as.matrix(data)) == 0))
        stop("The data matrix should not contain all zero columns. Use remove_zero_cols() method")
      first_set <-  is.null(self$st$data)
      private$reset_since("data")
      if (!inherits(data, "ExpressionSet")) data <- create_eset(data)
      self$st$data <- add_default_anno(data, gene_anno_lists, sample_anno_lists)
      self$st$scaling <- sinkhorn_scale(exprs(self$st$data), max_iter = sinkhorn_iterations)
      self$st$proj_ops <- calc_svd_ops(self$get_V_row(), max_dim = max_dim, tol = tol)
      self$st$proj <- efficient_svd_project(self$get_V_row(), self$get_V_column(), dims = NULL, ops = self$st$proj_ops)
      if (first_set) private$add_filtering_log_step("initial")
    },

    #' @description
    #' plot MAD distribution for the data
    #' In general it could be any matrix with names on columns and rows. Expression set will be created.
    #'
    #' @param ... parameters to pass to plot_feature method.
    plot_mad = function(...) {
      plot_feature(self$get_data(), "log_mad", ...)
    },
    #' @description
    #' Basic data filtering for gene expression datasets.
    #' Removes selected genes, filters by mad
    #'
    #' @param log_mad_gt log madthreshold to remove genes (we remove low mad genes).
    #' @param remove_true_cols_default FALSE if tou don't want to use default gene names filter.
    #' @param remove_true_cols_additional additional columns from annotation to use for "remove true" filter.
    #' @param keep_true_cols columns from annotation where we should keep instances with true value.
    #' @param genes true if want remove rows otherwise columns
    basic_filter = function(
      log_mad_gt = 0,
      remove_true_cols_default = NULL,
      remove_true_cols_additional = c(),
      keep_true_cols = c(),
      genes = T
    ) {
      if (genes && is.null(remove_true_cols_default)) {
        remove_true_cols_default <- c("RPLS", "LOC", "ORF", "SNOR")
      }
      private$set_data_first()
      # private$check_filtering_log()
      remove_true_cols <- c(remove_true_cols_default, remove_true_cols_additional)
      new_data <- self$get_data()
      new_data <- threshold_filter(new_data, "log_mad", log_mad_gt, genes = genes, keep_lower = F)
      new_data <- bool_filter(new_data, remove_true_cols, genes = genes, remove_true = T)
      new_data <- bool_filter(new_data, keep_true_cols, genes = genes, remove_true = F)
      self$set_data(new_data)
      private$add_filtering_log_step(
        "basic_filters",
        paste(
          paste0("remove_true_cols = ", paste(remove_true_cols, collapse = " ")),
          paste0("keep_true_cols = ", paste(keep_true_cols, collapse = " ")),
          paste0("log_mad_gt = ", log_mad_gt),
          sep = ", "
        )
      )
    },
    #' @description
    #' Add additional distance annotation based on KNN distances to selected annotations.
    #'
    #' @param annotation_names_list names of annotation columns with TRUE/FALSE.
    #' @param genes calculate for gene annotations or sample annotations.
    #' @param k_neighbors a number of neighbors to calculate the distance on for the annotation
    add_knn_distances_anno= function(annotation_names_list = NULL, genes = T, k_neighbors = 20) {
      self$st$data <- add_knn_distances_anno(
        self$st$data,
        self$st$proj_full,
        annotation_columns = annotation_names_list,
        genes = genes,
        k_neighbors = k_neighbors
      )
    },


    #' @description
    #' Interface to plot svd
    #' Will return the elbow plot of singular values.
    #'
    #' @param dims how many dimensions to plot
    #' @return plot to work with
    plot_svd = function(dims = self$st$dims) {
      private$set_data_first()
      private$check_max_dim(dims)
      plot_proj_svd(self$st$proj, dims)
    },

    #' @description
    #' Interface to plot svd history of filtering steps through the svd plot.
    #'
    #' @param steps_sel selected filtering steps, list.
    #' @param n_dims how many dimensions to plot
    #' @param cumulative wether plot should be cumulative
    #' @param variance plot variance explained
    #' @return plot to work with
    plot_svd_history = function(
      steps_sel = NULL,
      n_dims = NULL,
      cumulative = T,
      variance = T
    ) {
      svd_ds <- lapply(self$st$filtering_log$object_log, function(x) {
        diag(x$Sigma)
      })
      mlen <- max(sapply(svd_ds, length))
      svd_ds <- lapply(svd_ds, function(x) {
        c(x, rep(0, mlen - length(x)))
      })
      svd_ds <- matrix(unlist(svd_ds), ncol = mlen, byrow = T)
      rownames(svd_ds) <- self$st$filtering_log$stats_df$step_name
      if (!is.null(steps_sel)) {
        svd_ds <- svd_ds[steps_sel, ]
      }
      if (!is.null(n_dims)) {
        svd_ds <- svd_ds[, 1:n_dims]
      }
      return(plot_svd_ds_matrix(svd_ds, cumulative = cumulative, variance =  variance))
    },

    #' @description
    #' Calculate svd projection for current data. Will perform SVD.
    #'
    #' @param n_cell_types selected number of dimensions (K) to work with.
    project = function(n_cell_types) {
      private$set_data_first()
      private$reset_since("n_cell_types")
      self$st$n_cell_types <- n_cell_types
      self$st$dims <- if (!is.null(n_cell_types)) 1:n_cell_types else NULL
      self$st$proj <- efficient_svd_project(self$get_V_row(), self$get_V_column(), , dims = self$st$dims, ops = self$st$proj_ops)
      self$st$data <- add_distances_anno(
        self$st$data,
        list(
          "V_row" = self$get_V_row(),
          "V_column" = self$get_V_column()
        ),
        self$st$proj,
        self$st$n_cell_types
      )
    },


    #' @description
    #' A set of plots to extimate the projection.
    plot_projection_diagnostics = function() {
      plt1 <- self$plot_projected("zero_distance", "zero_distance")
      plt2 <- self$plot_projected("plane_distance", "plane_distance")
      plt3 <- self$plot_distances_distribution()
      plt4 <- self$plot_svd()
      return(list(plt1, plt2, plt3, plt4))
    },

    #' @description
    #' Plot scatterplot of plane distance to zero distance for each point
    plot_distances_distribution = function() {
      private$project_first()
      # TODO: maybe self$st$proj should be ExpressionSet with distances
      show(plot_numeric_features(self$get_data(), features = c("plane_distance", "zero_distance"), ncol = 2, labels = c("Plane distance", "Zero distance")))
      show(plot_numeric_features(self$get_data(), genes = F, features = c("plane_distance", "zero_distance"), ncol = 2, labels = c("Plane distance", "Zero distance"), bins = 30))
      # TODO: ordered plots
      plotlist <- list(
        plot_feature_pair(self$get_data(), "plane_distance", "zero_distance", T, size = 0.1),
        plot_feature_pair(self$get_data(), "plane_distance", "zero_distance", F, size = 0.1)
      )
      show(cowplot::plot_grid(plotlist = plotlist, ncol = 2))
    },

    #' @description
    #' Filter points based on distance thresholds. (keep lower).
    #'
    #' @param plane_d_lt threshold for plane distance.
    #' @param zero_d_lt threshold for zero distance.
    #' @param genes TRUE if filter rows, otherwise columns.
    distance_filter = function(
      plane_d_lt = NULL,
      zero_d_lt = NULL,
      genes = T
    ) {
      private$project_first()
      if (is.null(plane_d_lt) && is.null(zero_d_lt)) {
        stop("Choose at least one distance to filter by")
      }
      new_data <- self$get_data()

      if (!is.null(plane_d_lt))
        new_data <- threshold_filter(
          new_data,
          "plane_distance",
          plane_d_lt,
          genes
        )
      if (!is.null(zero_d_lt))
        new_data <- threshold_filter(
          new_data,
          "zero_distance",
          zero_d_lt,
          genes
        )
      new_data <- remove_zero_cols(new_data)
      new_data <- remove_zero_rows(new_data)

      self$set_data(new_data)
      private$add_filtering_log_step(
        "distance_filter",
        paste(
          paste0("zero_d_lt = ", zero_d_lt),
          paste0("plane_d_lt = ", plane_d_lt),
          sep = ", "
        )
      )
    },

    #' @description
    #' Do UMAP transformation for current projected data in both spaces.
    #'
    #' @param with_model specific umap model selection for Mac users. For Mac there is a known issue with this library.
    #' @param neighbors_X parameter for UMAP.
    #' @param neighbors_Omega parameter for UMAP.
    run_umap = function(with_model = Sys.info()[["sysname"]] != "Darwin", neighbors_X = 50, neighbors_Omega = 10) {
      private$project_first()
      self$st$proj <- add_proj_umap(self$st$proj, with_model, neighbors_X, neighbors_Omega)
    },

    #' @description
    #' Interface to plot points of the current object
    #'
    #'  color_genes / color_samples can be:
    #' - a set of names to be highlighted
    #' - a vector of values, the same length as the number of genes
    #' - a name of a column from annotation, default is zero_distance
    #' Important note: ggplot only allows to either draw history or color all the points
    #'
    #' @param color_genes how to color genes (see description of method).
    #' @param color_samples how to color samples (see description of method).
    #' @param use_dims which dimensions to use (e.g. 2:3).
    #' @param with_legend TRUE if want to add legends to plots.
    #' @param with_solution TRUE if want to  add current solution points.
    #' @param with_history TRUE if want to  add  solution history points/lines.
    #' @param wrap FALSE if want to have two separated plots, not single one.
    #' @param show_plots FALSE if don't want plots to be displayed.
    #' @param from_iter starting point for history of solutions.
    #' @param to_iter end point for history of solutions.
    #' @param ... any other params to be passed to plot_projected method.
    plot_projected = function(
      color_genes = "zero_distance", color_samples = "zero_distance",
      use_dims = private$display_dims, with_legend = NULL,
      with_solution = TRUE, with_history = TRUE,
      wrap = T, show_plots = T, from_iter = 1, to_iter = NULL, ...
    ) {
      if (inherits(use_dims, "list")) {
        plotlist <- lapply(use_dims, function(this_use_dims) {
          self$plot_projected(
            color_genes, color_samples, this_use_dims,
            with_legend, with_solution, with_history,
            wrap, show_plots, from_iter, to_iter, ...
          )
        })
        if (show_plots) {
          for (plot in plotlist) { show(plot) }
          return(invisible(NULL))
        } else {
          return(plotlist)
        }
      }
      private$project_first()

      color_genes <- private$resolve_color_col(color_genes, T)
      color_samples <- private$resolve_color_col(color_samples, F)

      plt_X <- plot_projection_points(self$st$proj, use_dims, "X", color = color_genes$color, color_name = color_genes$name, ...)
      plt_Omega <- plot_projection_points(self$st$proj, use_dims, "Omega", color = color_samples$color, color_name = color_samples$name, ...)

      if (!is.null(self$st$solution_proj)) {
        if (("optim_history" %in% names(self$st$solution_proj)) && with_history) {
          plt_X <- plt_X %>% add_solution_history(
            self$st$solution_proj,
            self$st$proj,
            use_dims = use_dims,
            spaces = "X",
            colored = is.null(color_genes$name),
            from_iter = from_iter,
            to_iter = to_iter
          )
          plt_Omega <- plt_Omega %>% add_solution_history(
            self$st$solution_proj,
            self$st$proj,
            use_dims = use_dims,
            spaces = "Omega",
            colored = is.null(color_samples$name),
            from_iter = from_iter,
            to_iter = to_iter
          )
        }

        if (with_solution) {
          plt_X <- plt_X %>% add_solution(
            self$st$solution_proj,
            self$st$proj,
            use_dims = use_dims,
            spaces = "X"
          )

          plt_Omega <- plt_Omega %>% add_solution(
            self$st$solution_proj,
            self$st$proj,
            use_dims = use_dims,
            spaces = "Omega"
          )
        }
      }

      if (!is.null(with_legend) && with_legend) {
          plt_X <- plt_X + theme(legend.position = "right")
        plt_Omega <- plt_Omega + theme(legend.position = "right")
      }

      plt_X <- plt_X + ggtitle(paste(self$st$proj$meta$M, "genes"))
      plt_Omega <- plt_Omega + ggtitle(paste(self$st$proj$meta$N, "samples"))

      plotlist <- list(plt_X, plt_Omega)
      return(if (wrap) cowplot::plot_grid(plotlist = plotlist) else plotlist)
    },

    #' @description
    #' Initialize current solution
    #'
    #' @param strategy strategy to use for initialization. valid values are "select_x", "select_omega", "random" and "marker_means"
    #' @param ... any other params to be passed to initialization methods (e.g. marker genes)
    init_solution = function(strategy = "select_x", ...) {
      private$project_first()
      private$reset_since("solution_proj")
      kwargs <- list(...)
      self$st$solution_proj <- initialize_solution(self$st$proj, strategy, kwargs)
    },

    #' @description
    #' Perform optimization for current solution
    #'
    #' @param iterations number of steps to perform
    #' @param config optimization config (result of optim_config method)
    optim_solution = function(
      iterations = 10000,
      config = OPTIM_CONFIG_DEFAULT
    ) {
      private$initialize_first()
      self$st$solution_proj <- optimize_solution(
        self$st$proj,
        self$st$solution_proj,
        iterations,
        config
      )
    },

    #' @description
    #' Plot errors log
    plot_error_history = function() {
      private$optimize_first()
      plot_errors(self$st$solution_proj)
    },

    #' @description
    #' Finalize solution.
    #' Return from projection to sinkhorn transformed matrices
    #' Perform reverse sinkhorn to get original matrices W and H
    finalize_solution = function() {
      private$initialize_first()
      solution_scaled <- reverse_solution_projection(self$st$solution_proj, self$st$proj)
      self$st$solution_no_corr <- reverse_solution_sinkhorn(solution_scaled, self$st$scaling)
      self$st$solution <- list(
        W = self$st$solution_no_corr$W,
        H = self$st$solution_no_corr$H
      )
      self$st$solution$W[self$st$solution$W < 0] <- 0
      self$st$solution$H[self$st$solution$H < 0] <- 0

      self$st$marker_genes <- get_signature_markers(self$st$solution$W)
      return(self$st$solution)
    },
    #' @description
    #' Interface to plot negativity changes in basis (Matrix W)
    plot_negative_basis_change = function() {
      private$optimize_first()
      return(plot_negative_basis_change(self$st$proj, self$st$solution_proj))
    },

    #' @description
    #' Interface to plot negativity changes in coefficients (Matrix H)
    plot_negative_proportions_change = function() {
      private$optimize_first()
      return(plot_negative_proportions_change(self$st$proj, self$st$solution_proj))
    },

    #' @description
    #' Interface to plot solution distribution as histogramm by main component
    plot_solution_distribution = function() {
      private$finalize_first()
      cowplot::plot_grid(plotlist = list(
        plot_proportions_distribution(self$get_solution()$H) + ggtitle(
          "Proportions, H",
          paste0(
            round(self$get_negative_ratios()$H * 100, 2),
            "% matrix entries were negative before setting to zero"
          )
        ),
        plot_basis_distribution(self$get_solution()$W) + ggtitle(
          "Basis, W",
          paste0(
            round(self$get_negative_ratios()$W * 100, 2),
            "% matrix entries were negative before setting to zero"
          )
        )
      ))
    },

    #' @description
    #' set and remember which dimensions to use for plotting functions
    #'
    #' @param display_dims  (e.g. 2:3, 3:4, NULL)
    set_display_dims = function(display_dims) {
      private$display_dims <- display_dims
    },

    ##### Saving to/Loading from files #####
    #' @description
    #' set and remember direvtory to save model state
    #'
    #' @param new_dir_path  path to save model
    set_save_dir = function(new_dir_path) {
      if (!dir.exists(new_dir_path)) {
        dir.create(new_dir_path, recursive = TRUE)
      }
      private$save_dir <- normalizePath(new_dir_path)
    },

    #' @description
    #' get current save directory for the model
    get_save_dir = function() {
      return(private$save_dir)
    },

    #' @description
    #' set and remember direvtory to save model state. returns directory name.
    #'
    #' @param new_dir_path  path to save model
    getset_save_dir = function(new_dir_path = NULL) {
      if (is.null(private$save_dir)) {
        if (is.null(new_dir_path)) stop("Specify save_dir or call set_save_dir")
        self$set_save_dir(new_dir_path)
      } else if (!is.null(new_dir_path) && !private$save_dir == new_dir_path) {
        stop("save_dir is not NULL and can only be changed via set_save_dir")
      }
      return(private$save_dir)
    },

    #' @description
    #' save model state to directory
    #'
    #' @param save_dir  path to save model
    save_state = function(save_dir = NULL) {
      out_dir <- self$getset_save_dir(save_dir)
      saveRDS(self$st, file.path(out_dir, "dualsimplex_state.rds"))
      if (("umap" %in% names(self$st$proj)) && (!is.null(self$st$proj$umap$model))) {
        fpath <- file.path(out_dir, "dualsimplex_umap_X.uwot")
        if (file.exists(fpath)) file.remove(fpath)
        uwot::save_uwot(self$st$proj$umap$model$X, fpath)

        fpath <- file.path(out_dir, "dualsimplex_umap_Omega.uwot")
        if (file.exists(fpath)) file.remove(fpath)
        uwot::save_uwot(self$st$proj$umap$model$Omega, fpath)
      }
      if (!is.null(self$st$solution)) {
        self$save_solution(out_dir)
      }
      return(invisible())
    },

    #' @description
    #' save current model solution (H and W) to separate tsv files
    #'
    #' @param save_dir  path to save matrices
    save_solution = function(save_dir = NULL) {
      private$finalize_first()
      out_dir <- self$getset_save_dir(save_dir)
      W <- cbind(gene_name = rownames(self$st$solution$W), self$st$solution$W)
      H <- cbind(sample_name = colnames(self$st$solution$H), t(self$st$solution$H))
      write.table(W, file.path(out_dir, "basis.tsv"), sep = "\t", quote = F, row.names = F)
      write.table(H, file.path(out_dir, "proportions.tsv"), sep = "\t", quote = F, row.names = F)
    },

    #' @description
    #' load model from the specified directory.
    #'
    #' @param input_dir  path to load model from
    load_state = function(input_dir = NULL) {
      if (is.null(input_dir)) {
        if (!is.null(private$save_dir)) {
          input_dir <- private$save_dir
        } else {
          stop("Current save_dir is null, specify input_dir or set_save_dir")
        }
      }
      self$st <- readRDS(file.path(input_dir, "dualsimplex_state.rds"))
      if (("umap" %in% names(self$st$proj)) && (!is.null(self$st$proj$umap$model))) {
        self$st$proj$umap$model$X <- uwot::load_uwot(file.path(input_dir, "dualsimplex_umap_X.uwot"))
        self$st$proj$umap$model$Omega <- uwot::load_uwot(file.path(input_dir, "dualsimplex_umap_Omega.uwot"))
      }
    },

    #' @description
    #' Generate automatically generated report for the method
    #'
    #' @param save_dir  path to save report
    #' @param seurat_obj  seurat object to use for markers visualization
    #' @param with_animated_optim TRUE to save gif with optimization process
    generate_summary = function(
      save_dir = NULL,
      seurat_obj = NULL,
      with_animated_optim = F
    ) {
      out_dir <- self$getset_save_dir(save_dir)

      rmd_path <- system.file("extdata", "solver_summary.Rmd", package = "DualSimplex")
      output_file <- file.path(out_dir, "dualsimplex_solver_summary.html")

      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(out_dir)

      rmarkdown::render(
        input = rmd_path,
        output_file = output_file,
        params = list(
          dso = self,
          with_animated_optim = with_animated_optim,
          seurat_obj = seurat_obj
        )
      )
      if (!is.null(self$st$solution)) {
        self$save_solution(out_dir)
      }
    },

    ##### Getters #####
    #' @description
    #' get current data.
    get_data = function() {
      private$set_data_first()
      return(self$st$data)
    },
    #' @description
    #' get current filtering log.
    get_filtering_stats = function() {
      return(self$st$filtering_log$stats_df)
    },
    #' @description
    #' get number of main components used (K)
    get_n_cell_types = function() {
      private$project_first()
      return(self$st$n_cell_types)
    },
    #' @description
    #' get number optimization iterations performed
    get_n_iters = function() {
      private$optimize_first()
      return(nrow(dso$st$solution_proj$optim_history$errors_statistics))
    },
    #' @description
    #' get proportionality of negative elements for H and W
    get_negative_ratios = function() {
      private$finalize_first()
      return(list(
        W = sum(self$st$solution_no_corr$W < 0) / length(self$st$solution_no_corr$W),
        H = sum(self$st$solution_no_corr$H < 0) / length(self$st$solution_no_corr$H)
      ))
    },

    #' @description
    #' get current solution
    get_solution = function() {
      private$finalize_first()
      return(self$st$solution)
    },
    #' @description
    #' get names for main components specified
    get_ct_names = function() {
      private$finalize_first()
      return(rownames(self$st$solution$H))
    },

    #' @description
    #' get marker genes for current solution
    get_marker_genes = function() {
      private$finalize_first()
      return(self$st$marker_genes)
    },

    #' @description
    #' Calculate V_row on the fly.
    #' TODO: It might be good to have a more efficient implementation, if this function will be frequently called.
    get_V_row = function(scaling) {
      private$set_data_first()

      dr <- matrixStats::rowProds(self$st$scaling$D_vs_row)
      dc <- matrixStats::rowProds(self$st$scaling$D_vs_col[, 1:(self$st$scaling$iterations-1)])
      
      sweep(
        dr * exprs(self$get_data()),  # R is column-oriented, direct multiplication is row-wise,
        MARGIN = 2,
        STATS = dc,
        FUN = `*`
      )
   },

    #' @description
    #' Calculate V_column on the fly.
    #' TODO: It might be good to have a more efficient implementation, if this function will be frequently called.
    get_V_column = function(scaling) {
      private$set_data_first()

      dr <- matrixStats::rowProds(self$st$scaling$D_vs_row)
      dc <- matrixStats::rowProds(self$st$scaling$D_vs_col)
      
      sweep(
        dr * exprs(self$get_data()),  # R is column-oriented, direct multiplication is row-wise,
        MARGIN = 2,
        STATS = dc,
        FUN = `*`
      )
    }
  )
)

#' DualSimplexSolver$from_state
#'
#' Static method to load model
#' @name DualSimplexSolver$from_state
#' @param input_dir directory to read from.
#' @param save_there TRUE to remember directory choise for object.
#' @return dso object
DualSimplexSolver$from_state <- function(input_dir, save_there = F) {
  dso <- DualSimplexSolver$new()
  dso$load_state(input_dir)
  if (save_there) {
    dso$set_save_dir(input_dir)
  }
  return(dso)
}
