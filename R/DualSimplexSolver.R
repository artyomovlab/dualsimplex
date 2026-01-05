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
#' @return Object of \code{\link[R6]{R6Class}} -- an interface to work with data.
#' @format \code{\link[R6]{R6Class}} object.
#'
#' @name DualSimplexObject
#' @import Rcpp
#' @import RcppArmadillo
#' @import dplyr
#' @import ggplot2
#' @import Matrix
#' @importFrom Biobase exprs AnnotatedDataFrame ExpressionSet fData pData
#' @import irlba
#' @import knitr
#' @import progress
#'
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
        Sigma = self$st$proj_ops$Sigma
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
        } else if (color %in% rownames(anno)) {
          name <-  "individual_highlight"
        } else {
          name <- "direct_single_color"
        }
      } else {
        name <- NULL
      }
      return(list(color = color, name = name))
    },
    check_max_dim = function(dims) {
      private$set_data_first()
      if (max(dims) > self$st$proj_ops$max_dim) stop("Not enough dimension in ops. Run `calc_svd_ops` with larger max_dim parameter")
    },

    update_variables = function(data,gene_anno_lists = NULL, sample_anno_lists = NULL, ...) {
      if (!inherits(data, "ExpressionSet")) data <- create_eset(data)
      if (any(rowSums(Biobase::exprs(data)) == 0))
        stop("The data matrix should not contain all zero rows. Use remove_zero_rows() method")
      if (any(colSums(Biobase::exprs(data)) == 0))
        stop("The data matrix should not contain all zero columns. Use remove_zero_cols() method")
      self$st$data <- add_default_anno(data, gene_anno_lists, sample_anno_lists)
      self$st$scaling <- sinkhorn_scale(Biobase::exprs(self$st$data), max_iter = self$st$max_sinkhorn_iterations, epsilon=self$st$sinkhorn_tol)
      self$st$proj_ops <- calc_svd_ops(self$get_V_row(), max_dim = self$st$max_dim, self$st$svd_method, ...)
    }
  ),
  public = list(
    #' @field st contain the "state" of the current object. (data, solution, projections etc..).
    st = list(
      data = NULL,                # Set by user
      filtering_log = NULL,       # Auto calculated
      max_dim = NULL,             # Can be set by user. Default is 50
      max_sinkhorn_iterations = NULL, # Can be set by user. Default is 20
      sinkhorn_tol = NULL,        # Can be set by user. Default is 1e-15
      svd_method = NULL,          # Can be set by user. Default is 'svd'
      scaling = NULL,             # Auto calculated
      proj_ops = NULL,            # Auto calculated
      n_cell_types = NULL,        # Set by user
      dims = NULL,                # Auto calculated
      proj = NULL,                # Auto calculated, proj$umap is triggered by user
      solution_proj = NULL,       # Triggered by user
      solution = NULL,            # Triggered by user
      solution_orig = NULL,       # Auto calculated
      marker_genes = NULL        # Auto calculated
    ),

    #' @description
    #' Set data to the object.
    #' In general it could be any matrix with names on columns and rows. Expression set will be created.
    #'
    #' @param data input data matrix
    #' @param gene_anno_lists named list of lists. Each sublist contains names of rows which should have TRUE value in annotaiton column.
    #' @param sample_anno_lists named list of lists. Each sublist contains names of columns which should have TRUE value in annotation column.
    #' @param max_sinkhorn_iterations number of sinkhorn iterations to perform.
    #' @param max_dim maximum dimention we want the projection operation. It is passed to `calc_svd_ops` function.
    #' @param sinkhorn_tol tolerance for Sinkhorn calculation. It is passed to `sinkhron_scale` function.
    #' @param svd_method which SVD algorithm to use.
    #' @param ... additional arguments passed to function `run_svd`
    set_data = function(
      data,
      gene_anno_lists = NULL,
      sample_anno_lists = NULL,
      max_sinkhorn_iterations=20,
      max_dim = 50L,
      sinkhorn_tol = 1e-12,
      svd_method = "svd",
      ...
    ) {
      # Sanity checks
      if (any(sapply(dimnames(data), is.null)))
        stop("Genes and samples should be named")
      if (any(sapply(dimnames(data), anyDuplicated)))
        stop("Gene and sample names should not contain duplicates")

      private$reset_since("data")
      self$st$max_sinkhorn_iterations <- max_sinkhorn_iterations
      self$st$max_dim <- max_dim
      if (self$st$max_dim > min(dim(data))) {
        self$st$max_dim <- min(dim(data))
        warning("Provided `max_dim` is bigger than smallest dimention of `data`. Setting `max_dim` to ", self$st$max_dim, ".")
      }
      self$st$sinkhorn_tol <- sinkhorn_tol
      self$st$svd_method <- svd_method

      private$update_variables(data, gene_anno_lists, sample_anno_lists, ...)
      private$add_filtering_log_step("initial")
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
      new_data <- remove_zero_cols(new_data)
      new_data <- remove_zero_rows(new_data)
      private$update_variables(new_data)
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
        self$st$proj,
        annotation_columns = annotation_names_list,
        genes = genes,
        k_neighbors = k_neighbors
      )
    },

    #' @description
    #' Add additional density annotation to points.
    #'
    #' @param radius radius to count neighbors within.
    #' @param genes calculate for gene annotations or sample annotations.
    add_density_anno= function(radius = NULL, genes = T) {
      self$st$data <- add_density_annotation(
        self$st$data,
        self$st$proj,
        radius = radius,
        genes = genes
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
      plot_svd_d(diag(self$st$proj_ops$Sigma), dims)
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
      self$st$proj <- efficient_svd_project(self$get_V_row(), self$get_V_column(), dims = self$st$dims, ops = self$st$proj_ops)
      self$st$data <- add_distances_anno(
        self$st$data,
        V_row = self$get_V_row(),
        self$st$proj
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

      private$update_variables(new_data)
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
    #' Filter points based on distance thresholds. (keep lower).
    #'
    #' @param plane_quantile quantile for plane distance.
    #' @param zero_quantile quantile for zero distance.
    #' @param genes TRUE if filter rows, otherwise columns.
    #' @param keep_lower TRUE if keep lower, FALSE to keep higher
    distance_quantile_filter = function(
      plane_quantile = NULL,
      zero_quantile = NULL,
      genes = T,
      keep_lower = T
    ) {
      private$project_first()
      if (is.null(plane_quantile) && is.null(zero_quantile)) {
        stop("Choose at least one distance to filter by")
      }
      new_data <- self$get_data()

      if (!is.null(plane_quantile))
        new_data <- quantile_filter(
          eset = new_data,
          feature = "plane_distance",
          quant =  plane_quantile,
          keep_lower = keep_lower,
          genes = genes
        )
      if (!is.null(zero_quantile))
        new_data <- quantile_filter(
          eset = new_data,
          feature = "zero_distance",
          quant =  zero_quantile,
          keep_lower = keep_lower,
          genes = genes
        )
      new_data <- remove_zero_cols(new_data)
      new_data <- remove_zero_rows(new_data)

      private$update_variables(new_data)
      private$add_filtering_log_step(
        "distance_quantile_filter",
        paste(
          paste0("zero_quant = ", zero_quantile),
          paste0("plane_quant = ", plane_quantile),
          sep = ", "
        )
      )
    },


    #' @description
    #' Iteratively filter by N sigma using all the features provided.
    #' aplication of  n_sigma_filter <- function(eset, feature, n_sigma = 3, genes = T)
    #' @param features feature names (columns of pData or fData) as a vector.
    #' @param n_sigma number of sigmas to keep.
    #' @param genes TRUE if filter rows, otherwise columns.
    #' @param max_filtering_iterations maximum fitering iterations to be performed
    iterative_n_sigma_filter = function(
      features = NULL,
      n_sigma = 3,
      max_filtering_iterations = 500,
      genes = T
    ) {
      private$project_first()
      if (is.null(features)) {
        stop("Choose feature names from fData and pData columns to filter by")
      }
      new_data <- self$get_data()

      if (!is.null(features)) {
        filtering_iteration <-  1
        previous_count <-  if(genes) dim(new_data)[[1]] else  dim(new_data)[[2]]
        new_count <- -1
        while((new_count < previous_count) && (filtering_iteration < max_filtering_iterations) ) {
          previous_count <-   if(new_count == -1) previous_count else  new_count
          intermediate_count <-  previous_count

          # Filter all features by selected sigma
          cell_types <-  self$st$n_cell_types
          for (current_feature in features) {
              print(paste("Feature:", current_feature))
              new_data <- n_sigma_filter(eset = new_data, feature = current_feature,  n_sigma = n_sigma, genes = genes)
              new_data <- remove_zero_cols(new_data)
              new_data <- remove_zero_rows(new_data)
              new_count <-  if(genes) dim(new_data)[[1]] else  dim(new_data)[[2]]
              print(paste("removed", intermediate_count - new_count, "points"))
          }
          new_count <-  if(genes) dim(new_data)[[1]] else  dim(new_data)[[2]]
          private$update_variables(new_data)
          self$project(cell_types)
          new_data <- self$get_data()
          filtering_iteration <-  filtering_iteration + 1
          print(paste("Total removed", previous_count - new_count, "points"))

        }
      private$add_filtering_log_step(
        "n_sigma_filter",
        paste(
          paste('features =', paste0(features, collapse=',')),
          paste0("n_sigma = ", n_sigma),
          paste0("iterations = ", filtering_iteration),
          sep = ", "
        )
      )
      }
    },


    #' @description
    #' Iteratively filter by N sigma using all the features provided.
    #' aplication of  n_sigma_filter <- function(eset, feature, n_sigma = 3, genes = T)
    #' @param threshold threshold to filter neighborhoods
    #' @param density_radius radius for density calculation.
    #' @param genes TRUE if filter rows, otherwise columns.
    #' @param max_filtering_iterations maximum fitering iterations to be performed
    iterative_density_filter = function(
      threshold = 0,
      max_filtering_iterations = 500,
      density_radius = NULL,
      genes = T
    ) {
      private$project_first()
      new_data <- self$get_data()
      feature <- 'density'
      filtering_iteration <-  1
      previous_count <-  if(genes) dim(new_data)[[1]] else  dim(new_data)[[2]]
      new_count <- -1
      while((new_count < previous_count) && (filtering_iteration < max_filtering_iterations) ) {
        previous_count <-   if(new_count == -1) previous_count else  new_count
        cell_types <-  self$st$n_cell_types
        new_data <- threshold_filter(eset = new_data, feature = feature,  threshold = threshold, genes = genes, keep_lower = F)
        new_data <- remove_zero_cols(new_data)
        new_data <- remove_zero_rows(new_data)
        new_count <-  if(genes) dim(new_data)[[1]] else  dim(new_data)[[2]]
        private$update_variables(new_data)
        self$project(cell_types)
        self$add_density_anno(radius=density_radius, genes=genes)
        new_data <- self$get_data()
        filtering_iteration <-  filtering_iteration + 1
        print(paste("removed", previous_count - new_count, "points"))
      }
      private$add_filtering_log_step(
        "iterative_density_filter",
        paste(
          paste('radius =', density_radius),
          paste0("threshold = ", threshold),
          paste0("iterations = ", filtering_iteration),
          sep = ", "
        )
      )
    },

    #' @description
    #' Iteratively filter by N sigma using all the features provided simultaniously.
    #' aplication of  n_sigma_filter <- function(eset, feature, n_sigma = 3, genes = T)
    #' @param features feature names (columns of pData or fData) as a vector.
    #' @param n_sigma number of sigmas to keep.
    #' @param genes TRUE if filter rows, otherwise columns.
    #' @param max_filtering_iterations maximum fitering iterations to be performed
    iterative_mahalanobis_filter = function(
      features = NULL,
      n_sigma = 3,
      max_filtering_iterations = 500,
      genes = T
    ) {
      private$project_first()
      if (is.null(features)) {
        stop("Choose feature names from fData and pData columns to filter by")
      }
      new_data <- self$get_data()

      if (!is.null(features)) {
        filtering_iteration <-  1
        previous_count <-  if(genes) dim(new_data)[[1]] else  dim(new_data)[[2]]
        new_count <- -1
        while((new_count < previous_count) && (filtering_iteration < max_filtering_iterations) ) {
        previous_count <-   if(new_count == -1) previous_count else  new_count
        cell_types <-  self$st$n_cell_types
        new_data <- mahalanobis_n_sigma_filter(eset = new_data, features = features,  n_sigma = n_sigma, genes = genes)
        new_data <- remove_zero_cols(new_data)
        new_data <- remove_zero_rows(new_data)
        new_count <-  if(genes) dim(new_data)[[1]] else  dim(new_data)[[2]]
        private$update_variables(new_data)
        self$project(cell_types)
        new_data <- self$get_data()
        filtering_iteration <-  filtering_iteration + 1
        print(paste("removed", previous_count - new_count, "points"))
      }
      private$add_filtering_log_step(
        "mahalanobis_filter",
        paste(
          paste('features =', paste0(features, collapse=',')),
          paste0("n_sigma = ", n_sigma),
          paste0("iterations = ", filtering_iteration),
          sep = ", "
        )
      )
      }
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
    #' This is best starting point to run optimization
    #' This is how we run optimization while performed comparison with other methods. you can use this method as a template for yourself
    #' @param config optimization config (coef_hinge_H, coef_hinge_W, coef_der_X, coef_der_Omega) will be overwritten.
    default_optimization = function(
    config = OPTIM_CONFIG_DEFAULT
    ) {
      private$initialize_first()
      LR_DECAY_STEPS = 15
      PARAMETERS_INCREASE_STEPS = 5
      lr_decay <- 0.5
      params_increase <- 10
      original_lambda_term <- config$coef_hinge_H #coef_hinge_H
      original_beta_term <- config$coef_hinge_W #coef_hinge_W
      lr_x <- config$coef_der_X
      lr_omega <- config$coef_der_Omega
      RUNS_EACH_STEP <- 1000
      pb <- progress_bar$new(total = LR_DECAY_STEPS * PARAMETERS_INCREASE_STEPS)
      for (lr_step in 1:LR_DECAY_STEPS) {
        lambda_term <-  original_lambda_term * lr_x * lr_x
        beta_term <- original_beta_term   * lr_omega * lr_omega
            for (x in 1:PARAMETERS_INCREASE_STEPS) {
                # Main training method, you can just run this
                config$coef_hinge_H <- lambda_term
                config$coef_hinge_W <- beta_term
                config$coef_der_X <- lr_x
                config$coef_der_Omega <- lr_omega
                self$optim_solution(RUNS_EACH_STEP, config)
        lambda_term <- lambda_term * params_increase
        beta_term <- beta_term * params_increase
        pb$tick()
      }
      lr_x <- lr_x * lr_decay
      lr_omega <- lr_omega * lr_decay
    }
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
    get_V_row = function() {
      private$set_data_first()

      res <- sinkhorn_sweep_c(
        V = Biobase::exprs(self$get_data()),
        D_vs_row = self$st$scaling$D_vs_row,
        D_vs_col = self$st$scaling$D_vs_col,
        iter = self$st$scaling$iterations,
        do_last_step = 0
      )

      rownames(res) <- rownames(self$get_data())
      colnames(res) <- colnames(self$get_data())

      res
   },

    #' @description
    #' Calculate V_column on the fly.
    get_V_column = function() {
      private$set_data_first()

      res <- sinkhorn_sweep_c(
        V = Biobase::exprs(self$get_data()),
        D_vs_row = self$st$scaling$D_vs_row,
        D_vs_col = self$st$scaling$D_vs_col,
        iter = self$st$scaling$iterations,
        do_last_step = 1
      )

      rownames(res) <- rownames(self$get_data())
      colnames(res) <- colnames(self$get_data())

      res
    },

    #' @description
    #' Get coordinates from external W and H using extended sinkhorn procedure.
    #'
    #' @param W  first factorization matrix for original V (V = WH). Rownames should match dso$st$data.
    #' @param H  second factorization matrix for original V (V = WH). Colnames should match dso$st$data.
    get_coordinates_from_external_matrices = function(W, H) {
      private$project_first()
      extended_scaling_result <- extended_sinkhorn_scale(V = Biobase::exprs(self$st$data),
                                                         W=W[rownames(dso$st$data),],
                                                         H=H[, colnames(dso$st$data)],
                                                         n_iter = dso$st$scaling$iterations)
      H_ss <-  extended_scaling_result$H_row
      W_gs <-  extended_scaling_result$W_col
      res <- get_coordinates_from_scaled_matrices(H_ss = H_ss, W_gs=W_gs, proj=self$st$proj)
      return(res)
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
