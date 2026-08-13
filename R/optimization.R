############ MAIN LOGIC ############

#' Make optimization config object for the training method
#'
#' Just to have centralized object to change
#' @param method method of optimization to use can be  basic/positivity.
#' @param debug_stats keep track on gradient norms
#' @param coef_der_X learning rate for X space
#' @param coef_der_Omega learning rate for Omega space
#' @param coef_hinge_H positiviy penalty for X space (lambda)
#' @param coef_hinge_W positiviy penalty for Omega space (beta)
#' @param coef_alignment alignment penalty (gamma)
#' @param total_regularization_weight regularization weight for optimization
#' @param reg_X regularization weight for X
#' @param reg_Omega regularization weight for Omega
#' @param convergence_tol the limit for the learnig rate decrasing
#' @param stop_criteria_window learning rate will be decreased if total error did not change within this window
#' @param x_center X rays around which to perform search in theta search.
#' @param omega_center Omega rays around which to perform search in theta search.
#' @param center_threshold constraint for the  step.
#' @param coef_pos_D_h EXPERIMENTAL: penalty for D_h value (how far X^TxDh if from A(R))). should be 0 since not tested.
#' @param coef_pos_D_w EXPERIMENTAL: penalty for D_w value (how far OmegaxDw if from B(S))). should be 0 since not tested.
#' @param limit_X EXPERIMENTAL: if you want to restrict X from changing to much. should be 0 since not tested.
#' @param limit_Omega EXPERIMENTAL: if you want to restrict Omega from changing to much. should be 0 since not tested.
#' @param cosine_thresh  EXPERIMENTAL: if you want to restrict derivative from changing to much. should be 0 since not tested.
#' @return ready to use list with algorithm configuration
#' @export
optim_config <- function(
  method = c("positivity", "coordinate_descent", "theta", "alignment"), # let R handel selection and checking
  debug_stats = FALSE,
  coef_der_X = 0.01,
  coef_der_Omega = 0.01,
  coef_hinge_H = 0.5,
  coef_hinge_W = 0.5,
  coef_alignment = 0.5, # [CJLee] first try default as 0.5
  # Positivity method with fair gradients and stopping criteria
  total_regularization_weight = 0,
  reg_X = 1,
  reg_Omega = 1,
  stop_criteria_window = 800,
  convergence_tol = 1e-9,
  patience = 10,
  decade_rate = 0.5,
  max_drop = 10,
  epsilon = 1e-8,
  # Theta optimization within angle
  x_center = NULL,
  omega_center = NULL,
  center_threshold = 0,
  # Experimental params
  limit_X = 0,
  limit_Omega = 0,
  cosine_thresh = 0,
  coef_pos_D_h = 0,
  coef_pos_D_w = 0
) {
  return(list(
    method = match.arg(method), # this will check if user provide proper method name, otherwise use the default `positivity`
    debug_stats = debug_stats,
    coef_der_X = coef_der_X,
    coef_der_Omega = coef_der_Omega,
    coef_hinge_H = coef_hinge_H,
    coef_hinge_W = coef_hinge_W,
    coef_alignment = coef_alignment,
    total_regularization_weight = total_regularization_weight,
    reg_X = reg_X,
    reg_Omega = reg_Omega,
    convergence_tol = convergence_tol,
    stop_criteria_window = stop_criteria_window,
    patience = patience,
    decade_rate = decade_rate,
    max_drop = max_drop,
    epsilon = epsilon,
    x_center = x_center,
    omega_center = omega_center,
    center_threshold = center_threshold,
    coef_pos_D_h = coef_pos_D_h,
    coef_pos_D_w = coef_pos_D_w,
    limit_X = limit_X,
    limit_Omega = limit_Omega,
    cosine_thresh = cosine_thresh
  ))
}

OPTIM_CONFIG_DEFAULT <- optim_config()


#' Main optimization entry point
#'
#' Will perform optimization with selected parameters
#'
#' @param proj dso$st$proj object with information about projections (projected points, vectors)
#' @param solution_proj dso$st$solution_proj object with current solution matrices
#' @param iterations number of iterations to run
#' @param config configuratio for algorithm (default is optim_config())
#' @param block_name name for this particular run in logs
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @export
optimize_solution <- function(
  # TODO: remove scaled data, this should work on proj level
  proj,
  solution_proj,
  iterations,
  config = OPTIM_CONFIG_DEFAULT,
  block_name = NULL
) {
  # Cleaning inputs
  if (!("X" %in% names(solution_proj)) && ("Omega" %in% names(solution_proj))) {
    spdl::error("Both X and Omega must be initialized first in solution_proj")
    stop("Both X and Omega must be initialized first in solution_proj")
  }

  n_cell_types <- proj$meta$K

  # Managing optimization history
  if (!"optim_history" %in% names(solution_proj)) {
    solution_proj$optim_history <- list(
      blocks_statistics = data.frame(matrix(0, nrow = 0, ncol = 13)),
      errors_statistics = NULL,
      points_statistics_X = NULL,
      points_statistics_Omega = NULL,
      points_statistics_Dw = NULL,
      points_statistics_X_dtilda = NULL,
      points_statistics_Omega_dtilda = NULL,
      points_statistics_X_dtilda_uncorrected = NULL,
      points_statistics_Omega_dtilda_uncorrected = NULL,
      history = NULL,
      metadata = NULL
    )
  }

  if (is.null(block_name)) {
    block_name <- paste0("block_", nrow(solution_proj$optim_history$blocks_statistics) + 1)
  }

  if (!is.null(solution_proj$optim_history$history)) {
    from_idx <- max(solution_proj$optim_history$history$iteration) + 2
  } else if (!is.null(solution_proj$optim_history$errors_statistics)) {
    from_idx <- nrow(solution_proj$optim_history$errors_statistics) + 1
  } else {
    from_idx <- 1
  }


  solution_proj$optim_history$blocks_statistics <- rbind(
    solution_proj$optim_history$blocks_statistics,
    c(
      block_name,
      from_idx,
      from_idx + iterations - 1,
      config$coef_der_X,
      config$coef_der_Omega,
      config$coef_hinge_H,
      config$coef_hinge_W,
      config$coef_pos_D_h,
      config$coef_pos_D_w,
      iterations,
      config$limit_X,
      config$limit_Omega,
      config$cosine_thresh
    )
  )

  colnames(solution_proj$optim_history$blocks_statistics) <- c(
    "block_name",
    "from",
    "to",
    "coef_der_X",
    "coef_der_Omega",
    "coef_hinge_H",
    "coef_hinge_W",
    "coef_pos_D_h",
    "coef_pos_D_w",
    "iterations",
    "limit_X",
    "limit_Omega",
    "cosine_thresh"
  )

  # Running optimization
  mean_radius_X <-
    mean(apply(proj$X[, -1], 1, function(x) {
      norm(x, "2")
    }))

  mean_radius_Omega <-
    mean(apply(proj$Omega[, -1], 1, function(x) {
      norm(x, "2")
    }))

  r_limits <- calc_r_limits(proj, config$limit_X, config$limit_Omega)

  optimization_params <- list(
    X = solution_proj$X,
    Omega = t(solution_proj$Omega),
    D_w = solution_proj$D_w,
    SVRt = proj$meta$Sigma, # Should be equal to SVRt
    R = proj$meta$R,
    S = proj$meta$S,
    coef_der_X = config$coef_der_X,
    coef_hinge_H = config$coef_hinge_H,
    coef_hinge_W = config$coef_hinge_W,
    cell_types = n_cell_types,
    N = proj$meta$N,
    M = proj$meta$M,
    iterations = iterations
  )
  optimization_result <- if (config$method == "positivity") {
    optimization_params$total_regularization_weight <- config$total_regularization_weight
    optimization_params$reg_X <- config$reg_X
    optimization_params$reg_Omega <- config$reg_Omega
    optimization_params$convergence_tol <- config$convergence_tol
    optimization_params$debug_stats <- config$debug_stats
    optimization_params$stop_criteria_window <- config$stop_criteria_window
    do.call(optimize_positivity, optimization_params)
  } else if (config$method == "coordinate_descent") {
    optimization_params$r_const_X <- r_limits$R_limit_X
    optimization_params$r_const_Omega <- r_limits$R_limit_Omega
    optimization_params$thresh <- config$cosine_thresh
    optimization_params$coef_pos_D_h <- config$coef_pos_D_h
    optimization_params$coef_pos_D_w <- config$coef_pos_D_w
    optimization_params$coef_der_Omega <- config$coef_der_Omega
    optimization_params$mean_radius_X <- mean_radius_X
    optimization_params$mean_radius_Omega <- mean_radius_Omega
    optimization_params$convergence_tol <- config$convergence_tol
    optimization_params$debug_stats <- config$debug_stats
    optimization_params$stop_criteria_window <- config$stop_criteria_window
    do.call(optimize_coordinate_descent, optimization_params)
  } else if (config$method == "theta") {
    optimization_params$r_const_X <- r_limits$R_limit_X
    optimization_params$r_const_Omega <- r_limits$R_limit_Omega
    optimization_params$thresh <- config$cosine_thresh
    optimization_params$coef_pos_D_h <- config$coef_pos_D_h
    optimization_params$coef_pos_D_w <- config$coef_pos_D_w
    optimization_params$coef_der_Omega <- config$coef_der_Omega
    optimization_params$mean_radius_X <- mean_radius_X
    optimization_params$mean_radius_Omega <- mean_radius_Omega

    # this optimization ensures that solution points are not going away to far from the predefined center points.
    # the distance is measured as cosine distance between rays originating from 0.
    optimization_params$X_center <- config$x_center # predefined center point for X space. could be NULL
    if (!is.null(config$omega_center)) {
      optimization_params$Omega_center <- t(config$omega_center) # predefined center point for Omega space. could be NULL
    } else {
      optimization_params$Omega_center <- config$omega_center
    }
    optimization_params$theta_threshold <- config$center_threshold # threshold for the angle
    do.call(optimize_theta, optimization_params)
  } else if (config$method == "alignment") {
    optimization_params$max_iteration <- optimization_params$iterations
    optimization_params$iterations <- NULL

    optimization_params$convergence_tol <- if (!is.null(config$convergence_tol)) config$convergence_tol else 1e-4
    optimization_params$patience <- if (!is.null(config$patience)) as.integer(config$patience) else 10L
    optimization_params$decade_rate <- if (!is.null(config$decade_rate)) config$decade_rate else 0.5
    optimization_params$max_drop <- if (!is.null(config$max_drop)) as.integer(config$max_drop) else 10L
    optimization_params$epsilon <- if (!is.null(config$epsilon)) config$epsilon else 1e-8

    optimization_params$debug_stats <- config$debug_stats

    # for compatibility, might have some overhead
    optimization_params[["initial_X"]] <- optimization_params$X
    optimization_params[["X"]] <- NULL
    optimization_params[["initial_Omega"]] <- optimization_params$Omega
    optimization_params[["Omega"]] <- NULL
    optimization_params[["initial_D_w"]] <- optimization_params$D_w
    optimization_params[["D_w"]] <- NULL
    optimization_params[["cell_types"]] <- NULL
    optimization_params[["M"]] <- NULL
    optimization_params[["N"]] <- NULL

    do.call(optimize_alignment_pgd, optimization_params)
  } else {
    spdl::warn("Unknown optimization method. Will do the basic one")
    optimization_params$r_const_X <- r_limits$R_limit_X
    optimization_params$r_const_Omega <- r_limits$R_limit_Omega
    optimization_params$thresh <- config$cosine_thresh
    optimization_params$coef_pos_D_h <- config$coef_pos_D_h
    optimization_params$coef_pos_D_w <- config$coef_pos_D_w
    optimization_params$coef_der_Omega <- config$coef_der_Omega
    optimization_params$convergence_tol <- config$convergence_tol
    optimization_params$debug_stats <- config$debug_stats
    optimization_params$stop_criteria_window <- config$stop_criteria_window
    do.call(optimize_coordinate_descent, optimization_params)
  }

  solution_proj$X <- optimization_result$new_X
  solution_proj$Omega <- t(optimization_result$new_Omega)
  solution_proj$D_w <- optimization_result$new_D_w
  solution_proj$D_h <- optimization_result$new_D_h

  colnames(solution_proj$Omega) <- rownames(proj$meta$R)
  colnames(solution_proj$X) <- rownames(proj$meta$R)

  target_iterations <- ifelse(from_idx == 1, iterations + 1, iterations)

  solution_proj$optim_history$errors_statistics <- rbind(
    solution_proj$optim_history$errors_statistics,
    tail(optimization_result$errors_statistics, target_iterations)
  )

  solution_proj$optim_history$points_statistics_X <- rbind(
    solution_proj$optim_history$points_statistics_X,
    tail(optimization_result$points_statistics_X, target_iterations)
  )
  solution_proj$optim_history$points_statistics_Omega <- rbind(
    solution_proj$optim_history$points_statistics_Omega,
    tail(optimization_result$points_statistics_Omega, target_iterations)
  )
  solution_proj$optim_history$points_statistics_Dw <- rbind(
    solution_proj$optim_history$points_statistics_Dw,
    tail(optimization_result$points_statistics_Dw, target_iterations)
  )
  if (!is.null(optimization_result$history)) {
    new_history <- optimization_result$history
    # Continue iteration numbers from previous run
    if (from_idx > 1 && !is.null(solution_proj$optim_history$history)) {
      # Drop the overlapping first iteration (iteration 0 of the new run)
      new_history <- new_history[new_history$iteration > 0, , drop = FALSE]

      # Offset the iteration numbers so they continue from the previous run
      last_iter <- max(solution_proj$optim_history$history$iteration)
      new_history$iteration <- new_history$iteration + last_iter
    }
    solution_proj$optim_history$history <- rbind(
      solution_proj$optim_history$history,
      new_history
    )
  }
  if (!is.null(optimization_result$metadata)) {
    solution_proj$optim_history$metadata <- optimization_result$metadata
  }
  if (config$method == "positivity") {
    #   solution_proj$optim_history$points_statistics_X_dtilda <- rbind(
    #     solution_proj$optim_history$points_statistics_X_dtilda,
    #     tail(optimization_result$points_statistics_X_dtilda, target_iterations)
    # )
    #   solution_proj$optim_history$points_statistics_X_dtilda_uncorrected <- rbind(
    #     solution_proj$optim_history$points_statistics_X_dtilda_uncorrected,
    #     tail(optimization_result$points_statistics_X_dtilda_uncorrected, target_iterations)
    # )
    #   solution_proj$optim_history$points_statistics_Omega_dtilda <- rbind(
    #     solution_proj$optim_history$points_statistics_Omega_dtilda,
    #     tail(optimization_result$points_statistics_Omega_dtilda, target_iterations)
    # )
    #   solution_proj$optim_history$points_statistics_Omega_dtilda_uncorrected <- rbind(
    #     solution_proj$optim_history$points_statistics_Omega_dtilda_uncorrected,
    #     tail(optimization_result$points_statistics_Omega_dtilda_uncorrected, target_iterations)
    # )
  }

  if (!is.null(solution_proj$optim_history$errors_statistics)) { # for compatibility
    colnames(solution_proj$optim_history$errors_statistics) <- head(
      c(
        "deconv_error",
        "lambda_error",
        "beta_error",
        "D_h_error",
        "D_w_error",
        "total_error",
        "scaled_total_error",
        "neg_props_count",
        "neg_basis_count",
        "sum_d_w",
        "average_norm",
        "learning_rate",
        "gradient_norm",
        "average_hinge_H_gradient_norm",
        "average_hinge_W_gradient_norm",
        "average_reg_X_gradient_norm",
        "average_reg_Omega_gradient_norm",
        "best_error_value",
        "best_error_iteration",
        "scaled_lambda_error",
        "scaled_beta_error",
        "average_final_gradient_norm",
        "total_shrink_iterations"
      ),
      ncol(solution_proj$optim_history$errors_statistics)
    )
  }

  return(solution_proj)
}

calc_r_limits <- function(
  data_proj,
  limit_X = 0,
  limit_Omega = 0
) {
  if (limit_X > 0) {
    zero_distances <- calc_partial_dist(data_proj$X, with_dims = c(2, ncol(data_proj$X)))
    limit_num_X <- floor(nrow(data_proj$X) * limit_X)
    R_limit_X <- norm(data_proj$X[names(zero_distances[limit_num_X]), -1], "2")
  } else {
    R_limit_X <- 0
  }

  if (limit_Omega > 0) {
    zero_distances <- calc_partial_dist(data_proj$Omega, with_dims = c(2, ncol(data_proj$Omega)))
    limit_num_Omega <- floor(nrow(data_proj$Omega) * limit_Omega)
    R_limit_Omega <- norm(data_proj$Omega[names(zero_distances[limit_num_Omega]), -1], "2")
  } else {
    R_limit_Omega <- 0
  }

  return(list(
    R_limit_X = R_limit_X,
    R_limit_Omega = R_limit_Omega
  ))
}


############ PLOTTING ############

#' Plot errors after optimization
#'
#' Will plot selected terms of the cost function
#'
#' @param solution_proj dso$st$solution_proj object with current solution matrices
#' @param variables which variables to plot (columns of solution_proj$optim_history$errors_statistics)
#' @return ggplot object
#' @import reshape2
#' @export
plot_errors <- function(
  solution_proj,
  variables = c(
    "deconv_error",
    "lambda_error",
    "beta_error",
    "D_h_error",
    "D_w_error",
    "total_error"
  )
) {
  # for compatibility
  if (is.null(solution_proj$optim_history$history)) { # for compatibility
    to_plot <- data.frame(solution_proj$optim_history$errors_statistics[, variables, drop = F])
    if (nrow(solution_proj$optim_history$errors_statistics) == 0) {
      spdl::warn("Nothing to plot, errors_statistics was empty")
    } else {
      to_plot$iteration <- 0:(nrow(solution_proj$optim_history$errors_statistics) - 1)
      to_plot <- reshape2::melt(to_plot, id.vars = "iteration", measure.vars = variables)
      plt <-
        ggplot(to_plot, aes(
          x = .data$iteration,
          y = log10(.data$value),
          color = .data$variable
        )) +
        geom_line() +
        theme_minimal()
      return(plt) # early return so we can skip else clause
    }
  }
  history <- solution_proj$optim_history$history
  ggplot(
    history[history$metric %in% variables, ],
    mapping = aes(
      x = iteration,
      y = log10(value),
      color = metric
    )
  ) +
    geom_line() +
    theme_minimal()
}


#' Plot changes in negativity of matrix H
#'
#' Will plot percentage of negative elements for current solution.
#'
#' @param proj dso$st$proj object with info about projection. (used to get N, M, K values) TODO: remove this.
#' @param solution_proj dso$st$solution_proj object with current solution matrices in projected space
#' @return ggplot object
#' @export
plot_negative_proportions_change <- function(proj, solution_proj) {
  N <- proj$meta$N
  M <- proj$meta$M
  K <- proj$meta$K
  total_H <- N * K

  if (is.null(solution_proj$optim_history$history)) { # for compatibility
    errors_statistics <- solution_proj$optim_history$errors_statistics
    toPlot <- as.data.frame(errors_statistics[, "neg_props_count", drop = F])
    last_prop_count <- round(toPlot[nrow(toPlot), "neg_props_count"] / total_H, 6) * 100
    toPlot$iteration <- 0:(nrow(toPlot) - 1)
    plt <- ggplot(toPlot, aes(y = .data$neg_props_count, x = .data$iteration)) +
      geom_line() +
      theme_minimal() +
      xlab("Iteration") +
      ylab("Negative proportions") +
      annotate("text", x = Inf, y = Inf, label = paste0(last_prop_count, "%"), vjust = 1, hjust = 1) +
      ggtitle("Number of negative proportions")
    return(plt)
  }

  history <- solution_proj$optim_history$history
  last_prop_percentage <- round(tail(history[history$metric == "neg_props", "value"], 1) / total_H, 6) * 100

  ggplot(
    history[history$metric == "neg_props", ],
    mapping = aes(
      x = iteration,
      y = value
    )
  ) +
    geom_line() +
    theme_minimal() +
    xlab("Iteration") +
    ylab("Negative proportions") +
    annotate("text", x = Inf, y = Inf, label = paste0(last_prop_percentage, "%"), vjust = 1, hjust = 1) +
    ggtitle("Number of negative proportions")
}

#' Plot changes in negativity of matrix W
#'
#' Will plot percentage of negative elements for current solution.
#'
#' @param proj dso$st$proj object with info about projection. (used to get N, M, K values) TODO: remove this.
#' @param solution_proj dso$st$solution_proj object with current solution matrices in projected space
#' @return ggplot object
#' @export
plot_negative_basis_change <- function(proj, solution_proj) {
  N <- proj$meta$N
  M <- proj$meta$M
  K <- proj$meta$K
  total_W <- M * K

  if (is.null(solution_proj$optim_history$history)) { # for compatibility
    errors_statistics <- solution_proj$optim_history$errors_statistics
    toPlot <- as.data.frame(errors_statistics[, "neg_basis_count", drop = F])
    last_basis_count <- round(toPlot[nrow(toPlot), "neg_basis_count"] / total_W, 6) * 100
    toPlot$iteration <- 0:(nrow(toPlot) - 1)
    plt <- ggplot(toPlot, aes(y = .data$neg_basis_count, x = .data$iteration)) +
      geom_line() +
      theme_minimal() +
      xlab("Iteration") +
      ylab("Negative basis") +
      annotate("text", x = Inf, y = Inf, label = paste0(last_basis_count, "%"), vjust = 1, hjust = 1) +
      ggtitle("Number of negative basis elements")
    return(plt)
  }

  history <- solution_proj$optim_history$history
  neg_basis_percentage <- round(tail(history[history$metric == "neg_basis", "value"], 1) / total_W, 6) * 100

  ggplot(
    history[history$metric == "neg_basis", ],
    mapping = aes(
      x = iteration,
      y = value
    )
  ) +
    geom_line() +
    theme_minimal() +
    xlab("Iteration") +
    ylab("Negative basis") +
    annotate("text", x = Inf, y = Inf, label = paste0(neg_basis_percentage, "%"), vjust = 1, hjust = 1) +
    ggtitle("Number of negative basis elements")
}
