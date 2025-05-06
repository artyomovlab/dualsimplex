############ MAIN LOGIC ############

#' Make optimization config object for the training method
#'
#' Just to have centralized object to change
#'
#' @param coef_der_X learning rate for X space
#' @param coef_der_Omega learning rate for Omega space
#' @param coef_hinge_H positiviy penalty for X space (lambda)
#' @param coef_hinge_W positiviy penalty for Omega space (beta)
#' @param coef_pos_D_h EXPERIMENTAL: penalty for D_h value (how far X^TxDh if from A(R))). should be 0 since not tested.
#' @param coef_pos_D_w EXPERIMENTAL: penalty for D_w value (how far OmegaxDw if from B(S))). should be 0 since not tested.
#' @param limit_X EXPERIMENTAL: if you want to restrict X from changing to much. should be 0 since not tested.
#' @param limit_Omega EXPERIMENTAL: if you want to restrict Omega from changing to much. should be 0 since not tested.
#' @param cosine_thresh  EXPERIMENTAL: if you want to restrict derivative from changing to much. should be 0 since not tested.
#' @param method method of optimization to use can be  basic/positivity.
#' @return ready to use list with algorithm configuration
#' @export
optim_config <- function(
  coef_der_X = 0.001,
  coef_der_Omega = 0.001,
  coef_hinge_H = 10,
  coef_hinge_W = 1,
  coef_pos_D_h = 0,
  coef_pos_D_w = 0,
  limit_X = 0,
  limit_Omega = 0,
  cosine_thresh = 0,
  method = "basic" # basic/positivity
) {
  return(list(
    coef_der_X = coef_der_X,
    coef_der_Omega = coef_der_Omega,
    coef_hinge_H = coef_hinge_H,
    coef_hinge_W = coef_hinge_W,
    coef_pos_D_h = coef_pos_D_h,
    coef_pos_D_w = coef_pos_D_w,
    limit_X = limit_X,
    limit_Omega = limit_Omega,
    cosine_thresh = cosine_thresh,
    method = method
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
    stop("Both X and Omega must be initialized first in solution_proj")
  }

  n_cell_types <- proj$meta$K

  # Managing optimization history
  if (!"optim_history" %in% names(solution_proj)) {
    solution_proj$optim_history <- list(
      blocks_statistics = data.frame(matrix(0, nrow = 0, ncol = 13)),
      errors_statistics = NULL,
      points_statistics_X = NULL,
      points_statistics_Omega = NULL
    )
  }

  if (is.null(block_name)) {
    block_name <- paste0("block_", nrow(solution_proj$optim_history$blocks_statistics) + 1)
  }

  if (is.null(solution_proj$optim_history$errors_statistics)) {
    from_idx <- 1
  } else {
    from_idx <- nrow(solution_proj$optim_history$errors_statistics) + 1
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
    solution_proj$X,
    t(solution_proj$Omega),
    solution_proj$D_w,
    proj$meta$Sigma, # Should be equal to SVRt
    proj$meta$R,
    proj$meta$S,
    config$coef_der_X,
    config$coef_der_Omega,
    config$coef_hinge_H,
    config$coef_hinge_W,
    config$coef_pos_D_h,
    config$coef_pos_D_w,
    n_cell_types,
    proj$meta$N,
    proj$meta$M,
    iterations,
    mean_radius_X,
    mean_radius_Omega,
    r_limits$R_limit_X,
    r_limits$R_limit_Omega,
    config$cosine_thresh
  )
  optimization_result <- if (config$method == "positivity") {
    do.call(alternative_derivative_stage2, optimization_params)
  } else if (config$method == "basic") {
    do.call(derivative_stage2, optimization_params)
  } else {
    print("Unknown optimization method. Will do the basic one")
    do.call(derivative_stage2, optimization_params)
  }

  solution_proj$X <- optimization_result$new_X
  solution_proj$Omega <- t(optimization_result$new_Omega)
  solution_proj$D_w <- optimization_result$new_D_w
  solution_proj$D_h <- optimization_result$new_D_h

  colnames(solution_proj$Omega) <- rownames(proj$meta$R)
  colnames(solution_proj$X) <- rownames(proj$meta$R)

  solution_proj$optim_history$errors_statistics <- rbind(
    solution_proj$optim_history$errors_statistics,
      optimization_result$errors_statistics
  )

  solution_proj$optim_history$points_statistics_X <- rbind(
    solution_proj$optim_history$points_statistics_X,
      optimization_result$points_statistics_X
  )
  solution_proj$optim_history$points_statistics_Omega <- rbind(
    solution_proj$optim_history$points_statistics_Omega,
           optimization_result$points_statistics_Omega
  )

  colnames(solution_proj$optim_history$errors_statistics) <-
    c(
      "deconv_error",
      "lamdba_error",
      "beta_error",
      "D_h_error",
      "D_w_error",
      "total_error",
      "neg_props_count",
      "neg_basis_count",
      "sum_d_w"
    )
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
    "lamdba_error",
    "beta_error",
    "D_h_error",
    "D_w_error",
    "total_error"
  )
) {
  to_plot <- data.frame(solution_proj$optim_history$errors_statistics[, variables, drop=F])
  to_plot$iteration <- 0:(nrow(solution_proj$optim_history$errors_statistics) - 1)
  to_plot <- reshape2::melt(to_plot, id.vars = "iteration", measure.vars = variables)
  plt <-
    ggplot(to_plot, aes(
      x = .data$iteration,
      y = log10(.data$value),
      color = .data$variable
    )) +
    geom_line() + theme_minimal()
  plt
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
  errors_statistics <- solution_proj$optim_history$errors_statistics
  total_H <- N * K
  toPlot <- as.data.frame(errors_statistics[, "neg_props_count",drop=F])
  last_prop_count <- round(toPlot[nrow(toPlot), "neg_props_count"] / total_H,6) * 100
  toPlot$iteration <- 0:(nrow(toPlot) - 1)
  plt <- ggplot(toPlot,aes(y=.data$neg_props_count,x=.data$iteration)) + geom_line() + theme_minimal() +
      xlab("Iteration") + ylab("Negative proportions")+
    annotate("text",  x=Inf, y = Inf, label = paste0(last_prop_count,"%"), vjust=1, hjust=1) +
    ggtitle("Number of negative proportions")
  return(plt)
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
  errors_statistics <- solution_proj$optim_history$errors_statistics
  total_W <- M * K
  toPlot <- as.data.frame(errors_statistics[,"neg_basis_count",drop=F])
  last_basis_count <- round(toPlot[nrow(toPlot),"neg_basis_count"] / total_W,6) * 100
  toPlot$iteration <- 0:(nrow(toPlot) - 1)
  plt <- ggplot(toPlot,aes(y=.data$neg_basis_count,x=.data$iteration)) + geom_line() + theme_minimal() +
      xlab("Iteration") + ylab("Negative basis")+
      annotate("text",  x=Inf, y = Inf, label = paste0(last_basis_count,"%"), vjust=1, hjust=1) +
    ggtitle("Number of negative basis elements")
  return(plt)
}


