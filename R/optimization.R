############ MAIN LOGIC ############

optim_config <- function(
  coef_der_X = 0.001,
  coef_der_Omega = 0.001,
  coef_hinge_H = 10,
  coef_hinge_W = 1,
  coef_pos_D_h = 0,
  coef_pos_D_w = 0,
  limit_X = 0,
  limit_Omega = 0,
  cosine_thresh = 0
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
    cosine_thresh = cosine_thresh
  ))
}

OPTIM_CONFIG_DEFAULT <- optim_config()

optimize_solution <- function(
  # TODO: remove scaled data, this should work on proj level
  scaling,
  proj,
  solution_proj,
  iterations,
  config = OPTIM_CONFIG_DEFAULT,
  block_name = NULL
) {
  # Cleaning inputs
  if (!("X" %in% solution_proj) && ("Omega" %in% solution_proj)) {
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

  step_errors_statistics <- matrix(0, nrow = iterations, ncol = 10)
  step_points_statistics_X <- matrix(0, nrow = iterations, ncol = n_cell_types ^ 2)
  step_points_statistics_Omega <- matrix(0, nrow = iterations, ncol = n_cell_types ^ 2)

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

  res_ <- derivative_stage2(
    solution_proj$X,
    t(solution_proj$Omega),
    solution_proj$D_w,
    scaling$V_row,
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
    step_errors_statistics,
    0,
    step_points_statistics_X,
    step_points_statistics_Omega,
    mean_radius_X,
    mean_radius_Omega,
    r_limits$R_limit_X,
    r_limits$R_limit_Omega,
    config$cosine_thresh
  )

  solution_proj$X <- res_[[1]]
  solution_proj$Omega <- t(res_[[2]])
  solution_proj$D_w <- res_[[3]]
  solution_proj$D_h <- res_[[4]]

  colnames(solution_proj$Omega) <- rownames(proj$meta$R)
  colnames(solution_proj$X) <- rownames(proj$meta$R)

  solution_proj$optim_history$errors_statistics <- rbind(
    solution_proj$optim_history$errors_statistics,
    res_[[5]]
  )
  solution_proj$optim_history$points_statistics_X <- rbind(
    solution_proj$optim_history$points_statistics_X,
    res_[[6]]
  )
  solution_proj$optim_history$points_statistics_Omega <- rbind(
    solution_proj$optim_history$points_statistics_Omega,
    res_[[7]]
  )

  colnames(solution_proj$optim_history$errors_statistics) <-
    c(
      "deconv_error",
      "lamdba_error",
      "beta_error",
      "D_h_error",
      "D_w_error",
      "total_error",
      "orig_deconv_error",
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
  to_plot <- data.frame(solution_proj$optim_history$errors_statistics[, variables])
  to_plot$iteration <- 0:(nrow(solution_proj$optim_history$errors_statistics) - 1)
  to_plot <- reshape2::melt(to_plot, id.vars = "iteration", measure.vars = variables)
  plt <-
    ggplot(to_plot, aes(
      x = iteration,
      y = log10(value),
      color = variable
    )) +
    geom_line() + theme_minimal()
  plt
}
