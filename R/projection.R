############ MAIN LOGIC ############

calc_svd_ops <- function(V_row, dims = NULL) {
  if (is.null(dims)) {
    dims <- 1:min(dim(V_row))
  }
  svd_ <- svd(V_row)
  S <- t(svd_$u[, dims])
  R <- t(svd_$v[, dims])
  Sigma <- diag(svd_$d[dims])
  
  # TODO: does this ever really happen?
  if (all(R[1, ] < 0)) {
    S[1, ] <- -S[1, ]
    R[1, ] <- -R[1, ]
  }
  
  # TODO: extract this common code
  rownames(R) <- paste0("dim_", 1:nrow(R))
  colnames(R) <- colnames(V_row)
  rownames(S) <- rownames(R)
  colnames(S) <- rownames(V_row)
  
  A <- matrix(apply(R, 1, sum),
              ncol = 1,
              nrow = length(dims))
  B <- matrix(apply(S, 1, sum),
              ncol = 1,
              nrow = length(dims))
  
  rownames(A) <- rownames(R)
  rownames(B) <- rownames(S)
  
  
  return(list(
    S = S,
    R = R,
    Sigma = Sigma,
    A = A,
    B = B
  ))
}

svd_project_with_ops <- function(scaling, ops) {
  proj <- list(
    X = scaling$V_row %*% t(ops$R),
    Omega = t(ops$S %*% scaling$V_column),
    meta = ops
  )
  
  proj$meta$M <- nrow(proj$X)
  proj$meta$N <- nrow(proj$Omega)
  proj$meta$K <- ncol(proj$X)

  return(proj)
}

svd_project <- function(scaling, dims) {
  ops <- calc_svd_ops(scaling$V_row, dims)
  proj <- svd_project_with_ops(scaling, ops)
  return(proj)
}

add_proj_umap <- function(proj, with_model = FALSE, neighbors_X = 15, neighbors_Omega = 15) {
  nn_method <- if (with_model) "annoy" else "fnn"
  
  data_X <- proj$X[, 2:proj$meta$K]
  data_Omega <- proj$Omega[, 2:proj$meta$K]

  proj$umap <- list(
    X = uwot::umap(data_X, nn_method = nn_method, ret_model = with_model, n_neighbors = neighbors_X),
    Omega = uwot::umap(data_Omega, nn_method = nn_method, ret_model = with_model, n_neighbors = neighbors_Omega)
  )
  if (with_model) {
    proj$umap$model <- list(
      X = proj$umap$X,
      Omega = proj$umap$Omega
    )
    proj$umap$X <- proj$umap$X$embedding
    proj$umap$Omega <- proj$umap$Omega$embedding
  } else {
    proj$umap$model <- NULL
  }
  
  rownames(proj$umap$X) <- rownames(proj$X)
  rownames(proj$umap$Omega) <- rownames(proj$Omega)
  colnames(proj$umap$X) <- c("umap_1", "umap_2")
  colnames(proj$umap$Omega) <- c("umap_1", "umap_2")
  
  return(proj)
}

transform_proj_umap <- function(points, proj) {
  if (is.null(proj$umap)) {
    stop(paste(
      "Projection umap is not calculated.",
      "Use add_proj_umap(proj, with_model = TRUE) first."
    ))
  } else if (is.null(proj$umap$model)) {
    stop(paste(
      "No model was calculated for umap.",
      "Use add_proj_umap(proj, with_model = TRUE) first."
    ))
  }
  transformed <- points
  transformed$X <- uwot::umap_transform(points$X[, 2:proj$meta$K], model = proj$umap$model$X)
  transformed$Omega <- uwot::umap_transform(points$Omega[, 2:proj$meta$K], model = proj$umap$model$Omega, n_epochs = 100, learning_rate = 0.001)
  colnames(transformed$X) <- c("umap_1", "umap_2")
  colnames(transformed$Omega) <- c("umap_1", "umap_2")
  return(transformed)
}

reverse_svd_projection <- function(X_space_pts, Omega_space_pts, proj) {
  return(list(
    rownorm = X_space_pts %*% proj$meta$R,
    colnorm = t(proj$meta$S) %*% Omega_space_pts
  ))
}

# TODO: peeks into solution layer
reverse_solution_projection <- function(solution_proj, proj) {
  solution_scaled <- reverse_svd_projection(solution_proj$X, t(solution_proj$Omega), proj)
  names(solution_scaled) <- c("H_row", "W_col")
  if (is.null(rownames(solution_scaled$H_row)) && is.null(colnames(solution_scaled$W_col))) {
    rownames(solution_scaled$H_row) <- paste0("cell_type_", 1:nrow(solution_scaled$H_row))
    colnames(solution_scaled$W_col) <- rownames(solution_scaled$H_row)
  }
  return(solution_scaled)
}


############ PLOTTING ############
plot_proj_svd <- function(proj, dims = NULL) {
  svd_d <- diag(proj$meta$Sigma)
  vars <- svd_d^2
  df <- data.frame(
    dimension = 1:length(vars),
    variance = cumsum(vars / sum(vars))
  )
  if (is.null(dims)) {
    dims = 1:length(svd_d)
  }
  ggplot(data = df[dims,], aes(x = dimension, y = variance)) +
    geom_point(
      aes(y = variance),
      size = 0.5,
      alpha = 1
    ) +
    geom_line(
      aes(y = variance),
      size = 0.5,
      alpha = 1
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.line.x = element_line(
        colour = 'black',
        size = 0.5,
        linetype = 'solid'
      ),
      axis.line.y = element_line(
        colour = 'black',
        size = 0.5,
        linetype = 'solid'
      ),
      legend.position = "none"
    ) +
    scale_x_continuous(minor_breaks = dims, limits = c(min(dims), max(dims)))
}
