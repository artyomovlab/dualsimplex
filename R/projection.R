############ MAIN LOGIC ############

#' Calculate SVD for the matrix
#'
#' produces parts of dso$st$proj object, containing vectors R,S matrix Sigma, vectors A(R) and B(S)
#' performs SVD on the matrix
#'
#' @param V_row input data matrix to perform SVD (in this method it should be Sinkhorn transformed matrix)
#' @param dims how many dimensions of SVD we leave
#' @return list with all calculated matrices
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

#' Project sinkhorn scaled matrix to SVD space
#'
#' produces parts dso$st$proj object,  ($meta and  projected points $X and $Omega)
#'
#' @param scaling dso$st$scaling object containing sinkhorn scaling result
#' @param ops svd result for the matirx
#' @return proj object
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

#' Prepare projection object
#'
#' Entry point to produce dso$st$proj object, containing all information about projection as well as projected points
#'
#' @param scaling dso$st$scaling object containing sinkhorn scaling result
#' @param dims how many dimensions we want to get
#' @return proj object
#' @export
svd_project <- function(scaling, dims) {
  ops <- calc_svd_ops(scaling$V_row, dims)
  proj <- svd_project_with_ops(scaling, ops)
  return(proj)
}

#' Adds UMAP representation to dso$st$proj object
#'
#' calls uwot::umap to make umap
#'
#' @param proj dso$st$proj object containing projected points and info about projection
#' @param with_model if TRUE uses "annoy" model else "fnn"
#' @param neighbors_X param for UMAP, play to get better picture
#' @param neighbors_Omega param for UMAP, play to get better picture
#' @return modified proj object
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

#' Transform additional poitns to the same UMAP
#'
#' @param points new points we want to map to same umap
#' @param proj dso$st$proj object containing projected points and umap info
#' @return transformed points
#' @export
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
  transformed$X <- uwot::umap_transform(points$X[, 2:proj$meta$K], model = proj$umap$model$X, n_epochs = 100, learning_rate = 0.001)
  transformed$Omega <- uwot::umap_transform(points$Omega[, 2:proj$meta$K], model = proj$umap$model$Omega, n_epochs = 100, learning_rate = 0.001)
  colnames(transformed$X) <- c("umap_1", "umap_2")
  colnames(transformed$Omega) <- c("umap_1", "umap_2")
  return(transformed)
}

#' Unproject points. Go from projected X coordinates to original row normalized matrix H_ss and
#' from projected Omega points to original column normalized matrix W_gs
#'
#' @param X_space_pts points in samples space (X space, samples space, left simplex)
#' @param Omega_space_pts   points in features space (Omega space, genes space, right simplex)
#' @param proj dso$st$proj object containing projected points and umap info
#' @return list of two matrices (H_ss, W_gs)
#' @export
reverse_svd_projection <- function(X_space_pts, Omega_space_pts, proj) {
  return(list(
    rownorm = X_space_pts %*% proj$meta$R,
    colnorm = t(proj$meta$S) %*% Omega_space_pts
  ))
}

# TODO: peeks into solution layer
#' Transform solution back from projected space to original space
#'
#' @param solution_proj dso$st$solution_proj object containing solution and optimization history
#' @param proj dso$st$proj object containing projected points and info about projection
#' @return solution_scaled object containing two matrices (H_ss, W_gs)
#' @export
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

# Util function to start plotting svd
plot_proj_svd <- function(proj, dims = NULL) {
  svd_d <- diag(proj$meta$Sigma)
  return(plot_svd_d(svd_d, dims))
}
# Util function to plot svd
plot_svd_d <- function(svd_d, dims = NULL, cumulative = T, variance = T) {
  vars <- svd_d
  if (variance) {
    vars <- vars ^ 2
  }
  df <- data.frame(
    dimension = 1:length(vars),
    variance = if (cumulative) cumsum(vars / sum(vars)) else vars
  )
  if (is.null(dims)) {
    dims = 1:length(svd_d)
  }
  return(ggplot(data = df[dims,], aes(x = dimension, y = variance)) +
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
    scale_x_continuous(minor_breaks = dims, limits = c(min(dims), max(dims))))
}

#' Util function to plot svd
#'
#' @param svd_ds svd Sigma (D) matrix to be plotted
#' @param cumulative TRUE if cumulative plot
#' @param variance TRUE if explained variance should be plotted
#' @return ggplot object
#' @importFrom tidyr pivot_longer
#' @export
plot_svd_ds_matrix <- function(svd_ds, cumulative = T, variance = T) {
  vars <- svd_ds
  if (variance) {
    vars <- vars ^ 2
  }
  if (cumulative) {
    vars <- t(apply(vars, 1, function(step_vars) {cumsum(step_vars / sum(step_vars))}))
  }
  vars <- as.data.frame(vars)
  colnames(vars) <- 1:ncol(vars)
  values <- colnames(vars)
  if (is.null(rownames(svd_ds)))
    rownames(svd_ds) <- 1:nrow(svd_ds)
  vars[, "step"] <- as.factor(rownames(svd_ds))
  to_plot <- tidyr::pivot_longer(
    vars,
    values,
    names_to = "component",
    values_to = "explained_variance"
  )
  to_plot$component <- as.integer(to_plot$component)
  
  return(ggplot(to_plot, aes(x = component, y = explained_variance, col = step)) + geom_line())
}