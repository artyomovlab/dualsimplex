#' Create initial data simulation for gene expression
#'
#' Simulation builder initialization: simulation of basis, proportions and gene expression data without noise.
#' The motivation behind this builder approach is that every add... method creates a copy of the input simulation,
#' so that it is possible to create a single base simulation with many different child simulations.
#' For example, with different noise applied to the same data.
#' The add... functions inherently support magrittr's pipe operator (%>%).
#'
#' @importFrom Biobase exprs ExpressionSet fData pData
#' @param n_genes number of rows in generated matrix (M)
#' @param n_samples number of columns in generated matrix (N)
#' @param n_cell_types number of hidden main components  in generated matrix (K)
#' @param with_marker_genes if TRUE will also generate marker genes
#' @param mean_coefs  generation param
#' @param lbl_dataset TRUE for specific generation scenario of lbl dataset
#' @return sim object
#' @examples
#' require(dplyr)
#' sim_initial <- simulation_gene_expression(12000, 100, 3)
#' sim_noise_3 <- add_noise(sim_initial, 3)
#' sim_noise_5 <- sim_initial %>% add_noise(5) %>% add_basis_samples()
#' @import dplyr
#' @export
simulation_gene_expression <- function(
  n_genes,
  n_samples,
  n_cell_types,
  with_marker_genes = 0,
  mean_coefs = NULL,
  lbl_dataset = NULL
) {
  if (with_marker_genes != 0) {
    n_genes <- n_genes - n_cell_types * as.integer(with_marker_genes)
  }
  if (is.null(lbl_dataset)) {
    basis <- generate_gene_expression_basis(n_genes, n_cell_types)
  } else {
    basis <- generate_basis_lbl(n_genes, n_cell_types, lbl_dataset)
  }
  if (with_marker_genes != 0) {
    tmp <- append_marker_genes_to_basis(basis, as.integer(with_marker_genes))
    basis <- tmp[["basis"]]
    marker_gene_names <- tmp[["marker_gene_names"]]
  }
  
  coefs <- generate_proportions(n_samples, n_cell_types, mean_coefs)
  data <- basis %*% coefs
  
  data[data < 0] <- 0
  
  res <- list(
    basis = basis,
    coefs = coefs,
    data = data,
    mixed_sample_names = colnames(data)
  )
  if (with_marker_genes) {
    res$marker_gene_names <- marker_gene_names
  }
  return(res)
}

#' Create initial data simulation for uniform basis
#'
#' Simulation builder initialization: simulation of basis, proportions and gene expression data without noise.
#' The motivation behind this builder approach is that every add... method creates a copy of the input simulation,
#' so that it is possible to create a single base simulation with many different child simulations.
#' For example, with different noise applied to the same data.
#' The add... functions inherently support magrittr's pipe operator (%>%).
#' @param m number of rows in generated matrix (M)
#' @param n number of columns in generated matrix (N)
#' @param k number of hidden main components  in generated matrix (K)
#' @param mean_coefs  generation param
#' @param coefs_type "simplex"/"random" wether to generate H matrix from simplex or randomly
#' @return sim object
#' @examples
#' require(dplyr)
#' sim_initial <- simulation_uniform_basis(12000, 100, 3)
#' sim_noise_3 <- add_noise(sim_initial, 3)
#' sim_noise_5 <- sim_initial %>% add_noise(5) %>% add_basis_samples()
#' @import dplyr
#' @export
simulation_uniform_basis <- function(m, n, k, mean_coefs = NULL, coefs_type="simplex") {
  basis <- matrix(stats::runif(k * m, min = 0, max = 10), ncol = k)
  rownames(basis) <- paste("feature_", seq_len(m))
  colnames(basis) <- paste("component_", seq_len(k))
  if (coefs_type == "simplex") {
    coefs <- generate_proportions(n, k, mean_coefs)
  } else if (coefs_type == "random") {
    coefs <- matrix(stats::runif(k * n), nrow = k)
  } else {
    stop("uknown coefs type")
  }
  colnames(coefs) <- paste0("mixture_", seq_len(n))
  rownames(coefs) <- paste0("component_", seq_len(k))
  data <- basis %*% coefs
  data[data < 0] <- 0
  res <- list(
    basis = basis,
    coefs = coefs,
    data = data,
    mixed_sample_names = colnames(data)
  )
  return(res)
}

#' Create initial data simulation for outside basis matrix
#'
#' Simulation builder initialization: proportions and gene expression data without noise for predefined basis.
#' The motivation behind this builder approach is that every add... method creates a copy of the input simulation,
#' so that it is possible to create a single base simulation with many different child simulations.
#' For example, with different noise applied to the same data.
#' The add... functions inherently support magrittr's pipe operator (%>%).
#' @param predefined_basis predefined matrix W (for example flattened pictures)
#' @param n number of columns in generated matrix (N)
#' @param k will use first k columns of predefined basis as basis elements.
#' @param mean_coefs  generation param
#' @param coefs_type "simplex"/"random" wether to generate H matrix from simplex or randomly
#' @return sim object
#' @examples
#' require(dplyr)
#' basis <- matrix(stats::runif(12000 * 3, min = 0, max = 10), ncol = 3)
#' sim_initial <- simulation_mixtures_of_basis(basis, 100)
#' sim_noise_3 <- add_noise(sim_initial, 3)
#' sim_noise_5 <- sim_initial %>% add_noise(5) %>% add_basis_samples()
#' @import dplyr
#' @export
simulation_mixtures_of_basis <- function(predefined_basis, n, k = NULL, mean_coefs = NULL, coefs_type="simplex") {
  n_features <- nrow(predefined_basis)
  n_components <-  if (is.null(k)) ncol(predefined_basis) else k
  basis <- predefined_basis[, seq_len(n_components)]

  if (is.null(rownames(predefined_basis))) {
    rownames(basis) <- paste("feature_", seq_len(nrow(basis)))
  }
  if (is.null(colnames(basis))) {
    colnames(basis) <- paste("component_", seq_len(ncol(basis)))
  }
  if (coefs_type == "simplex") {
    coefs <- generate_proportions(n, n_components, mean_coefs)
  } else if (coefs_type == "random") {
    coefs <- matrix(stats::runif(n_components * n), nrow = n_components)
  } else {
    stop("uknown coefs type")
  }
  colnames(coefs) <- paste0("mixture_", seq_len(n))
  rownames(coefs) <- paste0("component_", seq_len(n_components))
  data <- basis %*% coefs
  data[data < 0] <- 0
  
  res <- list(
    basis = predefined_basis,
    coefs = coefs,
    data = data,
    mixed_sample_names = colnames(data)
  )
  return(res)
}

#' Add N pure samples (basis + noise) per cell type to the simulation
#'
#' @param simulation original simulation object
#' @param samples_per_cell_type number of pure samples to add for each hidden main component
#' @param noise_deviation noise deviation value for pure samples
#' @return modified sim object
#' @export
add_pure_samples <- function(simulation, samples_per_cell_type, noise_deviation = 0.05) {
  pure_samples <- generate_pure_samples(
    simulation$basis,
    samples_per_cell_type,
    noise_deviation = noise_deviation
  )
  simulation$pure_sample_names <- colnames(pure_samples$data)
  simulation$data <- cbind(simulation$data, pure_samples$data)
  simulation$coefs <- cbind(simulation$coefs, pure_samples$coefs)
  simulation$pure_samples_noise_deviation <- noise_deviation
  simulation
}

#' Add samples, which are pure basis (without noise) to the simulation
#'
#' @param simulation original simulation object
#' @return modified sim object
#' @export
add_basis_samples <- function(simulation) {
  basis_samples <- generate_basis_samples(simulation$basis)
  simulation$basis_sample_names <- colnames(basis_samples$data)
  simulation$data <- cbind(simulation$data, basis_samples$data)
  simulation$coefs <- cbind(simulation$coefs, basis_samples$coefs)
  simulation
}

#' Add noise to the simulation's data
#'
#' @param simulation original simulation object
#' @param noise_deviation noise deviation
#' @param protect_genes do not modigy these rows
#' @param protect_samples do not modigy these samples
#' @param noise_type norm_exp/additive/proportional controls how noise is generated
#' @return modified sim object
#' @export
add_noise <- function(simulation, noise_deviation, protect_genes = c(), protect_samples = c(), noise_type="norm_exp") {
  simulation$data <- apply_noise_to_data(simulation$data, noise_deviation, protect_genes, protect_samples, noise_type)
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  simulation$data_noise_deviation <- noise_deviation
  simulation
}


# All underlying logic is below
generate_gene_expression_basis <- function(n_genes, n_cell_types, sd_ = 0.2) {
  data_1 <- c(
    stats::rnorm(n_genes * 2, mean = 4, sd = 0.75),
    stats::rnorm(n_genes * 3, mean = 10, sd = 1.5)
  )
  basis <- matrix(0, nrow = n_genes, ncol = n_cell_types)
  sds_ <- stats::rnorm(n_cell_types, mean = sd_, sd = 0.02)
  
  for (i in 1:n_genes) {
    basis[i, 1] <- data_1[sample(seq_len(length(data_1)), 1)] * stats::rnorm(1, mean = 1, sd = sds_[1])
    for (j in 2:n_cell_types) {
      basis[i, j] <- basis[i, 1] * stats::rnorm(1, mean = 1, sd = sds_[j])
    }
  }
  
  
  # Common part with generate_basis_lbl below
  basis <- limma::normalizeBetweenArrays(basis, method = "quantile")
  basis <- 2^basis
  
  rownames(basis) <- paste0("gene_", 1:n_genes)
  colnames(basis) <- paste0("cell_type_", 1:n_cell_types)
  
  basis
}

generate_basis_lbl <- function(n_genes, n_cell_types, lbl_dataset) {
  if (n_cell_types != 3) {
    spdl::error("LBL basis works only for 3 cell types")
    stop("LBL basis works only for 3 cell types")

  }
  basis <- lbl_dataset[sample(nrow(lbl_dataset), n_genes), c(1, 4, 7)]
  
  # Common part with generate_basis below
  basis <- limma::normalizeBetweenArrays(basis, method="quantile")
  basis <- 2^basis
  
  basis
}


append_marker_genes_to_basis <- function(basis, n_each_type = 1) {
  n_cell_types <- ncol(basis)
  marker_expressions <- do.call(
    rbind,
    replicate(n_each_type, diag(1, n_cell_types, n_cell_types), simplify = FALSE)
  )
  rownames(marker_expressions) <- paste0("marker_gene_", seq_len(n_cell_types * n_each_type))
  list(
    basis = rbind(basis, marker_expressions),
    marker_gene_names = rownames(marker_expressions)
  )
}

#' Generate n points, uniformly distributed on a k-dimensional standard simplex
#' @param n number of points
#' @param k dimensionality of generatred data
#' @param M number of generations
sample_from_simplex_uniformly <- function(n, k, M = 100000) {
  X <- matrix(0, nrow = k + 1, ncol = n)
  X[k + 1, ] <- M
  
  X[2:k, ] <- replicate(n, sample(1:(M - 1), k - 1))
  X <- apply(X, 2, sort)
  Y <- (X - X[c(k + 1, 1:k), ])[2:(k + 1), ]
  return(Y / M)
}

#' Sample n points, Dirichlet-distributed on a k-dimensional standard simplex
#' with means from 0 to 1, proportional to the numbers from the `alpha` argument.
#' Number of cell types is derived from `alpha`'s length. The `spread` argument
#' defines the deviation of generated values from one and zero.
#' @param n number of points
#' @param k dimensionality of generatred data
#' @param alpha proportional coefficient
#' @param spread deviation of generated values from one and zero.
sample_from_dirichlet <- function(n, k, alpha = NULL, spread = 7) {
  if (is.null(alpha)) {
    alpha <- rep(1, k)
  } else {
    stopifnot(length(alpha) == k)
  }
  l <- length(alpha)
  alpha_ <- alpha / sum(alpha) * spread
  x <- matrix(stats::rgamma(l * n, alpha_), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  t(x / as.vector(sm))
}

#' Generate proportions matrix H
#' @param n_samples number of samples (N)
#' @param n_cell_types number of cell types (K)
#' @param mean_proportions mean values for proportions
generate_proportions <- function(n_samples, n_cell_types, mean_proportions = NULL) {
  if (is.null(mean_proportions) || all(mean_proportions == mean_proportions[1])) {
    proportions <- sample_from_simplex_uniformly(n_samples, n_cell_types)
  } else {
    stopifnot(length(mean_proportions) == n_cell_types)
    proportions <- sample_from_dirichlet(n_samples, n_cell_types, mean_proportions)
  }
  colnames(proportions) <- paste0("sample_", 1:n_samples)
  rownames(proportions) <- paste0("cell_type_", 1:n_cell_types)
  
  proportions
}

apply_noise_to_data <- function(data, noise_deviation, protect_genes = c(), protect_samples = c(), noise_type= "norm_exp") {
  noise_mask <- matrix(stats::rnorm(length(data), mean = 0, sd = noise_deviation),
                       nrow = nrow(data), ncol = ncol(data))
  noise_mask[protect_genes, ] <- 0
  noise_mask[, protect_samples] <- 0
  rownames(noise_mask) <- rownames(data)
  colnames(noise_mask) <- colnames(data)

  if (noise_type == "norm_exp") {
    noise_mask <- 2 ^ noise_mask
    noisy_data <- data + noise_mask
  } else if (noise_type == "proportional") {
    noise_mask <- noise_mask + 1
    noisy_data <- data * noise_mask
  } else {
    noisy_data <- data + noise_mask
  }
  noisy_data[noisy_data < 0] <- 0
  noisy_data
}

generate_pure_samples <- function(basis, samples_per_cell_type, noise_deviation = 0.05) {
  data <- basis[, rep(seq_len(ncol(basis)), each = samples_per_cell_type)] * matrix(
    stats::rnorm(nrow(basis) * ncol(basis) * samples_per_cell_type, mean = 1, sd = noise_deviation),
    nrow = nrow(basis),
    ncol = ncol(basis) * samples_per_cell_type
  )
  proportions <- diag(1, ncol(basis), ncol(basis))[rep(seq_len(ncol(basis)), each = samples_per_cell_type),]
  colnames(proportions) <- colnames(basis)
  rownames(proportions) <- paste0("pure_sample_", seq_len(nrow(proportions)))
  proportions <- t(proportions)
  colnames(data) <- colnames(proportions)
  list(coefs = proportions, data = data)
}

generate_basis_samples <- function(basis) {
  proportions <- diag(1, ncol(basis), ncol(basis))
  rownames(proportions) <- colnames(basis)
  colnames(proportions) <- paste0("basis_sample_", seq_len(ncol(proportions)))
  data <- basis
  colnames(data) <- colnames(proportions)
  list(coefs = proportions, data = data)
}


