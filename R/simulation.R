#' Simulation builder initialization: simulation of basis, proportions and gene expression data without noise.
#' The motivation behind this builder approach is that every add... method creates a copy of the input simulation,
#' so that it is possible to create a single base simulation with many different child simulations.
#' For example, with different noise applied to the same data.
#' The add... functions inherently support magrittr's pipe operator (%>%).
#' @examples
#' sim_initial <- create_simulation(12000, 100, 3)
#' sim_noise_3 <- add_noise(sim_initial, 3)
#' sim_noise_5 <- sim_initial %>% add_noise(5) %>% add_basis_samples()
create_simulation <- function(
  n_genes,
  n_samples,
  n_cell_types,
  with_marker_genes = 0,
  mean_proportions = NULL,
  lbl_dataset = NULL
) {
  if (with_marker_genes != 0) {
    n_genes <- n_genes - n_cell_types * as.integer(with_marker_genes)
  }
  if (is.null(lbl_dataset)) {
    basis <- generate_basis(n_genes, n_cell_types)
  } else {
    basis <- generate_basis_lbl(n_genes, n_cell_types, lbl_dataset)
  }
  if (with_marker_genes != 0) {
    tmp <- append_marker_genes_to_basis(basis, as.integer(with_marker_genes))
    basis <- tmp[["basis"]]
    marker_gene_names <- tmp[["marker_gene_names"]]
  }
  
  proportions <- generate_proportions(n_samples, n_cell_types, mean_proportions)
  data <- basis %*% proportions
  
  data[data < 0] <- 0
  
  res <- list(
    basis = basis,
    proportions = proportions,
    data = data,
    mixed_sample_names = colnames(data)
  )
  if (with_marker_genes) {
    res$marker_gene_names <- marker_gene_names
  }
  
  res
}

#' Add N pure samples (basis + noise) per cell type to the simulation
add_pure_samples <- function(simulation, samples_per_cell_type, noise_deviation = 0.05) {
  pure_samples <- generate_pure_samples(
    simulation$basis,
    samples_per_cell_type,
    noise_deviation = noise_deviation
  )
  simulation$pure_sample_names <- colnames(pure_samples$data)
  simulation$data <- cbind(simulation$data, pure_samples$data)
  simulation$proportions <- cbind(simulation$proportions, pure_samples$proportions)
  simulation$pure_samples_noise_deviation <- noise_deviation
  simulation
}

#' Add samples, which are pure basis (without noise) to the simulation
add_basis_samples <- function(simulation) {
  basis_samples <- generate_basis_samples(simulation$basis)
  simulation$basis_sample_names <- colnames(basis_samples$data)
  simulation$data <- cbind(simulation$data, basis_samples$data)
  simulation$proportions <- cbind(simulation$proportions, basis_samples$proportions)
  simulation
}

#' Add noise to the simulation's data
add_noise <- function(simulation, noise_deviation, protect_genes = c(), protect_samples = c()) {
  simulation$data <- apply_noise_to_data(simulation$data, noise_deviation, protect_genes, protect_samples)
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  simulation$data_noise_deviation <- noise_deviation
  simulation
}


# All underlying logic is below
generate_basis <- function(n_genes, n_cell_types, sd_ = 0.2) {
  data_1 <- c(
    rnorm(n_genes * 2, mean = 4, sd = 0.75),
    rnorm(n_genes * 3, mean = 10, sd = 1.5)
  )
  basis <- matrix(0, nrow = n_genes, ncol = n_cell_types)
  sds_ <- rnorm(n_cell_types, mean = sd_, sd = 0.02)
  
  for (i in 1:n_genes) {
    basis[i, 1] <- data_1[sample(seq_len(length(data_1)), 1)] * rnorm(1, mean = 1, sd = sds_[1])
    for (j in 2:n_cell_types) {
      basis[i, j] <- basis[i, 1] * rnorm(1, mean = 1, sd = sds_[j])
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
  if (n_cell_types != 3) stop("LBL basis works only for 3 cell types")
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
sample_from_simplex_uniformly <- function(n, k, M = 100000) {
  X <- matrix(0, nrow = k + 1, ncol = n)
  X[k + 1,] <- M
  
  X[2:k,] <- replicate(n, sample(1:(M - 1), k - 1))
  X <- apply(X, 2, sort)
  Y <- (X - X[c(k + 1, 1:k),])[2:(k + 1),]
  return(Y / M)
}

#' Sample n points, Dirichlet-distributed on a k-dimensional standard simplex
#' with means from 0 to 1, proportional to the numbers from the `alpha` argument.
#' Number of cell types is derived from `alpha`'s length. The `spread` argument
#' defines the deviation of generated values from one and zero.
sample_from_dirichlet <- function(n, k, alpha = NULL, spread = 7) {
  if (is.null(alpha)) {
    alpha <- rep(1, k)
  } else {
    stopifnot(length(alpha) == k)
  }
  l <- length(alpha)
  alpha_ <- alpha / sum(alpha) * spread
  x <- matrix(rgamma(l * n, alpha_), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  t(x / as.vector(sm))
}

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

apply_noise_to_data <- function(data, noise_deviation, protect_genes = c(), protect_samples = c()) {
  noise_mask <- matrix(rnorm(length(data), sd = noise_deviation), nrow = nrow(data), ncol = ncol(data))
  rownames(noise_mask) <- rownames(data)
  colnames(noise_mask) <- colnames(data)
  noise_mask <- 2^noise_mask
  noise_mask[protect_genes,] <- 0
  noise_mask[, protect_samples] <- 0
  noisy_data <- data + noise_mask
  noisy_data[noisy_data < 0] <- 0
  noisy_data
}

generate_pure_samples <- function(basis, samples_per_cell_type, noise_deviation = 0.05) {
  data <- basis[, rep(seq_len(ncol(basis)), each = samples_per_cell_type)] * matrix(
    rnorm(nrow(basis) * ncol(basis) * samples_per_cell_type, mean = 1, sd = noise_deviation),
    nrow = nrow(basis),
    ncol = ncol(basis) * samples_per_cell_type
  )
  proportions <- diag(1, ncol(basis), ncol(basis))[rep(seq_len(ncol(basis)), each = samples_per_cell_type),]
  colnames(proportions) <- colnames(basis)
  rownames(proportions) <- paste0("pure_sample_", seq_len(nrow(proportions)))
  proportions <- t(proportions)
  colnames(data) <- colnames(proportions)
  list(proportions = proportions, data = data)
}

generate_basis_samples <- function(basis) {
  proportions <- diag(1, ncol(basis), ncol(basis))
  rownames(proportions) <- colnames(basis)
  colnames(proportions) <- paste0("basis_sample_", seq_len(ncol(proportions)))
  data <- basis
  colnames(data) <- colnames(proportions)
  list(proportions = proportions, data = data)
}
