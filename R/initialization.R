#' Initialize solution with selected initialization method
#'
#' In general it can use provided markers, solution for one of the simplex or random values
#'
#' @param proj dso$st$proj object containing containing all results of projection operation. (e.g. projected points and vectors)
#' @param strategy strategy to use for initialization. valid values are "select_x", "select_omega", "random", "random_invertible" and "marker_means"
#' @param kwargs put here marker gene names for each of the cell type if use  marker_means
#' @export
initialize_solution <- function(proj, strategy = "select_x", kwargs = NULL) {
  solution <- initializers[[strategy]](proj, kwargs)
  solution$init_strategy <- strategy
  return(solution)
}

# TODO: maybe not a good idea and it's better to just separate functions
initializers <- list(
  marker_means = function(proj, kwargs) {
    if (!"marker_list" %in% names(kwargs)) {
      spdl::error("Specify marker_list in kwargs for marker_means init.")
      stop("Specify marker_list in kwargs for marker_means init.")
    }
    mm <- lapply(kwargs$marker_list, function(ct_markers) {
      sel <- ct_markers[ct_markers %in% rownames(proj$X)]
      if (length(sel) == 0) {
        spdl::error("One or more cell types has zero markers present")
        stop("One or more cell types has zero markers present")
      }
    })
    K <- proj$meta$K
    X <- matrix(unlist(mm), ncol = K, byrow = T)

    rownames(X) <-  paste0("marker_derived", seq_len(nrow(X)))
    # Check if we have enough markers

    if (length(kwargs$marker_list) < K) {
        spdl::warn("Not enough markers were provided for the current K. Will set other vertices randomly with the highest center property after n attempts")
        if (!is.null(kwargs) && "n" %in% kwargs) {
            n <- kwargs$n
        } else {
         n <- 100
         }
        M <- proj$meta$M
        n_to_sample <- K - length(kwargs$marker_list)
        n <- 100
        idx_table_X <- matrix(0, ncol = n_to_sample + 1, nrow = n)
        for (i in 1:n) {
          #X
          ids_X <- sample(1:M, n_to_sample)
          extra_rows <- proj$X[ids_X, ]
          candidate_X <- rbind(X, extra_rows)
          metric_X <- sqrt(sum(apply(candidate_X[, -1], 2, mean) ^ 2))
          idx_table_X[i, ] <- c(ids_X, metric_X)
        }

        minrow <- function(mat) mat[which.min(mat[, ncol(mat)]), ]
        # X
        ids_X <- minrow(idx_table_X)[1:n_to_sample]
        X <- rbind(X, proj$X[ids_X, ])
    }
    return(set_solution_from_x(X, proj))
  },

  select_x = function(proj, kwargs = list()) {
    if ("features_subset" %in% kwargs) {
      points <- proj$X[kwargs$features_subset, ]
    } else {
      points <- proj$X
    }
    restored <- points %*% proj$meta$R
    n_cell_types <- proj$meta$K

    x <- t(points)
    u <- rowMeans(x)
    y <- x / matrix(kronecker(colSums(x * u), rep(1, n_cell_types)), nrow = n_cell_types)

    indice <- rep(0, n_cell_types)
    A <- matrix(0, nrow = n_cell_types, ncol = n_cell_types)
    A[n_cell_types, 1] <- 1
    for (i in 1:n_cell_types) {
      w <- matrix(runif(n_cell_types), ncol = 1)
      f <- w - A %*% MASS::ginv(A) %*% w
      f <- f / sqrt(sum(f^2))

      v <- t(f) %*% y
      indice[i] <- which.max(abs(v))
      A[, i] <- y[, indice[i]]
    }

    Ae <- restored[indice,]
    X <- Ae %*% t(proj$meta$R)
    return(set_solution_from_x(X, proj))
  },

  select_omega = function(proj, kwargs = list()) {
    if ("samples_subset" %in% kwargs) {
      points <- proj$Omega[kwargs$samples_subset, ]
    } else {
      points <- proj$Omega
    }

    restored <- t(proj$meta$S) %*% t(points)
    n_cell_types <- proj$meta$K

    x <- t(points)
    u <- rowMeans(x)
    y <- x / matrix(kronecker(colSums(x * u), rep(1, n_cell_types)), nrow=n_cell_types)

    indice <- rep(0, n_cell_types)
    A <- matrix(0, nrow=n_cell_types, ncol=n_cell_types)
    A[n_cell_types, 1] <- 1
    for (i in 1:n_cell_types) {
      w <- matrix(runif(n_cell_types), ncol=1)
      f <- w - A %*% MASS::ginv(A) %*% w
      f <- f / sqrt(sum(f^2))

      v <- t(f) %*% y
      indice[i] <- which.max(abs(v))
      A[, i] <- y[, indice[i]]
    }
    Ae <- restored[, indice]
    Omega <- proj$meta$S %*% Ae

    M <- proj$meta$M
    N <- proj$meta$N
    B <- matrix(apply(proj$meta$S, 1, sum), ncol=1, nrow=n_cell_types)
    D_w <- MASS::ginv(Omega) %*% B
    D_h <- D_w * (N/M)
    V__ <- proj$meta$S %*% proj$X
    X <- MASS::ginv(Omega %*% diag(D_w[, 1])) %*% V__
    return(list(
      X = X,
      Omega = t(Omega),
      D_w = D_w,
      D_h = D_h
    ))
  },

  random = function(proj, kwargs = NULL) {
    if (!is.null(kwargs) && "n" %in% kwargs) {
      n <- kwargs[["n"]]
    } else {
      n <- 100
    }
    n_cell_types <- proj$meta$K
    M <- proj$meta$M
    N <- proj$meta$N

    idx_table_X <- matrix(0, ncol = n_cell_types + 1, nrow = n)
    idx_table_Omega <- matrix(0, ncol = n_cell_types + 1, nrow = n)

    for (i in 1:n) {
      #Omega
      ids_Omega <- sample(1:N, n_cell_types)
      Omega <- t(proj$Omega)[, ids_Omega]
      metric_Omega <- sqrt(sum(apply(Omega[-1, ], 1, mean) ^ 2))
      idx_table_Omega[i, ] <- c(ids_Omega, metric_Omega)

      #X
      ids_X <- sample(1:M, n_cell_types)
      X <- proj$X[ids_X, ]
      metric_X <- sqrt(sum(apply(X[, -1], 2, mean) ^ 2))
      idx_table_X[i, ] <- c(ids_X, metric_X)
    }

    minrow <- function(mat) mat[which.min(mat[, ncol(mat)]), ]

    # Omega
    ids_Omega <- minrow(idx_table_Omega)[1:n_cell_types]
    Omega <- t(proj$Omega)[, ids_Omega]

    # X
    ids_X <- minrow(idx_table_X)[1:n_cell_types]
    X <- proj$X[ids_X, ]

    Ds <- get_Dwh_from_XOmega(X, Omega, proj)
    colnames(X) <- colnames(proj$X)
    rownames(X) <- paste("Init ", c(1: n_cell_types))
    rownames(Omega) <- colnames(proj$Omega)
    colnames(Omega) <- paste("Init ", c(1: n_cell_types))
    return(list(
      X = X,
      Omega = t(Omega),
      D_w = Ds$D_w,
      D_h = Ds$D_h
    ))
  },
  random_invertible = function(proj, kwargs = NULL) {
    n_cell_types <- proj$meta$K
    M <- proj$meta$M
    N <- proj$meta$N
    sigma_1 <-  sqrt(M / N)
    # We want to ensure that X_dtilda and Omega_dtilda are inverse
    # and have all positive first row/column respectively

    # First we generate d candidates  generationg vector whos elements sum to M
    d_elements <-  runif(n_cell_types, min = 1e-5, max = 1.0)
    d_elements <- M * d_elements / sum(d_elements)

    # r1 (row of omega) will be  sqrt(d) / sqrt(sigma_1) / sqrt(M)
    r1 <- sqrt(d_elements) / (sqrt(sigma_1) * sqrt(M))
    # c1 is predefined. its sqrt(d) / sqrt(sigma_1) / sqrt(N)
    c1 <- sqrt(d_elements) / (sqrt(sigma_1) * sqrt(N))
    if (!((sum(r1 * c1) - 1) < 1e-6)) {
      spdl::error("Error in first row/column initialization.  r1 and c1 should have specific relation")
    }
    # construct the null space vectors of the r1.
    # N_r orthonormal and orthogonal to r1
    N_r <- MASS::Null(matrix(r1, ncol = 1))
    # However orthonormality of N_r is not necessary for us
    # We break it with the random invertible matrix
    dim_null <- ncol(N_r)
    W <- matrix(runif(dim_null * dim_null, min = -1.0, max = 1.0), nrow = dim_null)
    # Force it to be diagonally dominant so it is computationally stable
    diag(W) <- rowSums(W) + 1.0
    # Scramble the null space. V is now linearly independent, in the null space,
    # but NO LONGER orthonormal.
    V <- N_r %*% W
    X_dtilda <- cbind(c1, V)
    # and omega should be the inverse of this matrix
    Omega_dtilda <- solve(X_dtilda)

    # Now let's combine all of this to the new matrix X and Omega
    # X should be   * 1/sqrt(D)  * X_d_tilda   * sqrt(Sigma)
    X <- diag(1 / sqrt(d_elements)) %*%  X_dtilda %*% sqrt(proj$meta$Sigma)
    # Omega should be sqrt(Sigma) * Omega_dtilda * 1/sqrt(D)
    Omega <-  sqrt(proj$meta$Sigma) %*% Omega_dtilda %*% diag(1/sqrt(d_elements))
    # --- ASSERTS ---
    identity_check <- Omega_dtilda %*% X_dtilda
    target_identity <- diag(nrow(X_dtilda))
    # Check inverse
    if (max(abs(identity_check - target_identity)) > 1e-6) {
      spdl::error("X_dtilda and Omega_dtilda are not mathematical inverses of each other.")
    }
    # Check r1 c1 property
    new_r1 <- Omega_dtilda[1, ]
    if (!(abs(sum(new_r1 * c1) - 1) < 1e-6)) {
      spdl::error("We lost r1 vs c1 properties in the middle")
    }
    colnames(X) <- colnames(proj$X)
    rownames(X) <- paste("Init ", c(1: n_cell_types))
    rownames(Omega) <- colnames(proj$Omega)
    colnames(Omega) <- paste("Init ", c(1: n_cell_types))
    # ------
    Dw <- as.matrix(d_elements)
    Dh <- Dw * (N / M)
    return(list(
      X = X,
      Omega = t(Omega),
      D_w = Dw,
      D_h = Dh
    ))

  },
  random_symmetric = function(proj, kwargs = NULL) {
    if (!is.null(kwargs) && "n" %in% kwargs) {
      n <- kwargs[["n"]]
    } else {
      n <- 100
    }
    n_cell_types <- proj$meta$K
    M <- proj$meta$M
    N <- proj$meta$N

    idx_table_X <- matrix(0, ncol = n_cell_types + 1, nrow = n)

    for (i in 1:n) {
      #X
      ids_X <- sample(1:M, n_cell_types)
      X <- proj$X[ids_X, ]
      metric_X <- sqrt(sum(apply(X[, -1], 2, mean) ^ 2))
      idx_table_X[i, ] <- c(ids_X, metric_X)
    }

    minrow <- function(mat) mat[which.min(mat[, ncol(mat)]), ]


    # X
    ids_X <- minrow(idx_table_X)[1:n_cell_types]
    X <- proj$X[ids_X, ]
    Omega <- t(proj$X[ids_X, ])

    Dw <- MASS::ginv(t(X)) %*% proj$meta$A
    Dh <- Dw *(N / M);

    Ds <- get_Dwh_from_XOmega(X, Omega, proj)
    return(list(
      X = X,
      Omega = t(Omega),
      D_w = Dw,
      D_h = Dh
    ))
  },
  random_X_within_theta_angle = function(proj, kwargs) {
    # Having center points provided we want to initialize in some
    # random point within theta angle
    if (!"init_centers" %in% names(kwargs)) {
      spdl::error("Specify init_centers with X and optionaly Omega in kwargs to set
      some starting point for this initialization.")
      stop("Specify init_centers with X and optionaly Omega in kwargs to set
      some starting point for this initialization.")
    }
    if (!"theta" %in% names(kwargs)) {
      spdl::error("Specify theta in kwargs to set some constraint on theta.")
      stop("Specify theta in kwargs to set some constraint on theta.")
    }
    max_length <- 1.5
    n_cell_types <- proj$meta$K
    init_centers <- kwargs$init_centers
    center_points <- init_centers$X
    theta_max_div <- kwargs$theta
    generated_vectors <- lapply(c(1:n_cell_types),  function(cell_type) {
      # original vector is a medoid for cell type (transposed)
      orig_vector <-  as.matrix(center_points[cell_type, ])
      # generate rotated vector
      cur_theta <- runif(1, min = -theta_max_div, max = theta_max_div)
      random_multiplier <-  runif(1, min = 0.4, max = max_length)
      u <-  as.matrix(orig_vector[2:n_cell_types])
      u <-  random_multiplier * u
      # random direction
      p <- as.matrix(rnorm(n_cell_types - 1, mean = 0, sd = 1))
      nom <- t(p) %*% u
      denom <-  t(u) %*% u
      multiplier <- (nom / denom)[1]
      w_prime <- p - multiplier * u
      w <- norm(u) * w_prime / norm(w_prime)
      v <- cos(cur_theta) * u + sin(cur_theta) * w   # new vector
      result_vector <- c()
      result_vector[2:n_cell_types] <- v
      result_vector[1] <- orig_vector[1]
      return(result_vector)
    })
    X <- t(do.call("cbind", generated_vectors))
    colnames(X) <- colnames(proj$X)
    # Set Omega, same or not same
    if ("Omega" %in% names(init_centers)) {
      center_points <- init_centers$Omega
      generated_vectors <- lapply(c(1:n_cell_types),  function(cell_type) {
        # original vector is a medoid for cell type (transposed)
        orig_vector <-  as.matrix(center_points[cell_type, ])
        # generate rotated vector
        cur_theta <- runif(1, min = -theta_max_div, max = theta_max_div)
        random_multiplier <-  runif(1, min = 0.4, max = max_length)
        u <-  as.matrix(orig_vector[2:n_cell_types])
        u <-  random_multiplier * u
        # random direction
        p <- as.matrix(rnorm(n_cell_types - 1, mean = 0, sd = 1))
        nom <- t(p) %*% u
        denom <-  t(u) %*% u
        multiplier <- (nom / denom)[1]
        w_prime <- p - multiplier * u
        w <- norm(u) * w_prime / norm(w_prime)
        v <- cos(cur_theta) * u + sin(cur_theta) * w   # new vector
        result_vector <- c()
        result_vector[2:n_cell_types] <- v
        result_vector[1] <- orig_vector[1]
        return(result_vector)
      })
      Omega <- t(do.call("cbind", generated_vectors))
      colnames(Omega) <- colnames(proj$Omega)
      Ds <- get_Dwh_from_XOmega(X, Omega, proj)
      new_solution_proj <- list(
        X = X,
        Omega = Omega,
        D_w = Ds$D_w,
        D_h = Ds$D_h
      )
    } else {
      new_solution_proj <- set_solution_from_x(X, proj)
    }
    return(new_solution_proj)
  }
)

#' Infer solution matrices Omega and D from predefined matrix X
#'
#' Uses properties of sinkhorn transformed matrices (their proportional relation)
#'
#' @param X KxK matrix X
#' @param proj dso$st$proj object, containing all necessary info about projection (e.g. vectors and projected points)
#' @export
set_solution_from_x <- function(X, proj) {
  D_h <- MASS::ginv(t(X)) %*% proj$meta$A
  M <- proj$meta$M
  N <- proj$meta$N
  D_w <- D_h * (M / N)
  V__ <- proj$meta$S %*% proj$X
  Omega <- V__ %*% MASS::ginv(diag(D_w[, 1]) %*% X)
  return(list(
    D_h = D_h,
    D_w = D_w,
    X = X,
    Omega = t(Omega)
  ))
}

#' Infer solution matrices D from predefined matrices X and Omega
#'
#' Uses properties of sinkhorn transformed matrices. Performs NNLS
#'
#' @param X KxK solution matrix X
#' @param Omega KxK solution matrix Omega
#' @param proj dso$st$proj object, containing all necessary info about projection (e.g. vectors and projected points)
#' @export
get_Dwh_from_XOmega <- function(X, Omega, proj) {
  V__ <- proj$meta$S %*% proj$X

  K <- proj$meta$K
  M <- proj$meta$M
  N <- proj$meta$N

  ## calculate D_w and D_h
  ## vectorizing deconvolution
  vec_mtx <- matrix(0, K ^ 2, K)
  for (col in 1:K) {
    vec_mtx[, col] <- cbind(c(t(t(Omega[, col])) %*% X[col, ]))
  }

  ## adding sum-to-one constraint and solving NNLS
  D_w <-
    matrix(nnls_C__(rbind(vec_mtx, Omega), rbind(cbind(c(
      V__
    )), proj$meta$B)),
    nrow = K,
    ncol = 1)
  D_h <- D_w * (N / M)
  return(list(D_w = D_w, D_h = D_h))
}
