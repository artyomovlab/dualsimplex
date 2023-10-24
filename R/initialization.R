initialize_solution <- function(proj, strategy = "select_x", kwargs = NULL) {
  solution <- initializers[[strategy]](proj, kwargs)
  solution$init_strategy <- strategy
  return(solution)
}

# TODO: maybe not a good idea and it's better to just separate functions
initializers <- list(
  marker_means = function(proj, kwargs) {
    if (!"marker_list" %in% names(kwargs)) {
      stop("Put marker_list in kwargs for marker_means init.")
    }
    mm <- lapply(kwargs$marker_list, function(ct_markers) {
      sel <- ct_markers[ct_markers %in% rownames(proj$X)]
      if (length(sel) == 0) stop("One or more cell types has zero markers present")
      colMeans(proj$X[sel, ])
    })
    X <- matrix(unlist(mm), ncol = length(kwargs$marker_list), byrow = T)
    return(set_solution_from_x(X, proj))
  },

  select_x = function(proj, kwargs = list()) {
    if ("genes_subset" %in% kwargs) {
      points <- proj$X[kwargs$genes_subset, ]
    } else {
      points <- proj$X
    }
    restored <- points %*% proj$meta$R
    n_cell_types <- proj$meta$K

    x <- t(points)
    u <- rowMeans(x)
    y <-
      x / matrix(kronecker(colSums(x * u), rep(1, n_cell_types)), nrow = n_cell_types)

    indice <- rep(0, n_cell_types)
    A <- matrix(0, nrow = n_cell_types, ncol = n_cell_types)
    A[n_cell_types, 1] <- 1
    for (i in 1:n_cell_types) {
      w <- matrix(runif(n_cell_types), ncol = 1)
      f <- w - A %*% MASS::ginv(A) %*% w

      f <- f / sqrt(sum(f ^ 2))

      v <- t(f) %*% y
      indice[i] <- which.max(abs(v))
      A[, i] <- y[, indice[i]]
    }

    Ae <- restored[indice,]
    X <- Ae %*% t(proj$meta$R)
    return(set_solution_from_x(X, proj))
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
    return(list(
      X = X,
      Omega = t(Omega),
      D_w = Ds$D_w,
      D_h = Ds$D_h
    ))
  }
)

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

  ## adding sum-to-one constraint
  D_w <-
    matrix(nnls::nnls(rbind(vec_mtx, Omega), rbind(cbind(c(
      V__
    )), proj$meta$B))$x,
    nrow = K,
    ncol = 1)
  D_h <- D_w * (N / M)
  return(list(D_w = D_w, D_h = D_h))
}
