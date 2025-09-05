# Required functions:
# 2. Plot specified dims of projected data, optionally colored by a passed column or highlighted genes
# 3. Plot solution for both spaces
# 4. Plot solution history
# 5. Plot solution history animated (with requireNamespace())
# Functions should accept both umap and dims
# Final function should work with a data.frame with x_col, y_col, color_col and group_col specified
# Function, accepting proj, solution_proj, proj_umap should be present
# Solution should be in form of colored points, added to initial plots
# Should be a function to perform one of two projections, and both of them

############ EXTRACTING 2D DATASET ############

#' Get required two-dimensional subset of the data
#'
#' Conditionally this method returns 2D subset of points in original SVD space or umap transformed.
#' If extra_points_proj specified it will use these points for visualization and proj is just information about projection
#'
#' @param proj dso$st$proj object containing all projected points and projection info
#' @param use_dims dimensions to extract (e.g. 2:3). NULL to get umap points
#' @param extra_points_proj additional points (should contain two matrices $X and $Omega)
#' @return subset of projected points (columns)
#' @export
get_2d_subset <- function(proj, use_dims, extra_points_proj = NULL) {
  if (is.null(use_dims)) {
    if ("umap" %in% names(proj)) {
      if (!is.null(extra_points_proj)) {
        return(transform_proj_umap(extra_points_proj, proj))
      } else {
        return(proj$umap[c("X", "Omega")])
      }
    } else {
      if (proj$meta$K == 3) {
        stop("Specify use_dims (e.g. 2:3) or calculate umap")
      } else {
        use_dims <- 2:3
      }
    }
  }
  if (!is.null(extra_points_proj)) {
    proj <- extra_points_proj
  }
  return(list(
    X = proj$X[, use_dims],
    Omega = proj$Omega[, use_dims]
  ))
}

#' Util function to merge list of dataframes to a single one with column to separate original dataframes
#'
#' @param data_list the list of dataframes
#' @param group_colname new column to create in both. column values will be taken from names(data_list)
#' @return merged data
#' @export
concat_data <- function(data_list, group_colname) {
  merged <- lapply(1:length(data_list), function(i) {
    this_data <- as.data.frame(data_list[[i]])
    this_data[, group_colname] <- names(data_list)[[i]]
    return(this_data)
  })
  merged <- as.data.frame(do.call(rbind, merged))
  merged[, group_colname] <- factor(
    merged[, group_colname],
    levels = sort(unique(merged[, group_colname]), decreasing = T)
  )
  return(merged)
}

#' Get history of values for X and Omega from optimization logs
#'
#' Points will be taken uniformly in the
#'
#' @param solution_proj dso$st$solution_proj object containing all optimization logs and result
#' @param step step to skip some values
#' @param from_iter starting point
#' @param to_iter end point
#' @return X and Omega values for desired timestamps
#' @export
get_solution_history <- function(solution_proj, step, from_iter = 1, to_iter = NULL) {
  stats <- list()
  stats$X <- solution_proj$optim_history$points_statistics_X
  stats$Omega <- solution_proj$optim_history$points_statistics_Omega
  nct <- sqrt(ncol(stats$Omega))
  nit <- nrow(stats$Omega)
  if (is.null(to_iter)) {
    to_iter <- nit
  }
  # TODO: correct order initially
  correct_order <- order((1:ncol(stats$Omega) - 1) %% nct)
  stats$Omega <- stats$Omega[, correct_order]

  solution_history <- lapply(stats, function(mat) {
    tran <- matrix(c(t(mat)), ncol = nct, byrow = T)
    tran <- cbind(tran, rep(1:nct, times = nit))
    tran <- cbind(tran, rep(1:nit, each = nct))
    colnames(tran) <- c(paste0("dim_", 1:nct), "point", "iter")
    tran <- tran[(tran[, "iter"] >= from_iter) & (tran[, "iter"] <= to_iter), ]
    tran <- tran[((tran[, "iter"] %% step) == 0) | (tran[, "iter"] == to_iter) | (tran[, "iter"] == from_iter), ]
    tran <- as.data.frame(tran)
    tran[, "point"] <- as.factor(tran[, "point"])
    return(tran)
  })
  return(solution_history)
}

############ EXTRACTING COLOR VECTOR & SCHEME ############

#' Prepare a list of color values for points we want to plot
#'
#' @param to_plot dataframe to plot
#' @param color color column or values
#' @param color_name color column or values
#' @return list of values to be used for coloring
get_color_params <- function(to_plot, color, color_name) {
  base_color <-  "black"
  if (is.null(color)) {
    return(list(
      color_vec = NULL,
      color_name = NULL,
      color_scheme = "default",
      base_color = base_color
    ))
  }

  if (is.null(color_name)) color_name <- "annotation"
  if (length(color[[1]]) > 1) color <- unlist(color)

  if (length(color) == nrow(to_plot)) {
    if (is.factor(color) || is.logical(color)) {
      return(list(
        color_vec = color,
        color_name = color_name,
        color_scheme = "factor",
        base_color = NULL
      ))
    } else {
      return(list(
        color_vec = color,
        color_name = color_name,
        color_scheme = "direct_multiple_colors",
        base_color = NULL
      ))
    }
  } else {
    if (color_name == "direct_single_color") {
      base_color <-  color
      return (list(
        color_vec = NULL,
        color_name = NULL,
        color_scheme = "direct_single_color",
        base_color = base_color
      ))
    }
    if (color_name == "individual_highlight") {
      color_name <- "highlight"
      return(list(
        color_vec = rownames(to_plot) %in% color,
        color_name = color_name,
        color_scheme = "highlight",
        base_color = base_color
    ))
    }
    if (color_name == "annotation")
    return(list(
      color_vec = rownames(to_plot) %in% color,
      color_name = "highlight",
      color_scheme = "highlight",
      base_color = base_color
    ))
  }
}


############ PLOTTING ############

#' The starting point of plotting projected points.
#'
#' This funtion can be used to visualize projected points. It handles turning projection into 2-dim data
#' and grouping the plots by X and Omega uses ggplot facet_wrap, can add solution easily
#' Any ggplot arguments can be added after main arguments
#'
#'
#' @param proj dso$st$proj object containing all projected points
#' @param use_dims dimensions to use (e.g. 2:3), NULL means UMAP if umap is calculated
#' @param spaces which space to visualize "X" or "Omega" or c("X", "Omega")
#' @param ... any additional geom_point_params
#' @return requested plot
#' @export
plot_projection_points <- function(
  proj,
  use_dims = NULL,
  spaces = c("X", "Omega"),
  ...
) {
  points_2d <- get_2d_subset(proj, use_dims)[spaces]
  if (length(spaces) > 1) {
    points_2d <- concat_data(points_2d, "space")
  } else {
    points_2d <- as.data.frame(points_2d[[1]])
  }
  plt <- plot_points_2d(points_2d, ...)
  if (length(spaces) > 1) {
    plt <- plt + facet_wrap("space", scales = "free")
  }
  return(plt)
}

#' Add solution points (X, Omega) to the plot
#'
#' @param plt plot to add
#' @param solution_proj dso$st$solution_proj containing optimization history and solution
#' @param proj dso$st$proj containing info about projection (e.g. vectors, projected points)
#' @param use_dims dimensions of the proj to use (e.g. 2:3 or NULL for umap)
#' @param pt_size size of solution points
#' @param spaces for which space to add c("X", "Omega")
#' @return modified plot(s)
#' @export
add_solution <- function(
  plt,
  solution_proj,
  proj,
  use_dims = NULL,
  pt_size = 3,
  spaces = c("X", "Omega")
) {
  points_2d <- get_2d_subset(proj, use_dims, solution_proj)[spaces]
  points_2d <- lapply(
    points_2d,
    function(pts) cbind(pts, point = 1:nrow(pts))
  )
  if (length(spaces) > 1) {
    points_2d <- concat_data(points_2d, "space")
  } else {
    points_2d <- as.data.frame(points_2d[[1]])
  }
  points_2d$point <- as.factor(points_2d$point)
  x_col <- colnames(points_2d)[[1]]
  y_col <- colnames(points_2d)[[2]]
  plt <- plt + geom_point(
    data = points_2d,
    aes_string(x_col, y_col, fill = "point"),
    color = "black",
    pch = 21,
    size = pt_size
  ) + theme(legend.position = "none")
  return(plt)
}

#' Add solution history lines/points to the plot
#'
#' @param plt plot to add
#' @param solution_proj dso$st$solution_proj containing optimization history and solution
#' @param proj dso$st$proj containing info about projection (e.g. vectors, projected points)
#' @param use_dims dimensions of the proj to use (e.g. 2:3 or NULL for umap)
#' @param pt_size size of history points
#' @param pt_opacity alpha for lines
#' @param step how many values to skip
#' @param from_iter starting point
#' @param to_iter end point
#' @param spaces for which space to add c("X", "Omega")
#' @param colored color according to solution points
#' @param path if TRUE do lines, else points
#' @return modified plot(s)
#' @export
add_solution_history <- function(
  plt,
  solution_proj,
  proj,
  use_dims = NULL,
  pt_size = 0.5,
  pt_opacity = 0.95,
  step = 100,
  from_iter = 1,
  to_iter = NULL,
  spaces = c("X", "Omega"),
  colored = TRUE,
  path = TRUE
) {
  solution_history <- get_solution_history(solution_proj, step, from_iter = from_iter, to_iter = to_iter)
  points_2d <- get_2d_subset(proj, use_dims, solution_history)[spaces]
  if (length(spaces) > 1) {
    solution_history <- concat_data(solution_history[spaces], "space")
    points_2d <- concat_data(points_2d[spaces], "space")
  } else {
    solution_history <- solution_history[[1]]
    points_2d <- as.data.frame(points_2d[[1]])
  }
  points_2d <- cbind(points_2d, solution_history[c("point", "iter")])
  x_col <- colnames(points_2d)[[1]]
  y_col <- colnames(points_2d)[[2]]

  geom <- if (path) geom_path else geom_point

  if (colored) {
    plt <- plt + geom(
      data = points_2d,
      aes_string(x_col, y_col, col = "point"),
      size = pt_size,
      alpha = pt_opacity
    ) + theme(legend.position = "none")
  } else {
    plt <- plt + geom(
      data = points_2d,
      aes_string(x_col, y_col, group = "point"),
      color = "black",
      size = pt_size,
      alpha = pt_opacity
    )
  }

  plt <- plt +
    labs(title = paste0(max(solution_history$iter), " iterations"))

  return(plt)
}

#' Intermediate step of plotting
#' This function deals with color and arguments cleaning
#'
#' @param points_2d 2d subset of the data
#' @param color color argument (row names, or column from annotation or maybe vector of values)
#' @param color_name color column name in the preprocessed data (can be column from annotation or "highlight" or "annotation"
#' @param order_by_color should we order by color_name column?
#' @param ... any additional geom_point params
#' @return ggplot object
plot_points_2d <- function(
  points_2d,
  color = NULL,
  color_name = NULL,
  order_by_color = TRUE,
  ...
) {
  to_plot <- as.data.frame(points_2d)
  x_col <- colnames(to_plot)[1]
  y_col <- colnames(to_plot)[2]

  cp <- get_color_params(points_2d, color, color_name)

  if ((cp$color_scheme != "default") & (cp$color_scheme != "direct_single_color")) {
    to_plot[, cp$color_name] <- cp$color_vec
    if (order_by_color)
      to_plot <- to_plot[order(to_plot[, cp$color_name]), ]
  }

  return(plot_points_2d_clean(
    to_plot,
    x_col,
    y_col,
    cp$color_name,
    cp$color_scheme,
    cp$base_color,
    ...
  ))
}

#' Final step of plotting
#' This function just does ggplot object
#' You can add any geom_point params after main arguments
#'
#' @param to_plot clean 2D dataframe to plot
#' @param x_col x axe column
#' @param y_col y axe column
#' @param color_col which column is color
#' @param color_scheme color scheme identified by color params parsing
#' @param base_color for points
#' @param pt_opacity opasity for points
#' @param pt_size point size
#' @param ... any additional geom_point params
#' @return ggplot object
plot_points_2d_clean <- function(
  to_plot,
  x_col,
  y_col,
  color_col,
  color_scheme,
  base_color,
  pt_opacity = 0.2,
  pt_size = 1,
  ...
) {
  plt <- ggplot(to_plot, aes_string(x = x_col, y = y_col, color = color_col))


  if (color_scheme == "default" || color_scheme == "highlight" || color_scheme == "direct_single_color" ) {
    plt <- plt + rasterize_if_needed(geom_point(size = pt_size, color = base_color, alpha = pt_opacity))
    if (color_scheme == "highlight") {
      plt <- plt +
      rasterize_if_needed(geom_point(
        data = to_plot[to_plot[[color_col]], ],
        size = pt_size * 2,
        alpha = 1,
        color = "green",
        ...
      ))
    }
  } else if (color_scheme == "direct_multiple_colors") {
    plt <- plt +
      rasterize_if_needed(geom_point(size = pt_size, alpha = min(1, pt_opacity * 4), ...)
      + scale_color_distiller(palette = "Spectral", name = color_col)
      )

  } else if (color_scheme == "factor") {
    plt <- plt +
    rasterize_if_needed(geom_point(size = pt_size, color = "grey70", alpha = pt_opacity, ...)) +
    rasterize_if_needed(geom_point(
        data = to_plot[!is.na(to_plot[[color_col]]), ],
        alpha = 1,
        size = pt_size * 2,
        ...
      ))
  } else {
    stop("Invalid color_scheme")
  }

  plt <- plt + theme_minimal()

  return(plt)
}

#' Produce animation gif of optimization with gganimate (you need to install it)
#'
#' @param proj dso$st$proj containing projected points and info about projection
#' @param solution_proj dso$st$solution_proj containing optimiation history
#' @param step how many steps we skip
#' @param height height of result animation
#' @param width width of result animation
#' @param gif_filename filname
#' @param gif_dir directory
#' @param nframes how many frames in total
#' @param use_dims which dimensions to show (e.g. 2:3 or 3:4)
#' @return ganimate::animate or nothing if file is specified
#' @export
plot_solution_history_anim <- function(
  proj,
  solution_proj,
  step = 100,
  height = 600,
  width = 1200,
  gif_filename = NULL,
  gif_dir = NULL,
  nframes = 300,
  use_dims = NULL
) {
  plt <- plot_projection_points(proj, use_dims = use_dims) %>%
    add_solution_history(
      solution_proj, proj, step = 100, pt_size = 4, use_dims = use_dims, path = FALSE
    ) +
    labs(title = 'Iteration: {as.integer(closest_state)}') +
    gganimate::transition_states(.data$iter, transition_length = 2, wrap = F) +
    gganimate::ease_aes("linear")

  if (is.null(gif_filename) || is.null(gif_dir)) {
    gganimate::animate(plt, height = height, width = width, nframes = nframes, renderer = gganimate::gifski_renderer())
  } else {
    gganimate::anim_save(gif_filename, plt, height = height, width = width, nframes = nframes, renderer = gganimate::gifski_renderer())
  }
}
