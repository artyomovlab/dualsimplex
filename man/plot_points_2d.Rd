% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projected_plots.R
\name{plot_points_2d}
\alias{plot_points_2d}
\title{Intermediate step of plotting
This function deals with color and arguments cleaning}
\usage{
plot_points_2d(
  points_2d,
  color = NULL,
  color_name = NULL,
  order_by_color = TRUE,
  ...
)
}
\arguments{
\item{points_2d}{2d subset of the data}

\item{color}{color argument (row names, or column from annotation or maybe vector of values)}

\item{color_name}{color column name in the preprocessed data (can be column from annotation or "highlight" or "annotation"}

\item{order_by_color}{should we order by color_name column?}

\item{...}{any additional geom_point params}
}
\value{
ggplot object
}
\description{
Intermediate step of plotting
This function deals with color and arguments cleaning
}
