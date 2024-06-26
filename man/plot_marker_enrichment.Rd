% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_analysis.R
\name{plot_marker_enrichment}
\alias{plot_marker_enrichment}
\title{Plot enrichement for markers}
\usage{
plot_marker_enrichment(
  so,
  marker_names,
  ncol = 4,
  limits = NULL,
  ggadd = function(plt, i) plt,
  wrap = T,
  ...
)
}
\arguments{
\item{so}{Seurat object}

\item{marker_names}{marker genes for each cell type (list of lists)}

\item{ncol}{number of columns in result plot (cowplot::plot_grid)}

\item{limits}{color limits (to adjust color)}

\item{ggadd}{add text to each plot (function taking plot and index of plot)}

\item{wrap}{whould we do grid or not}

\item{...}{any valid Seurat::FeaturePlot params}
}
\value{
single plot or multiple plots
}
\description{
Using Seurat::FeaturePlot to plot
}
