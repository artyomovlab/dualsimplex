% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projection.R
\name{svd_project}
\alias{svd_project}
\title{Prepare projection object}
\usage{
svd_project(scaling, dims, ops = NULL)
}
\arguments{
\item{scaling}{dso$st$scaling object containing sinkhorn scaling result}

\item{dims}{how many dimensions we want to get}

\item{ops}{= SVD projection operations. Default is NULL and calculated on the fly for backward compatability.}
}
\value{
proj object
}
\description{
Entry point to produce dso$st$proj object, containing all information about projection as well as projected points
}
