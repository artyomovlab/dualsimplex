% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{add_noise}
\alias{add_noise}
\title{Add noise to the simulation's data}
\usage{
add_noise(
  simulation,
  noise_deviation,
  protect_genes = c(),
  protect_samples = c()
)
}
\arguments{
\item{simulation}{original simulation object}

\item{noise_deviation}{noise deviation}

\item{protect_genes}{do not modigy these rows}

\item{protect_samples}{do not modigy these samples}
}
\value{
modified sim object
}
\description{
Add noise to the simulation's data
}
