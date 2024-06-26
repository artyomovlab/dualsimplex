% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{add_pure_samples}
\alias{add_pure_samples}
\title{Add N pure samples (basis + noise) per cell type to the simulation}
\usage{
add_pure_samples(simulation, samples_per_cell_type, noise_deviation = 0.05)
}
\arguments{
\item{simulation}{original simulation object}

\item{samples_per_cell_type}{number of pure samples to add for each hidden main component}

\item{noise_deviation}{noise deviation value for pure samples}
}
\value{
modified sim object
}
\description{
Add N pure samples (basis + noise) per cell type to the simulation
}
