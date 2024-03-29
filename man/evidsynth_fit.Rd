% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evidsynth_fit.R
\name{evidsynth_fit}
\alias{evidsynth_fit}
\title{Evidence synthesis model fit}
\usage{
evidsynth_fit(
  prevalence_data,
  progression_data,
  model = c("cut", "independent", "joint"),
  n_iter = 2000,
  n_burnin = 100,
  n_thin = 10
)
}
\arguments{
\item{prevalence_data}{Latent infection prevalence data}

\item{progression_data}{Disease progression data}

\item{model}{Which model to use, from "cut", "independent", or "joint"}

\item{n_iter}{Number of iterations}

\item{n_burnin}{Number of iterations in burn-in}

\item{n_thin}{Number of iterations to thin}
}
\value{
List
}
\description{
Evidence synthesis model fit
}
