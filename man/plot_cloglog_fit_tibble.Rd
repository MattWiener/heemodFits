% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_fit_output.R
\name{plot_cloglog_fit_tibble}
\alias{plot_cloglog_fit_tibble}
\title{Assemble data from a fit tibble and plot a log cumulative hazard plot}
\usage{
plot_cloglog_fit_tibble(fit_tibble, treatments = NULL, set_name, type, ...)
}
\arguments{
\item{fit_tibble}{the output of \code{\link[=survival_fits_from_tabular]{survival_fits_from_tabular()}}.}

\item{treatments}{the treatments to include; if left \code{NULL},
all treatments will be included.}

\item{set_name}{subset name}

\item{type}{PFS or OS}

\item{...}{additional parameters to pass to \code{\link[=plot_fit_data]{plot_fit_data()}}}
}
\value{
a \code{ggplot} plot
}
\description{
Assemble data from a fit tibble and plot a log cumulative hazard plot
}
