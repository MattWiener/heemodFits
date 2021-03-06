% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience.R
\name{find_least_cost_partition}
\alias{find_least_cost_partition}
\title{Find lowest-cost way to assemble units of different sizes to a certain total}
\usage{
find_least_cost_partition(desired_dose, available_units, subset_col = NULL,
  subset_val = NULL)
}
\arguments{
\item{desired_dose}{the dose(s) required.}

\item{available_units}{available sizes and costs; must have column names
'size' and 'cost'}

\item{subset_col}{optionally, a column to select on.}

\item{subset_val}{if subset_col is provided, the value to select on in that column.}
}
\value{
a data frame, with columns \code{desired_dose}, \code{used_dose},
\code{waste}, \code{cost}, and \code{cost.no.waste}
}
\description{
Find lowest-cost way to assemble units of different sizes to a certain total
}
\examples{
units <- data.frame(size = c(1000, 250, 200, 100, 50),
    cost = c(40, 11.5, 8.5, 5.5, 4.4))
find_least_cost_partition(450, available_units = units)
temp <- find_least_cost_partition(sample(250:450, 10, replace = TRUE), 
   available_units = units)
}
