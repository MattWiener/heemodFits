% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_define.R
\name{define_wtd_surv_table}
\alias{define_wtd_surv_table}
\alias{define_wtd_surv_table.character}
\alias{define_wtd_surv_table.data.frame}
\title{combine multiple columns into a survival table}
\usage{
define_wtd_surv_table(x, time_col, weights, by_age)

\method{define_wtd_surv_table}{character}(x, time_col, weights, by_age)

\method{define_wtd_surv_table}{data.frame}(x, time_col, weights, by_age)
}
\arguments{
\item{x}{a data frame, or a string pointing to such a data frame}

\item{time_col}{the column of x measuring time}

\item{weights}{a named vector of weights.   Names correspond to
the columns of \code{x} to which the weights should be applied.}

\item{by_age}{is time in the table measured by age instead of cycle
(as is typical for a mortality table, for example).}
}
\value{
a \code{surv_pooled} object (see \code{\link[=mix]{mix()}}).
}
\description{
combine multiple columns into a survival table
}
\details{
\code{weights} are not checked for adding to 1.
}
\examples{
df <- data.frame(age = c(50, 55, 60, 65), 
                 male = c(1, 0.9, 0.8, 0.7),
                 female = c(1, 0.8, 0.7, 0.6)
                 )
wtd_table <- 
   define_wtd_surv_table(df, time_col = "age",
                         weights = c(male = 0.52, female = 0.48),
                         by_age = TRUE)                
}
