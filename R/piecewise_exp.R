H_piecewise_exp <- function(t, pe_def, log = FALSE){
  ## pe.def must be a data frame with columns start and hazard
  ## (see r_piecewise_exponential_builder)
  ## we're assuming no interaction with covariates here
 if(!is.data.frame(pe_def)){
   use_nrow <- length(pe_def)/2
   pe_def <- data.frame(start = pe_def[1:use_nrow],
                        hazard = pe_def[use_nrow + 1:use_nrow])
 }

  pe_def$cum_hazard <- 
    c(0, cumsum(diff(pe_def$start) * pe_def$hazard[-length(pe_def$hazard)]))

  ## because start itself will get turned into a factor
  pe_def$start2 <- pe_def$start 

  cum_until_last <- heemod::look_up(pe_def, start = t, 
                            bin = "start", value = "cum_hazard")
  hazard_after_last <- heemod::look_up(pe_def, start = t,
                               bin = "start", value = "hazard")
  time_after_last <- t - heemod::look_up(pe_def, start = t,
                                 bin = "start", value = "start2")
  res <- cum_until_last + time_after_last * hazard_after_last

  if(log == TRUE) res <- log(res)
  res
}

h_piecewise_exp <- function(t, pe_def, log = FALSE){
  ## pe.def must be a data frame with columns start and hazard
  ## (see r_piecewise_exponential_builder)
     if(!is.data.frame(pe_def)){
      use_nrow <- length(pe_def)/2
      pe_def <- data.frame(start = pe_def[1:use_nrow],
                           hazard = pe_def[use_nrow + 1:use_nrow])
     }
   # print(pe_def)
    res <- heemod::look_up(pe_def, start = t, 
                 bin = "start", value = "hazard")
  if(log == TRUE) res <- log(res)
  return(res)
}

#' Construct piecewise exponential distribution function with use with survreg
#'
#' @param fixed_bp Should breakpoints be considered fixed?
#' @param nr number of rows in the parameter matrix
#' @param nc number of columns in the parameter matrix
#'
#' @return A custom distribution function for use with survreg.  We don't
#'   include a function to determine initial values, which will be
#'   provided in the call to flexsurvreg.

get_custom_fss_piecewise_exp <-
  function(fixed_bp, nr, nc) {
    
    
    custom_fss <-
      list(
        name = "piecewise_exp",
        pars = c(paste0("pe_def", 0:(nr * nc - 1))),
        location = c("pe_def0")
      )
    
    if (fixed_bp) {
      custom_fss$transforms <- rep(c(identity, log), each = nr)
      custom_fss$inv.transforms <- rep(c(identity, exp), each = nr)
    }
    else{
      custom_fss$transforms <- c(identity, rep(c(log), nr * nc - 1))
      custom_fss$inv.transforms <-
        c(identity, rep(c(exp), nr * nc - 1))
    }
    custom_fss
  }

#' Create piecewise exponential definition
#'
#' @param fixed_bp Should breakpoints be considered fixed?  (Otherwise optimized over.)
#' @param used_bp breakpoints
#' @param best_fit the fit corresponding to those breakpoints
#'
#' @return a piecewise exponential specification, a data frame with
#'   columns start (the breakpoints) and hazard.
#'
get_best_pe <- function(fixed_bp, used_bp, best_fit) {
  nr <- length(used_bp)
  if (fixed_bp) {
    best_pe_def <- data.frame(start = used_bp,
                              hazard = exp(best_fit$opt$par[1:nr]))
  }
  else
    best_pe_def <-
      data.frame(start = cumsum(c(0, exp(
        best_fit$opt$par[1:(nr - 1)]
      ))),
      hazard = exp(best_fit$opt$par[nr - 1 + 1:nr]))
  rownames(best_pe_def) <- NULL
  best_pe_def
}


#' Create a random number generator for a piecewise exponential distribution
#'
#' @param piecewise_exponential a data frame with columns \code{start} and
#'   \code{hazard}.  \code{hazard[i]} is the hazard from \code{start[i]} to
#'   \code{start[i+1]} (or all remaining time in the case of the 
#'   last \code{hazard}).
#'
#' @return A function that takes in a number \code{n} and returns
#'    \code{n} random draws from the piecewise exponential distribution
#'    defined by \code{piecewise_exponential}
#' @export
#'
#' @examples
#' my_pe <- data.frame(start = c(0, 200), hazard = c(0.005, 0.025))
#' r_my_pe <- r_piecewise_exponential_builder(my_pe)
#' plot(survival::survfit(survival::Surv(r_my_pe(1000), rep(1, 1000)) ~ 1), 
#'      fun = "cumhaz")
#' abline(0, 0.005, col = "red")
#' abline(-200 * (0.025 - 0.005), 0.025, col = "red", lty = 2)
r_piecewise_exponential_builder <- 
  function(piecewise_exponential){
    
    cum.hazard <- 
      piecewise_exponential$start *piecewise_exponential$hazard
    
    ## get the time corresponding to cumulative hazard of 100
    ##   (giving effectively 0 survival probability)
    last_cum_hazard <- cum.hazard[length(cum.hazard)]
    last_exp_haz <- 
      piecewise_exponential$hazard[length(piecewise_exponential$hazard)]
    last_time <- (100 - last_cum_hazard)/last_exp_haz
    
    if(piecewise_exponential$start[length(piecewise_exponential$start)] < last_time){
      piecewise_exponential <- 
        data.frame(start = c(piecewise_exponential$start, last_time),
                   hazard = c(0,piecewise_exponential$hazard)
        )
      cum.hazard <- 
        cumsum(piecewise_exponential$start *piecewise_exponential$hazard)
    }
    ## get a function that takes negative cumulative hazard
    ##    (log survival) to times; so it gets fed the logarithm
    ##    of uniform (0,1)
    this_fun <-  stats::approxfun(x = -cum.hazard, 
                                  y =piecewise_exponential$start)
    function(n){this_fun(log(stats::runif(n)))}
  }


