H_piecewise_multi <- function(t, pieces_def, dists, log = FALSE){
  ## pieces_def must be a data frame with columns start,
  ##  dist.name, and and columns for the arguments for the 
  ##  hazard functions for those distributions
  ## (see r_piecewise_exponential_builder)
  if(!is.data.frame(pieces_def)){
    pieces_def <- convert_vector_to_struct(pieces_def, dists)
    }
  pieces_def$start2 <- pieces_def$start
  start_times <- look_up(pieces_def, start = t, #pieces_def$start,
                                    bin = "start", value = "start2")
  indices <- match(start_times, pieces_def$start)
  pieces_def$Haz_fn <- paste0("H", pieces_def$dist_name)
  pieces_def$finish <- c(pieces_def$start[-1], Inf)
  
  start_data <- data.frame(x = t, 
                           pieces_def[indices,])

  ## there is probably an easier way to do this with dplyr or purrr
  collect_cum_hazards <- function(start_data){
    res <- matrix(0, nrow = nrow(start_data), 
                  ncol = nrow(pieces_def))
    
    for(this_piece in 1:nrow(pieces_def)){
      start_data$piece_start <- 
        pieces_def[this_piece, "start"]
      ## if t is in an earlier piece, finish at the start of this piece
      ##    (so no contribution from this piece)
      ## if t is in this piece, finish at t;
      ## if t is in a later piece, finish at the end of this piece
      ##    (so get the entire cumulative hazard from this piece)
      start_data$piece_end <- 
        pmax(pieces_def[this_piece, "start"], 
             pmin(t, pieces_def[this_piece, "finish"]))
      this_Haz <- pieces_def[this_piece, "Haz_fn"]
      ## construct argument lists 
      use_arg_names <- setdiff(names(formals(get(this_Haz))), 
                               c("log", "x"))
      base_args <-   pieces_def[this_piece, use_arg_names]
      names(base_args) <- use_arg_names
      start_arg_list <- c(base_args, x = list(start_data[, "piece_start"]))
      finish_arg_list <- c(base_args, x = list(start_data[, "piece_end"]))
      ## and call the Hazard functions
      res[, this_piece] <- 
        do.call(this_Haz, args = finish_arg_list) - 
          do.call(this_Haz, args = start_arg_list)
         }
    rowSums(res)
  }
  if(nrow(start_data) > 0)
    res <- collect_cum_hazards(start_data)
  else
    res <- numeric(0)
  if(log == TRUE) res <- log(res)
  return(res)
}

h_piecewise_multi <- function(t, pieces_def, dists, log = FALSE){
  ## pe.def must be a data frame with columns start and hazard
  ## (see r_piecewise_exponential_builder)
     if(!is.data.frame(pieces_def)){
       pieces_def <- convert_vector_to_struct(pieces_def, dists)
       }
  pieces_def$start2 <- pieces_def$start
  start_times <- look_up(pieces_def, start = t,
                         bin = "start", value = "start2")

  indices <- match(start_times, pieces_def$start)
  pieces_def$haz_fn <- paste0("h", pieces_def$dist_name)  
  
    start_data <- data.frame(x = t, 
                           pieces_def[indices,])
  
  ## there is probably an easier way to do this with dplyr or purrr
  collect_hazards <- function(start_data){
    res <- numeric(nrow(start_data))
    for(this_start_time in unique(start_data$start)){
      indices <- which(this_start_time == start_data$start)
      this_haz_fn <- unique(start_data[indices, "haz_fn"])
      use_arg_names <- setdiff(names(formals(get(this_haz_fn))), "log")
      res[indices] <- 
        do.call(this_haz_fn, args = start_data[indices, use_arg_names])
      }
    res
  }
  res <- collect_hazards(start_data)
  
  if(log == TRUE) res <- log(res)
  return(res)
}

convert_vector_to_struct <-
  function(pieces_def_vec, dists){
    num_rows <- length(dists)
    break_points <- as.numeric(pieces_def_vec[1:num_rows])
    dist_args <- as.numeric(pieces_def_vec[-c(1:num_rows)])
    build_pieces_def_structure(dists, break_points, dist_args)
  }


#' Build up general piecewise distribution specification
#'   from the component pieces
#'
#' @param dists Character vector naming the distributions for each piece.
#' @param breaks Breakpoints at which distributions change.
#' @param params parameters of the hazard functions for the distributions.
#'
#' @return  A data frame with columns \code{breaks}, \code{dist_name},
#'   and one additional column for each argument (other than \code{x}
#'   or \code{log}) of the arguments to the hazard functions for the
#'   distributions.
#' @examples
#' \dontrun{
#' ## requires package flexsurv to be loaded
#' build_pieces_def_structure(c("exp", "weibull", "exp"),
#'                            c(0, 100, 200),
#'                            c(0.005, 1, 100, 0.002))
#'                            }
build_pieces_def_structure <-
  function(dists, breaks, params) {
    ## dists is a character vector of distribution types,
    ##   one for each piece
    hfn_list <- haz_fn_args(dists)
    hfn_args <- hfn_list$hfn_args
    num_hfn_args <- hfn_list$num_hfn_args
    cum_hfn_args <- c(0, cumsum(num_hfn_args))
    ## params:  1 for number of pieces, then specify the distributions,
    ##          then all the arguments for the distributions
    num_params <- length(dists) + sum(num_hfn_args)
    
    arg_names <- unique(unlist(hfn_args))
    pieces_def <- data.frame(start = breaks,
                             dist_name = dists,
                             matrix(NA, nrow = length(dists),
                                    ncol = length(arg_names)),
                             stringsAsFactors = FALSE
    )
    names(pieces_def)[-c(1,2)] <- arg_names
    ## put the parameters in the correct places in
    ##   the pieces_def data frame
    for(i in seq(along = dists)){
      param_ind <- cum_hfn_args[i] + 1:num_hfn_args[i]
      arg_cols <- match(hfn_args[[i]], names(pieces_def))
      pieces_def[i, arg_cols] <- params[param_ind]
    }
    pieces_def
  }


#' Construct piecewise distribution function with use with survreg
#'
#' @param dists The distributions for the consecutive pieces
#' @param fixed_bp Should breakpoints be considered fixed?
#'
#' @return A custom distribution function for use with survreg.  We don't
#'   include a function to determine initial values, which will be
#'   provided in the call to flexsurvreg.
#'
get_custom_fss_piecewise_multi <-
  function(dists, fixed_bp) {
   
    hfn_list <- haz_fn_args(dists)
    hfn_args <- hfn_list$hfn_args
    num_hfn_args <- hfn_list$num_hfn_args
    ## params:  1 for number of pieces, then specify the distributions,
    ##          then all the arguments for the distributions
    num_params <- length(dists) + sum(num_hfn_args)
    custom_fss <-
      list(
        name = "piecewise_multi",
        pars = paste0("pieces_def", 0:(num_params - 1)),
        location = paste0("pieces_def",  length(dists))
      )
    
    if (fixed_bp) {
      transform_nums <- c(length(dists), sum(num_hfn_args))
          }
    else{
      transform_nums <- c(1, length(dists) + sum(num_hfn_args) - 1)
         }
    custom_fss$transforms <- rep(c(identity, log), transform_nums)
    custom_fss$inv.transforms <- rep(c(identity, exp), transform_nums)
    custom_fss
  }

haz_fn_args <- function(dists){
  hfns <- paste0("h", dists)
  hfn_args <- lapply(hfns, function(this_fn){
    setdiff(names(formals(this_fn)), c("x", "log"))
  }
  )
  num_hfn_args <- sapply(hfn_args, length)
  list(hfns = hfns, 
       hfn_args = hfn_args, 
       num_hfn_args = num_hfn_args)
}

get_multi_dist_inits <- function(these_bp, dists, t){
  hfn_list <- haz_fn_args(dists)$hfns
  ## go to flexsurv::flexsurv.dists to get initialization functions
  flexsurv_dists <- match(dists, names(flexsurv::flexsurv.dists))
  init_vals <- 
    lapply(seq(along = dists), 
         function(i){
           init_fn <- flexsurv::flexsurv.dists[[flexsurv_dists[i]]]$inits
           if("aux" %in% names(formals(init_fn)))
             init_fn(t, aux = list(counting = TRUE))
           else
             init_fn(t)
         }
    )
  c(these_bp, unlist(init_vals))
  }
  
  

#' Create a random number generator for a piecewise exponential distribution
#'
#' @param piecewise_def a data frame with columns \code{start} and
#'   \code{dist_name}, and additional columns for parameters of distributions.
#'     \code{dist_name[i]} along with its parameters define the hazard 
#'     from \code{start[i]} to
#'   \code{start[i+1]} (or all remaining time in the case of the 
#'   last distribution).
#'
#' @return A function that takes in a number \code{n} and returns
#'    \code{n} random draws from the piecewise distribution
#'    defined by \code{piecewise_def}
#' @export
#'
#' @examples
#' my_pe <- data.frame(start = c(0, 200), hazard = c(0.005, 0.025))
#' r_my_pe <- r_piecewise_exponential_builder(my_pe)
#' plot(survival::survfit(survival::Surv(r_my_pe(1000), rep(1, 1000)) ~ 1), 
#'      fun = "cumhaz")
#' abline(0, 0.005, col = "red")
#' abline(-200 * (0.025 - 0.005), 0.025, col = "red", lty = 2)
r_piecewise_builder <- 
  function(piecewise_def){
    
    eval_dist <- 
      function(x){
        dist_name <- x$dist_name
        args <- x[-match(c("start", "dist_name"), names(x))]
        names(args) <- setdiff(names(x), c("start", "dist_name"))
        dont.drop <- !is.na(args)
        name.vec <- names(args)[dont.drop]
        args <- args[dont.drop]
        names(args) <- name.vec
        function(x, letter){
          fn_name <- paste0(letter, dist_name)
          new_args <- as.list(args)
          old_names <- names(new_args)
          new_args[[length(new_args) + 1]] <- "x"
          new_args[[length(new_args)]] <- x
          names(new_args)[length(new_args)] <- formals(get(fn_name))[1]
          do.call(fn_name, new_args)}
      }
    
    eval_fns <- lapply(1:nrow(piecewise_def),
                       function(this_row){
                         eval_dist(piecewise_def[this_row,])
                       }
    )
    cum_hazard_at_start <- 
      lapply(1:nrow(piecewise_def),
             function(i){
               eval_fns[[i]](piecewise_def[i, "start"], "H")
             }
      )
    
    last_time <- eval_fns[[length(eval_fns)]](0.9999999, "q")
    last_haz <- eval_fns[[length(eval_fns)]](last_time, "H")
    
    cum_hazard <- cumsum(c(cum_hazard_at_start, last_haz))
    time <- c(piecewise_def$start, last_time)
    browser()
    ## hazards are not constant, so can't use linear interpolation,
    ##   when distribution is not exponential.
    this_fun <-  stats::approxfun(x = -cum_hazard, 
                                  y = time)
    function(n){this_fun(log(stats::runif(n)))}
  }
