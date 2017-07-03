
#' @rdname eval_surv
#' @export
eval_surv.surv_shift <- function(x, time, ...) {
  
  time_ <- time
  time_ <- time_ - x$shift
  ret <- rep(1, length(time_))
  keep_me <- time_ >= 0
  if(any(keep_me)){
    time_ <- time_[keep_me]
    ##check_cycle_inputs(time_, cycle_length)
    ret[keep_me] <- eval_surv(
      x$dist,
      time = time_,
      ...
    ) 
  }

  ret
}

#' @rdname eval_surv
#' @export
eval_surv.surv_dist <- function(x, time, ...) {

  if (! requireNamespace("flexsurv")) {
    stop("'flexsurv' package required.")
  }
  
  pf <- get(paste0("p", x$distribution),
            envir = asNamespace("flexsurv"))
  
  args <- x[- match("distribution", names(x))]
  args[["q"]] <- time
  args[["lower.tail"]] <- FALSE
  ret <- do.call(pf, args)
  
  ret
}

#' @rdname eval_surv
#' @export
eval_surv.surv_table <- function(x, time, ...){
  heemod::look_up(data = x, time = time, bin = "time", value = "survival")
}

eval_surv.lazy <- function(x, ...){
  dots <- list(...)
  use_data <- list()
  if("extra_env" %in% names(dots))
    use_data <- as.list.environment(dots$extra_env)
  eval_surv(lazyeval::lazy_eval(x, data = use_data), ...)
}

eval_surv.character <- function(x, ...){
  eval_surv(eval(parse(text = x)), ...)
}