
#' Apply a time shift
#' 
#' Shift a survival distribution in time.
#' 
#' @param dist A survival distribution.
#' @param shift A time shift to be applied.
#'   
#' @return A `surv_shift` object.
#' 
#' @details A positive shift moves the fit backwards in time.   That is,
#'   a shift of 4 will cause time 5 to be evaluated as time 1, and so on.
#'   If `shift == 0`, `dist` is returned unchanged.
#' @export
#' 
#' @examples
#' 
#' dist1 <- heemod::define_survival(distribution = "gamma", rate = 0.25, shape = 3)
#' shift_dist <- apply_shift(dist1, 4)
#' heemod::compute_surv(dist1, 1:10)
#' heemod::compute_surv(shift_dist, 1:10)
apply_shift = function(dist, shift) {
  
  stopifnot(
    length(shift) == 1,
    is.finite(shift)
  )
  if(shift == 0) return(dist)
  if(inherits(dist, "surv_shift")){
      dist$shift <- dist$shift + shift
      if(dist$shift == 0) return(dist$dist)
      else return(dist)
  }  
  structure(
      list(
        dist = dist,
        shift = shift
      ),
      class = "surv_shift"
    )
}

#' Plot general survival models
#'
#' @param x a survival object of class `surv_aft`, `surv_add_haz`,
#'   `surv_ph`, `surv_po`, `surv_model`, `surv_pooled`, or `surv_projection`.
#' @param times Times at which to evaluate and plot the survival object.
#' @param type either `surv` (the default) or `prob`, depending on whether
#'   you want to plot survival from the start or conditional probabilities.
#' @param join_col,join_pch,join_size graphical parameters for points
#'   marking points at which different survival functions are joined.
#' @param ... additional arguments to pass to `ggplot2` functions.
#'   
#' @details The function currently only highlights join points that are at
#'   the top level; that is, for objects with class `surv_projection`.
#'   
#'   To avoid plotting the join points, set join_size to a negative number.  
#'
#' @return a [ggplot2::ggplot()] object.
#' @export
#'
plot.surv_obj <- function(x, times, type = c("surv", "prob"), 
                          join_col = "red", join_pch = 20,
                          join_size = 3, ...){
  type <- match.arg(type)
  y_ax_label <- c(surv = "survival", prob = "probability")[type]
  res1 <- data.frame(times = times,
                     res = heemod::compute_surv(x, times, ..., type = type))
  
  this_plot <- 
    ggplot2::ggplot(res1, ggplot2::aes_string(x = "times", y = "res")) + 
    ggplot2::geom_line() + 
    ggplot2::scale_x_continuous(name = "time") + 
    ggplot2::scale_y_continuous(name = y_ax_label)
  
  if("at" %in% names(x))
    this_plot <- this_plot +
    ggplot2::geom_point(data = dplyr::filter_(res1, ~ times == x$at),
                        ggplot2::aes_string(x = "times", y = "res"),
                        pch = "join_pch", size = "join_size", 
                        col = "join_col")
  
  this_plot
  
}

plot.surv_projection <- plot.surv_obj
plot.surv_ph <- plot.surv_obj
plot.surv_add_haz <- plot.surv_obj
plot.surv_model <- plot.surv_obj
plot.surv_po <- plot.surv_obj
plot.surv_aft <- plot.surv_obj
plot.surv_pooled <- plot.surv_obj
plot.surv_shift <- plot.surv_obj


#' Summarize surv_shift objects
#'
#' @param object a `surv_shift` object 
#' @param summary_type "standard" or "plot" - "standard"
#'   for the usual summary of a `survfit` object,
#'   "plot" for a fuller version
#' @param ... other arguments
#' 
#' @return a data frame with columns time, est, lcl, ucl
#'
summary.surv_shift <- 
  function(object, summary_type = c("plot", "standard"), ...){
    summary_type <- match.arg(summary_type)
    res <- summary(object$dist, ...)
    if(inherits(res, "summary.survfit")){
      if(summary_type == "plot"){
        res <- data.frame(res[c("time", "surv", "lower", "upper")])
        names(res) <- c("time", "est", "lcl", "ucl")
      }
    }
    if(length(res) == 1) res <- res[[1]]
    res$time <- res$time + object$shift
    res
    }

summary.surv_projection <-
  function(object, summary_type = c("plot", "standard"), 
           type = c("survival", "cumhaz"), ...){
    summary_type <- match.arg(summary_type)
    res1 <- summary(object$dist1, ...)
    res2 <- summary(object$dist2, ...)
    if(inherits(res1, "summary.survfit")){
      res1 <- data.frame(res1[c("time", "surv", "lower", "upper")])
      names(res1) <- c("time", "est", "lcl", "ucl")
    }
    if(inherits(res2, "summary.survfit")){
      res2 <- data.frame(res1[c("time", "surv", "lower", "upper")])
      names(res2) <- c("time", "est", "lcl", "ucl")
    }
    surv1_p_at <- heemod::compute_surv(
      object$dist1,
      time = object$at 
    )
    surv2_p_at <- heemod::compute_surv(
      object$dist2,
      time = object$at,
      .internal = TRUE)
    res1 <- res1[!is.na(res1$est),]
    res2 <- res2[!is.na(res2$est),]
    res1 <- dplyr::filter_(res1, ~time < object$at)
    res2 <- dplyr::filter_(res2, ~time >= object$at)
    time_col <- match("time", names(res2))
    res2[, -time_col] <-
      res2[, -time_col] * surv1_p_at/surv2_p_at
    res <- dplyr::bind_rows(res1, res2) %>%
      dplyr::arrange_(~time)
    if(type == "cumhaz")
      res[, -time_col] <- -log(res[, -time_col])
    res
          }
  
#' Convert an age-based table to a cycle-based table
#'
#' @param x a survival object
#' @param age_init initial age the table will be applied to
#' @param cycle_length length of cycles
#'
#' @return a table with `time` renormalized to be in terms of cycles,
#'   with time 0 corresponding to age `age_init`
#' @details This is used to convert a mortality table for use as a
#'   survival curve in terms of cycles of the model.  Typically,
#'   age_surv will have been constructed using [define_surv_table()]
#'   or [define_wtd_surv_table()].
#' @export
#'
#' @examples
#' 
#'  df <- data.frame(time = c(50, 51, 55, 60), survival = c(1, 0.9, 0.7, 0.5))
#'  st <- define_surv_table(df, by_age = TRUE)
#'  age_to_time(st, age_init = 50, cycle_length = 1)

age_to_time <- function(x, age_init, cycle_length){
  UseMethod("age_to_time")
}

#' @rdname age_to_time
#' @export
age_to_time.surv_table <- function(x, age_init, cycle_length){
  if(!inherits(x, "data.frame"))
    stop("first argument must be a data frame.")
  if(!("time" %in% names(x)))
    stop("'time' must be one of the names of the first argument.")
  if(!("survival" %in% names(x)))
    stop("'survival' must be one of the names of the first argument.")
  if(age_init < 0)
    stop("age_init has value ",
         age_init,
         " ; but negative age makes no sense")
  if(cycle_length <= 0)
    stop("cycle_length has value ",
         cycle_length,
         " ; but zero or negative cycle length makes no sense")
  x$time <- x$time - age_init
  x <- x[x$time >= 0, ]
  new_time <- 
    seq(from = min(x$time), 
        to = max(x$time), 
        by = cycle_length)
  evenly_spaced_hazard <- 
    stats::approx(x$time, log(x$survival),
         xout = new_time)
  res <- data.frame(time = seq(along = new_time) - 1, 
                    survival = exp(evenly_spaced_hazard$y))
  class(res) <- class(x)
  res
}

#' @rdname age_to_time
#' @export
age_to_time.surv_object <- function(x, age_init, cycle_length){
  pieces <- grep("dist", names(x))
  if(length(pieces))
    x[pieces] <-
      lapply(x[pieces], function(y){age_to_time(y, age_init, cycle_length)})
  x
}
#' @rdname age_to_time
#' @export
age_to_time.surv_pooled <- 
  function(x, age_init, cycle_length){
    x$dists <- lapply(x$dists,
                      age_to_time,
                      age_init = age_init,
                      cycle_length = cycle_length)
    x
  }

