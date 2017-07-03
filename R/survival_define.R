
#' Define a survival distribution based on explicit survival probabilities
#'
#' @param x a data frame with columns `time` and `survival` 
#' @param by_age is time actually age?  If so, don't require time to start wtih 0.
#'
#' @return a `surv_table` object, which can be used with [heemod::compute_surv()].
#' @export
#'
#' @examples
#'  x <- data.frame(time = c(0, 1, 5, 10), survival = c(1, 0.9, 0.7, 0.5))
#'  define_surv_table(x)
#'  
define_surv_table <- function(x, by_age = FALSE){
  UseMethod("define_surv_table")
}

#' @rdname define_surv_table
#' @export
define_surv_table.data.frame <- function(x, by_age = FALSE){
  required_names <- c("time", "survival")
  names_present <- required_names %in% names(x)
  if(any(!names_present)){
    stop("missing column",
         plur(sum(!names_present)),
         " in surv_table object: ",
         paste(required_names[!names_present], collapse = ", ")
    )
  }
  x$time <- as.numeric(x$time)
  x <- x[order(x$time),]
  dup_time <- duplicated(x$time)
  if(any(dup_time))
    stop("any time can appear only once in explicit survival data. ",
         "Duplicated time",
         plur(sum(dup_time)),
         ": ",
         paste(x$time[dup_time], collapse = ", ")
    )
  
  if(x$survival[1] != 1)
    stop("surv_table data must start with survival 1")
  if(x$time[1] != 0 & !by_age)
    stop("surv_table data must start with time 0 (unless by_age is TRUE)")
  increasing_survival <- diff(x$survival) > 0
  if(any(increasing_survival)){
    problem_times <- matrix(x$time[which(increasing_survival) + c(0,1)],
                            ncol = 2, byrow = TRUE)
    stop("survival cannot increase over time; see times:\n",
         paste("(", 
               problem_times[,1],
               ", ",
               problem_times[,2],
               ")",
               sep = "", collapse = ", ")
    )
  }
  class(x) <- c("surv_table", "surv_object", "data.frame")
  x
}
#' @rdname define_surv_table
#' @export
define_surv_table.character <- function(x, by_age = FALSE){
  define_surv_table(read_file(x), by_age)
}

#' combine multiple columns into a survival table
#'
#' @param x a data frame, or a string pointing to such a data frame
#' @param time_col the column of x measuring time
#' @param weights a named vector of weights.   Names correspond to
#'   the columns of `x` to which the weights should be applied.
#' @param by_age is time in the table measured by age instead of cycle
#'   (as is typical for a mortality table, for example).
#' @details `weights` are not checked for adding to 1.
#' @return a `surv_pooled` object (see [mix()]).
#' @export
#'
#' @examples
#' df <- data.frame(age = c(50, 55, 60, 65), 
#'                  male = c(1, 0.9, 0.8, 0.7),
#'                  female = c(1, 0.8, 0.7, 0.6)
#'                  )
#' wtd_table <- 
#'    define_wtd_surv_table(df, time_col = "age",
#'                          weights = c(male = 0.52, female = 0.48),
#'                          by_age = TRUE)                
define_wtd_surv_table <- function(x, time_col, weights, by_age){
  UseMethod("define_wtd_surv_table")
}
#' @rdname define_wtd_surv_table
#' @export
define_wtd_surv_table.character <- function(x, time_col, weights, by_age){
  define_wtd_surv_table(read_file(x), time_col, weights, by_age)
}

#' @rdname define_wtd_surv_table
#' @export
define_wtd_surv_table.data.frame <- function(x, time_col, weights, by_age){
  if(!(time_col %in% names(x)))
    stop("column '",
         time_col, 
         "' not in table")
  if(time_col == "age" & !by_age)
    warning("time_col is 'age', but 'by_age' is FALSE - is that intentional?")
  other_cols <- names(weights)
  missing_cols <- !(other_cols %in% names(x))
  if(any(missing_cols))
    stop("column",
         plur(sum(missing_cols)),
         " ",
         paste(other_cols[missing_cols], collapse = " ,"),
         " specified in weights but not present in table"
         )
  individual_tables <- 
    lapply(other_cols, function(this_col){
      y <- x[, c(time_col, this_col)]
      names(y) <- c("time", "survival")
      define_surv_table(y, by_age)
    }
  )
  heemod::mix_(individual_tables, weights = weights)
  }
  