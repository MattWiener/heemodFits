#' Allow assignment of a cost at a particular cycle or age
#'
#' @param cost cost, in units when incurred
#' @param age age at a given cycle
#' @param cycle the cycle
#' @param at_age at what age should the cost be incurred?
#' @param at_cycle in what cycle should cost be incurred?
#'
#' @return  0 except when the age or cycle condition is met; 
#'   \code{cost} when at least one condition is met.
#' @export
#'
#' @details intended for use in parameters object 
#'   for \code{heemod} models.
#' 
#' @examples
#' cost_at_right_time(500, age = 1:20, at_age = 10)
#' cost_at_right_time(500, cycle = 1:20, at_cycle = 10)
cost_at_right_time <-
  function(cost, age, cycle, at_age, at_cycle){
    if(missing(cost))
      stop("'cost' must be specified")
    if(missing(at_age) & missing(at_cycle))
      stop("must specify at least one of 'at_cost' and 'at_age'")
    if(missing(cycle) & missing(age))
      stop("must specify at least one of 'cost' and 'age'")
    if((!missing(cycle) & !missing(age)) &&
       length(cycle) != length(age))
        stop("if both cycle and age are defined, ",
             "they must have the same length")
    have_cycle_info <- !missing(cycle) & !missing(at_cycle)
    have_age_info <- !missing(age) & !missing(at_age)
    if(!have_cycle_info & !have_age_info)
      stop("must give at least one of both 'cycle' and 'at_cycle' ", 
           "or both 'age' and 'at_age'")
    
    cond_length <- ifelse(missing(cycle), length(age), length(cycle))
            
    if(missing(at_age) | missing(age))
      cost_age <- rep(FALSE, cond_length)
    else
      cost_age <- age %in% at_age
    
    if(missing(at_cycle) | missing(cycle))
      cost_cycle <- rep(FALSE, cond_length)
    else
      cost_cycle <- cycle %in% at_cycle

    return(cost * (cost_age | cost_cycle))
            
  }


#' Calculate utilities based on time before death
#' @param x an object of of type \code{cycle_counts}, 
#'   \code{eval_model}, or \code{run_models}.  x is either
#'   a matrix of counts for from a model.
#'   (of type \code{cycle_counts}) or an \code{eval_model} or
#'   a \code{run_models} list.
#' @param m name or number of a model, only if \code{x} is
#'   of type \code{run_models}.
#' @param death_state the name of the state representing death.
#' @param util_before_death a matrix with two columns, \code{through_cycle} 
#'   and \code{util}.  From \code{through_cycle[i - 1]} (or 1, for i = 1)
#'   to \code{through_cycle[i+1]} before death, 
#'   an individual will accrue \code{util[i]} units of
#'   utility.  
#' @param util_long_before_death utility more than 
#'   \code{max(util_before_death[, through_cycle])} cycles before death.
#' @param discount_rate discount rate to be applied to utilities.   
#' @param ... other arguments to be passed through

#'
#' @return  a vector of utilities, of length \code{nrow(counts)}.
#' @export
#'
#' @examples
#' ProgFree <- round(1000 * exp(-0.4 * 0:24))
#' Progressive <- round((1000 - ProgFree) * exp(-0.4 * 0:24))
#' Death <- 1000 - ProgFree - Progressive
#' state_names <- rep(c("ProgFree", "Progressive", "Death"), each = 25)
#' 
#' counts <- data.frame(.strategy = rep("s1", 25),
#'                      markov_cycle = 0:24,
#'                      state_names = state_names,
#'                      count = c(ProgFree, Progressive, Death)
#' )
#' class(counts) <- c("cycle_counts", class(counts))
#' aa1 <- data.frame(until_lag = c(1,2,3), util = c(0.1, 0.3, 0.5))
#'  utility_by_time_from_death(counts, util_before_death = aa1, 
#'    util_long_before_death = 1)

utility_by_time_from_death <- 
  function(x, ...){
    UseMethod("utility_by_time_from_death")
  }

#' @export
#' @rdname utility_by_time_from_death
utility_by_time_from_death.cycle_counts <- 
  function(x, util_before_death,
           util_long_before_death, death_state = "Death",
           discount_rate = 0, ...){
    counts <- x
    stopifnot(discount_rate >= 0)
    stopifnot(is.character(death_state))
    if(! death_state %in% counts$state_names)
      stop("death_state must be the name of one of the states")
    if(!all(c("until_lag", "util") %in% names(util_before_death)))
      stop("util_before_death must have columns until_lag and util")
    subjects_by_cycle <- 
      counts %>% 
      dplyr::group_by_( ~ markov_cycle) %>%
      dplyr::summarize_(num_subjects = ~sum(count))
    num_subjects <- unique(round(subjects_by_cycle$num_subjects, 6))
    stopifnot(length(num_subjects) ==  1)
    if(any(util_before_death[, "until_lag"] <= 0))
      stop("problem with util_before_death: can't specify values for Markov cycles <= 0")
    if(max(util_before_death[, "until_lag"]) < 1)
      stop("problem with util_before_death:  must specify utility for some cycle > 1")
    
    death_counts <- 
      dplyr::filter_(counts, ~state_names == death_state)$count
    if(max(death_counts) != num_subjects){
      stop("not all subjects reach the death state: (", 
           round(max(death_counts), 3),
           " out of ",
           num_subjects,
           "); \n",
           "to compute utility by time before death, ",  
           "you need to know time of death for each subject.")
    }
    alive_counts_by_cycle <- 
      counts %>% 
      dplyr::filter_( ~ state_names != death_state) %>%
      dplyr::group_by_(~ markov_cycle) %>%
      dplyr::summarize_(alive_counts = ~ sum(count))
    
    alive_counts <- alive_counts_by_cycle$alive_counts
    new_deaths <- c(0, diff(death_counts))
    discounted_new_deaths <- 
      heemod::discount(new_deaths, discount_rate, first = FALSE) 
    
    ## feels like there should be a better way to do this by having
    ## util_long_before_death instead 0 at the end of the modified
    ## util_df below.   but then miss the long before death utils.
    util_df <- 
      data.frame(cycles = c(0, util_before_death[, "until_lag"] + 1e-5),
                 util = c(util_before_death[, "util"], 0))
    
    ## utilities associated with different times before death
    util_vec <- heemod::look_up(
      util_df,
      cycles = 1:nrow(counts) - 1,
      bin = "cycles",
      value = "util"
    )
    ## filter vector to just count up the people who are
    ##   the specified times before death
    use_length <- 
      min(max(util_before_death[, "until_lag"] + 1),
          length(util_vec))
    
    near_death_vec <- rep(1, max(util_before_death[, "until_lag"]) + 1)
    util_vec <- util_vec[pmin(1:length(near_death_vec),
                              length(util_vec))]
    ## adjust for lag
    util_vec[1] <- near_death_vec[1] <- 0        
    
    
    padding <- rep(0, length(near_death_vec))
    drop_padding <- -c(1:length(near_death_vec))
    num_near_death <- 
      stats::filter(c(padding, rev(new_deaths)), 
                    near_death_vec, sides = 1)
    num_near_death <- rev(num_near_death[drop_padding])
    
    utils_before_death <- 
      stats::filter(c(padding, rev(discounted_new_deaths)), 
                    util_vec, sides = 1)
    utils_before_death <- 
      rev(utils_before_death[drop_padding])
    earlier_utils <- (alive_counts - num_near_death) * util_long_before_death
    
    utils_before_death + earlier_utils
    
  }

#' @rdname utility_by_time_from_death
#' @export
utility_by_time_from_death.eval_strategy <- 
  function(x, ...){
    utility_by_time_from_death(heemod::get_counts(x), ...)
  }

#' @rdname utility_by_time_from_death
#' @export
utility_by_time_from_death.run_model <- 
  function(x, m, ...){
    temp <- dplyr::filter_(heemod::get_counts(x), ~ .strategy_names == m)
    utility_by_time_from_death(temp, ...)
  }

utility_by_time_from_death.data.frame <-
  function(x, ...){
    x$markov_cycle <- 1:nrow(x)
    tidied <- tidyr::gather(x, key = "state_names", 
                            value = "count", -dplyr::matches("markov_cycle"))
    utility_by_time_from_death.cycle_counts(tidied, ...)
  }

#' Model piecewise linear decline in vaccine efficacy
#' or other variables
#'
#' @param marked_cycles cycles at which effectiveness
#'   fractions are specified.
#' @param marked_values effectiveness fractions
#'   at the corresponding cycles.
#' @param init_value a point at (0, init_value) will be
#'  added
#' @param change_per_cycle can be specified instead of
#'   \code{marked_cycles} and \code{marked_times}.   If both
#'   are specified, \code{marked_cycles} and \code{marked_values} 
#'   take precedence.
#' @param max_change_cycles used with change_per_cycle
#' @param age_start_change_per_cycle at what age should the decline in
#'   efficacy start (if not specified, there's no threshold).
#' @param cycle_length used to calculate age from cycles if necessary
#' @param cycles The cycles for which values should be returned.  
#'   If \code{NULL}, the default, a function is returned.
#' @param age_at_cycle_0 used in conjunction with
#'   \code{age_start_change_per_cycle}.
#' @param min_val minimum value the interpolation can take on.
#'
#' @return a function that calculates the decline in efficacy
#'   based on cycle and, if specified, patient age at the
#'   start of a simulation.
#' @details intended for use in parameters object 
#'   for \code{heemod} models.
#' @export
#'
#' @examples
#' linear_interpolation(c(1, 10, 20), c(1, .5, .25))
#' linear_interpolation(change_per_cycle = -0.02, max_change_cycles = 50)
#' linear_interpolation(change_per_cycle = -0.02, max_change_cycles = 100,
#'    min_val = -0.4)
#' linear_interpolation(age_start_change_per_cycle = 65, age_at_cycle_0 = 60,
#'   change_per_cycle = -0.05, cycles = 1:15, max_change_cycles = 5)
#' 
linear_interpolation <- function(marked_cycles,
                                 marked_values,
                                 init_value = 1,
                                 change_per_cycle,
                                 max_change_cycles,
                                 age_start_change_per_cycle = NULL,
                                 cycle_length = 1,
                                 cycles = NULL,
                                 age_at_cycle_0,
                                 min_val = -Inf){
  
  if(missing(marked_cycles) & missing(marked_values)){
    if(missing(change_per_cycle))
      stop(paste("must specify either marked_cycles and marked_values",
                 "or change_per_cycle"))
    else{
      marked_cycles <- max_change_cycles
      marked_values <- init_value + max_change_cycles * change_per_cycle
    }
  }
  
  if(is.character(marked_cycles))
    marked_cycles <- eval(parse(text = marked_cycles))
  if(is.character(marked_values))
    marked_values <- eval(parse(text = marked_values))
  
  if(length(marked_cycles) != length(marked_values))
    stop(paste("marked_cycles and marked_values",
               "must have the same length"))
  
  marked_cycles <- c(0, marked_cycles)
  marked_values <- c(init_value, marked_values)
  
  
  
  decline_fn <- stats::approxfun(x = marked_cycles, 
                          y = marked_values, 
                          rule= c(2,2))
  
  if(is.null(age_start_change_per_cycle))
    final_fn <- function(cycle, start_age = NULL) {
      pmax(decline_fn(cycle), min_val)
      }
  else
    final_fn <- function(cycle, start_age = NULL){
      age_at_time <- cycle * cycle_length + start_age
      pmax(ifelse(age_at_time < age_start_change_per_cycle, 1, 
                  decline_fn((age_at_time - age_start_change_per_cycle - cycle_length)/cycle_length)),
           min_val)
    }
  
  if(is.null(cycles)) return(final_fn)
  else return(final_fn(cycles, age_at_cycle_0))
}

#' Determine whether to dose during certain cycles.
#'
#' @param N cycles or periods to check.
#' @param init Non-repeating dosing indicator at beginning.
#' @param pattern Repeating dosing pattern after initial.
#' @param first if \code{init} is not specified, how many periods
#'   at the beginning should be dosed.  Must be non-negative.
#' @param then_every if \code{pattern} is not specified, make it
#'   (then_every - 1) zeroes followed by a 1.   If \code{then_every}
#'   has any negative value, there is no dosing after the initial
#'   period (equivalent to setting pattern = 0).
#' @param cap largest dosing period
#'
#' @return A logical vector indicating whether or not a patient
#'   is dosed during the relevant period (markov cycle).
#' @export
#' @details \code{init} takes precedence over \code{first}; that is,
#'   if \code{init} is defined, then \code{first} is ignored.  Similarly,
#'   \code{pattern} takes precedence over \code{then_every}.
#' @examples
#'   is_dosing_period(N = 1:13, first = 4, then_every = 3, cap = 40)
#'   is_dosing_period(N = 37:46, first = 4, then_every = 3, cap = 40)
#'   is_dosing_period(N = 1:100, init = c(1,0,1,0,1,0,1,1), pattern = c(1, 0, 1, 1, 0), cap = 120)
#'   ## stop after initial period
#'   is_dosing_period(N = 1:8, first = 4, pattern = 0, cap = 40)
#'   ## demonstrating argument precedence rules
#'   is_dosing_period(N = 1:10, init = c(1,0,1), first = 3, then_every = 5)
#'   is_dosing_period(N = 1:10, init = numeric(0), pattern = c(1, 1, 0, 1, 0), then_every = 2)
#' 
is_dosing_period <- function(N, init, pattern, first, then_every, cap = Inf){
  if(missing(init)){
    if(missing(first)) stop("must specify either init or first")
    if(first < 0) stop("first must be 0 or positive")
    init <- rep(1, first)
  }
  if(missing(pattern)){
    if(missing(then_every)) stop("must specify either pattern or then_every")
    if(then_every < 0) stop("then_every cannot be negative")
    else pattern <- c(rep(0, then_every - 1), 1)
  }
  if(!all(c(init, pattern) %in% c(0,1, TRUE, FALSE)))
    stop("all elements of init and pattern must be FALSE or 0 or TRUE or 1")
  cond1 <- N <= length(init) & as.logical(init[N])
  cond2 <- N > length(init) & N <= cap
  cond3 <- as.logical(pattern[(N - length(init) - 1) %% length(pattern) + 1])
  cond1 | (cond2 & cond3)
}

#' Find lowest-cost way to assemble units of different sizes to a certain total
#'
#' @param desired_dose the dose(s) required.
#' @param available_units available sizes and costs; must have column names
#'   'size' and 'cost'
#' @param subset_col optionally, a column to select on.
#' @param subset_val if subset_col is provided, the value to select on in that column.
#'
#' @return a data frame, with columns `desired_dose`, `used_dose`,
#'  `waste`, `cost`, and `cost.no.waste`
#' @export
#'
#' @examples
#' units <- data.frame(size = c(1000, 250, 200, 100, 50),
#'     cost = c(40, 11.5, 8.5, 5.5, 4.4))
#' find_least_cost_partition(450, available_units = units)
#' temp <- find_least_cost_partition(sample(250:450, 10, replace = TRUE), 
#'    available_units = units)
find_least_cost_partition <-
  function(desired_dose,
           available_units,
           subset_col = NULL,
           subset_val = NULL) {
    if(any(is.na(desired_dose)))
      stop("'desired_dose' input cannot have any NA values")
    if(!requireNamespace("lpSolve")){
      stop("must have package 'lpSolve' for finding least cost dose strategies")
    }
    ## input checking
    if (!all(c("size", "cost") %in% names(available_units))){
      stop("available_units must have columns 'size' and 'cost'")
    }
    
    if (is.null(subset_col) + is.null(subset_val) == 1)
      stop("subset_col and subset_val should either both be NULL, or both be set")
    if (!is.null(subset_col) & !is.null(subset_val)) {
      if (length(subset_col) != 1)
        stop("must specify exactly one column to select on")
      if (length(subset_val) != 1)
        stop("must specify exactly one value to select")
      if (!(subset_col %in% names(available_units)))
        stop(paste(subset_col, "is not a column of available_units"))
      if (!(subset_val %in% available_units[, subset_col]))
        stop(paste(
          "column",
          subset_col,
          "does not contain the value",
          subset_val
        ))
      available_units <-
        available_units[available_units[[subset_col]] == subset_val, , 
                        drop = FALSE]
    }
    
    ## set up the problem for lpSolve
    constraint_coefs <- matrix(available_units[, "size"],
                               nrow = 1, byrow = TRUE)
    original_dose <- desired_dose
    
    output_list <- lapply(unique(desired_dose), 
                          function(this_dose){
                            opt_fn(this_dose, available_units, constraint_coefs)
                          }
    )
    
    #   ## combine outputs for all doses
    output_list <- do.call("rbind", output_list)
    names(output_list) <- c("desired_dose", "used_dose", 
                            "waste", "cost") 
    output_list$cost.no.waste <-
      round(with(output_list, cost * desired_dose/used_dose), 2)
    output_list[match(original_dose, output_list$desired_dose),]
  }



opt_fn_ <- function(this_dose, available_units, constraint_coefs){
  lp_soln <-
    lpSolve::lp("min",
                available_units[, "cost"], 
                constraint_coefs,      ## vial sizes
                ">=",                  ## must get at least the dose
                this_dose,
                all.int = TRUE)
  actual_cost <- lp_soln$objval
  used_dose <- lp_soln$solution %*% available_units[, "size"]
  waste <- used_dose - this_dose
  output_list <-
    data.frame(this_dose, used_dose, waste, actual_cost)
  output_list
}

opt_fn <- memoise::memoise(opt_fn_)






#' Cost of administration for an intravenous treatment
#'
#' @param iv_time The infusion time, in units compatible with the costs
#' @param cost_first_unit cost for the first time unit
#' @param cost_addl_units cost for the second time unit
#' @param prorate_first can we charge for less than 1 unit of time?
#'   See details.  
#' @return the cost
#' @details Simple formula:  cost for the first time unit (prorated if 
#'   total time is less than one unit and `prorate_first == TRUE`) 
#'   plus the prorated cost for additional units.
#'   
#'   If `prorate_first == FALSE`, then the smallest time that
#'   can be charged is one time unit.   For example, if time = 0.5
#'   and `prorate_first == TRUE`, the function will charge for one
#'   half of a time unit.   If time = 0.5 and `prorate_first == FALSE`,
#'   the function will charge for a full time unit.
#' @export
#' @examples
#' cost_iv_administration(0.5, 100, 20, prorate_first = TRUE) # = 50
#' cost_iv_administration(1.5, 100, 20, prorate_first = TRUE) # = 110
#' cost_iv_administration(0.5, 100, 20, prorate_first = FALSE) # = 100
#' cost_iv_administration(1.5, 100, 20, prorate_first = FALSE) # = 110
cost_iv_administration <- 
  function(iv_time, cost_first_unit, cost_addl_units, prorate_first = FALSE){
    if(prorate_first)
      first_time <- pmin(iv_time, 1)
    else
      first_time <- 1
    
    first_time * cost_first_unit + 
      pmax(iv_time - 1, 0) * cost_addl_units
    
  }


#' Cost of administration for an intravenous treatment
#'
#' @param data_table a data frame; see details.
#' @param compound the name of the compound.
#' @param time_col,first_cost_col,addl_cost_col parameter names
#'   in data_table.
#' @param prorate_first can we charge for less than 1 unit of time?
#'   See details in [cost_iv_administration()].  
#' @param ... additional arguments, filtering criteria 
#'   that will be passed to [heemod::look_up()].

#' @details `data_table` must have columns `compound`, `param`,
#'   and `value`.   `time_col`, `first_cost_col` and `addl_cost_col`
#'   must be names in the `param` column.  
#'   The required values are found in `data_table`
#'   using [heemod::look_up()].   Then the values are summed using
#'   [cost_iv_administration()].
#'   
#' @return  the cost
#' @export
#'
#' @examples
#' exampleParams <- 
#'   data.frame(compound = "X", 
#'   param = c("cost_admin_first_hr", "cost_admin_addl_hr",
#'             "iv_time_hr"),
#'   value = c(100, 20, 1.5))
#'   cost_iv_compound_administration(exampleParams, "X")
cost_iv_compound_administration <- 
  function(data_table, compound,
           time_col = "iv_time_hr",
           first_cost_col = "cost_admin_first_hr",
           addl_cost_col = "cost_admin_addl_hr",
           prorate_first = FALSE, ...)
    {
    stopifnot(all(c("compound", "param", "value") %in%
                names(data_table)))
    if(!(compound %in% unique(data_table$compound)))
      stop(compound, 
           " not present in 'compound' column of data_table")
    all_cols <- c(time_col, first_cost_col, addl_cost_col)
    cols_present <- all_cols %in% unique(data_table$param)
    if(!all(cols_present))
      stop("parameter",
           plur(sum(!cols_present)),
           " ",
           paste(all_cols[!cols_present], sep = ","),
           " not in data_table"
      )
  iv_time <- 
    heemod::look_up(data_table, compound = compound, 
            param = time_col, value = "value", ...
            )
  cost_first_unit <- 
    heemod::look_up(data_table, compound = compound, 
            param = first_cost_col, value = "value", ...
            )
  cost_addl_units <- 
    heemod::look_up(data_table, compound = compound, 
            param = addl_cost_col, value = "value", ...
            )
  cost_iv_administration(iv_time, 
                         cost_first_unit, 
                         cost_addl_units,
                         prorate_first)
  }


compute_vals_for_adv_ev_ <- function(ae_table){
  required_names <- c("treatment", "ae", "prob")
  if(!all(required_names %in% names(ae_table)))
    stop("adverse events table must have columns ",
         paste(required_names, collapse = ", ")
    )
  val_names <- setdiff(names(ae_table),
                       required_names)
  if(length(val_names) == 0)
    stop("must have some column for a value ", 
         "other than required columns: ",
         paste(required_names, collapse = ", ")
         )
  if(any(missing_prob <- is.na(ae_table$prob))){
    ae_table$prob[is.na(ae_table$prob)] <- 0
    warning("some AE probabilities are missing in the table; ",
            "setting to 0")
  }
  ## if there are missing values for some treatments that
  ##   are defined for other treatments, fill them in.
  ##   if some are actually undefined (that is, not given
  ##   for any treatment), make them 0
  other_names <- setdiff(names(ae_table), required_names)
  values <- dplyr::distinct(ae_table[, c("ae", other_names)])
  values <- values[stats::complete.cases(values), ]
  if(any(is.na(ae_table[, -match(required_names, names(ae_table))]))){
    ae_table <- dplyr::left_join(ae_table[, required_names], 
                                 values, 
                                 by = "ae")
  }
  ## now check whether any AEs have multiple values defined
  values <- values[!duplicated(values),]
  mult_def <- table(values$ae) > 1
  if(any(mult_def))
    stop("multiple values defined for ",
         paste(names(mult_def[mult_def]), collapse = ", ")
    )
  if(any(na_spots <- is.na(ae_table[, other_names]))){
     ae_table[na_spots, other_names] <- 0
     warning("NA values defined for AEs; converting to 0")
   }
  
  ae_table[, other_names] <- 
    ae_table[, other_names, drop = FALSE] * ae_table[, "prob"]
  res <- ae_table[, c("treatment", "prob", other_names)] %>% 
    dplyr::group_by_(~treatment) %>%
      dplyr::summarize_at(.cols = other_names, .funs = "sum") %>%
        dplyr::ungroup()
  res
}
compute_vals_for_adv_ev <- memoise::memoise(compute_vals_for_adv_ev_)

#' Get (computed) values related to adverse events.
#'
#' @param ae_table a table with columns `treatment`, 
#'   `ae` (the names of different adverse events), `prob`
#'   (the probability of the given adverse event occurring under the
#'   given treatment),
#'   and others naming values that can be represented, for example
#'   cost or disutility.
#' @param treatment the treatment for which the value is desired
#' @param value the specific value desired
#'
#' @return the requested value for the requested treatment
#' 
#' @details The underlying function is memoised for efficiency.
#' 
#' @export
#'
#' @examples 
#' AEs <- data.frame(treatment = c("A", "A", "B", "B"),
#'                   ae = c("ae1", "ae2", "ae1", "ae2"),
#'                   prob = c(0.1, 0.1, 0.2, 0),
#'                   cost = c(100, 200, 100, 200))
#' ae_val(AEs, "A", "cost")
#' ae_val(AEs, "B", "cost")
#' 
ae_val <- function(ae_table, treatment, value){
  if(!(value %in% names(ae_table)))
    stop("no column '", value, "' in ae_table")
  this_table <- compute_vals_for_adv_ev(ae_table)
  filter_str <- paste("treatment == '", treatment, "'", sep = "")
  res<- this_table %>% dplyr::filter_(filter_str) %>%
    dplyr::select_(value)
  res <- as.numeric(res[[1]])
  if(length(res) == 0) 
    stop("no AE information returned for treatment ",
         treatment, 
         " and value ",
         value,
         ".\n",
         "Treatments available: ",
         paste(unique(ae_table$treatment), collapse = ", "),
         "\n",
         "Values available: ",
         paste(setdiff(names(ae_table), 
                       c("treatment", "ae", "prob")),
               collapse = ", ")
    )
  res
}

#' Compute weighted drug cost over a distribution of characteristics
#'
#' @param dist distribution for the dose-determining variable
#' @param params arguments for the quantile function for `dist`
#' @param var_base see details
#' @param dose_base see details
#' @param dose_multiplier see details
#' @param available_units available sizes and costs; must have column names
#'   'size' and 'cost'
#' @param subset_col optionally, a subset to select on in available_units
#' @param subset_val if `subset_col` exists, the value to select
#' @param share_vials should vials be shared (or, looked at from the 
#'    other side, should any drug be allowed to go unused?)
#' @param qmin,qmax,by determine the sequence of quantiles of the variable
#'   at which the dose function will be evaluated
#'
#' @details 
#' Typical use for these functions are when dose varies by some
#'   patient characteristic, for example weight or body surface area.
#'   
#' The general formula is 
#'   dose = base_dose + dose_multiplier * pmax(0, var - var_base),
#' and the number of vials is then computed using 
#' [find_least_cost_partition()].
#' 
#' Changing the percentiles at which the doses are calculated by
#' changing `qmin`, `qmax`, and `by` will produce slightly different
#' values, but the values will converge if 
#' 
#' The functions [weighted_norm_dose_costs()] and 
#' [weighted_lognorm_dose_costs()] are provided to cover two common cases.
#' @return the weighted dose cost
#' @export
#'
#' @examples
#' vialCost <- data.frame(treatment = "fake", 
#'                       size = c(50, 250), 
#'                       cost = c(2521, 12297)
#'                       )
#' weighted_dose_costs("norm", params = c(mean = 1.85, sd = 0.25),
#'                     var_base = 0, 
#'                     dose_base = 0, 
#'                     dose_multiplier = 320, 
#'                     vialCost, "treatment", "fake",
#'                     share_vials = FALSE)
weighted_dose_costs <- function(dist, params, var_base, dose_base,
                                dose_multiplier, available_units, subset_col, subset_val,
                                share_vials,
                                qmin = 0.01, qmax = 0.99, by = 0.01){

  ## the function originally assumed a single value
  ##   of each parameter, but can be called with one 
  ##   value per cycle when evaluating heemod parameters,
  ##   which would throw off the averaging;
  ##   add this code to get the unique value or throw an error
  for(param_name in names(params)){
    unique_val <- unique(params[[param_name]])
    if(length(unique_val) > 1)
      stop("can only have one value for parameter '",
           param_name,
           "'", 
           sep = "")
    params[[param_name]] <- unique_val
  }
  
  qfn <- paste("q", dist, sep = "")
  if(qmin <= 0)
    stop("it does not make sense for 'qmin' to be <= 0")
  if(qmax >= 1)
    stop("it does not make sense for 'qmax' to be >= 1")
  prob_pts <- seq(from = qmin, to = qmax, by = by)
  formal_args <- setdiff(names(formals(qfn)), c("p", "lower.tail", "log.p"))
  missing_args <- setdiff(formal_args, names(params))
  extra_args <- setdiff(names(params), formal_args)
  if(length(missing_args) | length(extra_args))
    stop("mismatched arguments for function '",
         qfn,
         "':\n",
         ifelse(length(missing_args), 
                paste("missing args: ", 
                      paste(missing_args, collapse = ", "),
                      ";\n"),
                ""),
         ifelse(length(extra_args),
                paste("unused arguments: ",
                      paste(extra_args, collapse = ", "),
                      ";\n"),
                "")
         )
  arg_list <- c(params, p = list(prob_pts))
  quantiles <- do.call(qfn, arg_list)
  dose <- dose_base + pmax(0, (quantiles - var_base) * dose_multiplier)
  costs <- find_least_cost_partition(dose, available_units, subset_col = subset_col, subset_val = subset_val)
  if(share_vials) 
    costs <- costs$cost.no.waste
  else 
    costs <- costs$cost 
  sum(costs)/length(prob_pts)
}

#' @rdname weighted_dose_costs
#' @param mean mean of normal distribution
#' @param sd sd of normal distribution
#' @export
#' @examples
#' vialCost <- data.frame(treatment = "fake", 
#'                       size = c(50, 250), 
#'                       cost = c(2521, 12297)
#'                       )
#' weighted_norm_dose_costs(1.85, 0.25, 
#'                          var_base = 0, dose_base = 0, 
#'                          dose_multiplier = 320, vialCost, "treatment", "fake",
#'                          share_vials = FALSE)
#' 
weighted_norm_dose_costs <- function(mean, sd, var_base, dose_base, 
                                     dose_multiplier, available_units, subset_col, subset_val,
                                     share_vials,
                                     qmin = 0.01, qmax = 0.99, by = 0.01){
  
weighted_dose_costs("norm", params = list(mean = mean, sd = sd),
                      var_base, dose_base, dose_multiplier, 
                      available_units, subset_col, subset_val,
                      share_vials = share_vials,
                      qmin, qmax, by)
}

#' @rdname weighted_dose_costs
#' @param meanlog mean of lognormal distribution
#' @param sdlog sd of lognormal distribution
#' @export
#' @examples
#' vialCost <- data.frame(treatment = "fake", size = 50, cost = 1000)
#' weighted_lognorm_dose_costs(4.3423559, 0.2116480, 
#'                             var_base = 50, dose_base = 100, 
#'                             dose_multiplier = 2, vialCost, "treatment", "fake",
#'                             share_vials = FALSE)
#' weighted_lognorm_dose_costs(4.3423559, 0.2116480, 
#'                             var_base = 50, dose_base = 100, 
#'                             dose_multiplier = 2, vialCost, "treatment", "fake",
#'                             share_vials = TRUE)

weighted_lognorm_dose_costs <- function(meanlog, sdlog, var_base, dose_base, 
                                        dose_multiplier, available_units, subset_col, subset_val,
                                        share_vials,
                                        qmin = 0.01, qmax = 0.99, by = 0.01){
  weighted_dose_costs("lnorm", params = list(meanlog = meanlog, sdlog = sdlog),
                      var_base, dose_base, dose_multiplier, 
                      available_units, subset_col, subset_val,
                      share_vials = share_vials, qmin, qmax, by)
}

