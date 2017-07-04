allowed_fit_distributions <- c("exp", "weibull", "lnorm", "llogis", 
  "gamma", "gompertz", "gengamma")

#' Create a partitioned survival model object from tabular input
#'
#' @param base_dir directory where the reference file is located
#' @param ref_file Name of the reference file.
#' @param ref data frame read from the reference file.
#' @param state_names names of states in the strategies
#' @param save_fits should the models be saved to disk?
#' @param just_load should we just load the models instead of fitting?
#' @param ... additional arguments to be passed to 
#'   [survival_fits_from_ref_struc()]
#' @export
#' @return a tibble of partitioned survival objects 
#'   (for `part_surv_from_tabular`)
#'   or a list (for `survival_fits_from_tabular`), 
#'   with each element being a similar
#'   tibble of flexsurvreg fits for either overall survival 
#'   or progression-free survival.
#'   There are columns to distinguish the different distributions used
#'   used in survival fitting, different treatments, and different subsets.
#' @examples 
#' fit_tibble <- partitioned_survival_from_tabular(system.file("tabular/surv",
#'                                                             package = "heemod"), 
#'                                   "example_oncSpecs.csv", 
#'                                   c("ProgressionFree", "Progressive", 
#'                                   "Terminal", "Death"), 
#'                                   save_fits = FALSE,
#'                                   just_load = FALSE)
survival_fits_from_tabular <- function(base_dir, ref_file, ...) {
  ref <- read_file(file.path(base_dir, ref_file))
  ref$full_file <- file.path(base_dir, ref$file)
  survival_fits_from_ref_struc(ref, ...)
}

#' @rdname survival_fits_from_tabular
#' @export
survival_fits_from_ref_struc <- function(ref, 
                                         save_fits = FALSE, just_load = FALSE, 
                                         ...){  
  surv_ref_full_file <- ref[ref$data == "tm", "full_file"]
  surv_ref_file <- ref[ref$data == "tm", "file"]
  ## slightly roundabout way of getting the base location back
  location <- gsub(paste0(surv_ref_file, "$"),
                   "",
                   surv_ref_full_file)
  
  survival_specs <- read_file(file.path(surv_ref_full_file))
  ##dists = c("exp", "weibull", "lnorm", "llogis", 
  ##          "gamma", "gompertz", "gengamma")
  survival_from_data(location,
                     survival_specs,
                     dists = allowed_fit_distributions, 
                     save_fits = save_fits,
                     just_load = just_load,
                     ...)
  
  
}

#' @export
#' @rdname survival_fits_from_tabular
partitioned_survival_from_tabular <- 
  function(base_dir, ref_file,
           state_names,
           save_fits = FALSE, just_load = FALSE) {
  
  ref <- read_file(file.path(base_dir, ref_file))
  
  surv_inputs <- survival_fits_from_tabular(base_dir, ref_file, 
                             save_fits, just_load)
  heemod::part_survs_from_surv_inputs(surv_inputs, state_names)
}

#' @export
#' @rdname survival_fits_from_tabular
partitioned_survival_from_ref_struc <- function(ref, 
                                              state_names,
                                              save_fits = FALSE, 
                                              just_load = FALSE) {
  surv_inputs <- survival_fits_from_ref_struc(ref, 
                                            save_fits, just_load)
  heemod::part_survs_from_surv_inputs(surv_inputs[[1]], state_names)
}

#' Title Get survival analysis curves from data
#'
#' @param location base directory for the analysis.
#' @param survival_specs a data.frame containing information
#'   on how to perform each fit - see @details.
#' @param dists the distributions to use to fit the survival function
#' @param save_fits should fits be saved to disk?  Can be useful for testing.
#' @param just_load If TRUE, data files are ignored in favor of loading fits
#'    from `fit_files`
#' @param report_file_name base name for survival output - Excel and the
#'    diagnostic tables and plots on the fits.
#' @param set_definitions definitions of different subsets to fit
#' @param write_report should the call produce an html report with
#'    diagnostic tables and plots on the fits?
#' @param write_excel should the call produce a csv file summarizing fits?
#' @param report_template If `write_report == TRUE`, 
#'   a template for the report.
#' @return a tibble with columns `type`, `treatment`, `set_name`,
#    `dist`, `fit`, `set_def`, and `time_subtract`
#' @details By default, the function fits with seven 
#'   different distribution functions:
#'   exponential,  Weibull,  lognormal, log-logistic, 
#'   Gompertz, gamma, and generalized gamma.
#'   
#'   `write_report` and `write_excel` are ignored if `just_load = TRUE`.
#' 
#' survival_specs contains information about how to create the fits:
#'   the directory (or directories) in which
#'   files are located, the file names, the names of the columns
#'   that should be used for time, events, and censorship status,
#'   and where and under what names to save the fits.
#' If data_files is NULL and fit files exists, then fit_files
#'   should have the names of files in which survival models are kept.
#'   If data_files is not NULL, then survival models will be fit from
#'   the data in data_files (using the `flexsurvreg` package), and if fit_files
#'   is also not NULL, the survival models will be saved in fit_files.
#'   
#' If `data_directory` begins with "\\", with the platform file separator
#' (`.Platform$file.sep`) or with a letter and a colon (to accomodate Windows),
#' it's considered absolute and not appended to `location`.

survival_from_data <-
  function(location,
           survival_specs,
           dists = dists,
           save_fits = TRUE,
           just_load = FALSE,
           set_definitions = "set_definitions.csv",
           report_file_name = "survival_report",
           write_report = FALSE,
           write_excel = write_report,
           report_template) {
    survival_specs <- check_survival_specs(survival_specs)
    
    if (just_load) {
      return(heemod::load_surv_models(location, survival_specs))
    }
    else{
      ## going to check whether we have an absolute directory
      ##   for data directory, or relative.
      loadNamespace("flexsurv")
      is_absolute <-
        grepl("^[A-Z]:", survival_specs$data_directory) |
        grepl("^(/|\\\\)", survival_specs$data_directory)
      data_files <- file.path(location,
                              survival_specs$data_directory,
                              survival_specs$data_file)
      data_files[is_absolute] <-
        file.path(survival_specs$data_directory,
                  survival_specs$data_file)[is_absolute]
      fit_files <- file.path(location,
                             survival_specs$fit_directory,
                             survival_specs$fit_file)
      
      
      surv_models <-
        lapply(seq(1:nrow(survival_specs)),
               function(this_row) {
                 ## read in the data and filter on whether patient treated
                 ##   and, if so, which treatment
                 filter_strings <- make_filter_strings(survival_specs, this_row)
                 this_data <-
                   read_file(data_files[this_row]) %>%
                    dplyr::filter_(filter_strings$treatment_flag_str) %>%
                      dplyr::filter_(filter_strings$choose_treatment_str)
                 
                 ## set up the event values: 1 for event, 0 for censored
                  this_data <-
                    fix_censor_col(this_data,
                                  survival_specs[this_row, "censor_col"],
                                  survival_specs[this_row, "censor_code"],
                                  survival_specs[this_row, "event_code"],
                                  data_files[this_row])
                 
                 ## get set definitions, if there are any
                 ## (if not, will return a data frame with no rows)
                 set_definitions <-
                   get_set_definitions(file.path(location,
                                                 survival_specs$fit_directory[this_row]),
                                       file_name = set_definitions)
                 
                 ## get the set definitions associated with this treatment
                 these_sets <-
                   dplyr::filter_(set_definitions,
                                  ~ treatment == survival_specs[this_row, "treatment"])
                 ## if relevant, restrict to the set definitions associated
                 ##   with this type (pfs or os)
                 if ("type" %in% names(these_sets)) {
                   these_sets <-
                     dplyr::filter_(these_sets,
                                    ~ type == survival_specs[this_row, "type"])
                 }
                 
                 if (nrow(these_sets) == 0)
                   these_sets <- data.frame(
                     set_name = "all",
                     condition = "TRUE",
                     time_subtract = 0,
                     stringsAsFactors = FALSE
                   )
                 surv_fits <-
                   lapply(1:nrow(these_sets), function(set_index) {
                     subset_data <-
                       this_data %>%
                       dplyr::filter_(.dots = strsplit(these_sets[set_index,
                                                           "condition"],
                                                ",")[[1]])
                     
                     subset_data <-
                       subtract_times(
                         subset_data,
                         to_subtract = these_sets[[set_index, "time_subtract"]],
                         this_time_col = survival_specs[this_row, "time_col"],
                         set_name = these_sets[[set_index, "set_name"]],
                         condition = these_sets[[set_index, "condition"]]
                       )
                     
                     these_surv_fits <-
                       f_fit_survival_models(
                         subset_data,
                         dists = dists,
                         time_col_name = survival_specs[this_row, "time_col"],
                         censor_col_name = survival_specs[this_row, "censor_col"],
                         treatment_col_name = survival_specs[this_row, "treatment_col"],
                         covariate_col_names = NULL,
                         fit_indiv_groups = TRUE
                       )
                     these_surv_fits$set_name <-
                       these_sets[[set_index, "set_name"]]
                     these_surv_fits$set_def <-
                       these_sets[[set_index, "condition"]]
                     these_surv_fits$time_subtract <-
                       these_sets[[set_index, "time_subtract"]]
                     
                     ## modify the fits to take into account
                     ##   the time subtraction
                     new_fits <-
                       these_surv_fits %>%
                       dplyr::rowwise() %>%
                       dplyr::do_(fit = ~ heemod::apply_shift(dist = .$fit,
                                                      shift = .$time_subtract)) %>%
                       dplyr::ungroup()
                     
                     
                     these_surv_fits$fit <- new_fits$fit
                     
                     these_surv_fits
                   })
                 names(surv_fits) <- these_sets$set_name
                 surv_fits_tib <- do.call("rbind", surv_fits)
                 surv_fits_tib$type <-
                   survival_specs[this_row, "type"]
                 surv_fits_tib <-
                   surv_fits_tib[, c("type",
                                     "treatment",
                                     "set_name",
                                     "dist",
                                     "fit",
                                     "set_def",
                                     "time_subtract"),
                                 drop = FALSE]
                 assign(survival_specs[this_row, "fit_name"],
                        surv_fits_tib)
                 if (save_fits) {
                   save(list = survival_specs[this_row, "fit_name"],
                        file = paste(fit_files[this_row], ".RData", sep = ""))
                 }
                 
                 surv_fits_tib
               })
      
      names(surv_models) <- survival_specs$fit_name
      surv_models <- do.call("rbind", surv_models)

      
      if (write_report) {
        if(!file.exists(report_template))
          stop("cannot find report template file ",
               report_template)
        rmarkdown::render(
          report_template,
          params = list(fit_tibble = surv_models,
                        extra_args = NULL),
          output_file = file.path(location,
                                  paste(report_file_name, ".html",
                                        sep = ""))
        )
      }
      ## we don't put out graphs for multi-treatment models,
      ##  but we do save them and put the parameters into Excel
      surv_models <- add_cross_treatment_fits(surv_models, dists)
      
      if (write_excel) {
        write_fits_to_excel_from_tibble(
          fit_tibble = surv_models,
          wb = file.path(location,
                         paste(report_file_name, ".xlsx", sep = "")),
          alignment = "horizontal"
        )
      }
      
    }
    return(surv_models)
  }


get_set_definitions <-
  function(data_dir, file_name = "set_definitions") {
    set_definitions <- data.frame(
      treatment = character(0),
      set_name = character(0),
      condition = character(0)
    )
    set_definition_file_name <-
      list.files(data_dir, pattern = file_name, full.names = TRUE)
    if (length(set_definition_file_name) > 1)
      stop("multiple files matching set_definition file name")
    if (length(set_definition_file_name) == 1)
      set_definitions <- read_file(set_definition_file_name)
    missing_names <-
      setdiff(c("treatment", "set_name", "condition"),
              names(set_definitions))
    extra_names <- setdiff(names(set_definitions),
                           c("treatment", "set_name",
                             "condition", "time_subtract"))
    extra_names <- extra_names[!grep("^.comment", extra_names)]
    if (length(missing_names) > 0)
      stop("set_definitions file missing column(s): ",
           paste(missing_names, collapse = ", "))
    if (length(extra_names) > 0)
      stop("unrecognized column name(s): ",
           paste(extra_names, collapse = ", "))
    ## just in case we have only logicals
    set_definitions$condition <-
      as.character(set_definitions$condition)
    
    if (!("time_subtract" %in% names(set_definitions))) {
      set_definitions <-
        dplyr::mutate_(set_definitions, time_subtract = 0)
    }
    else{
      neg_offset <- which(set_definitions$time_subtract < 0)
      if (length(neg_offset))
        stop("bad offset in set_definitions line(s) ",
             paste(neg_offset, collapse = ", "))
      
      ## check whether time_subtracts match numbers in conditions,
      ##  warning if not
      time_matches_condition <-
        sapply(1:nrow(set_definitions),
               function(i) {
                 if (is.na(set_definitions$time_subtract[i]))
                   TRUE
                 else
                   grepl(set_definitions$time_subtract[i],
                         set_definitions$condition[i])
               })
      problem_conditions <- which(!time_matches_condition)
      if (length(problem_conditions)) {
        bad_conds <-
          sapply(problem_conditions, function(this_row) {
            paste(
              this_row,
              set_definitions$condition[this_row],
              set_definitions$time_subtract[this_row],
              sep = "\t"
            )
          })
        bad_conds <- paste(bad_conds, collapse = "\n")
      
      
      message <-
        paste(
          "\n",
          bad_conds,
          "\n",
          "value of 'time_subtract' does not appear in 'condition' for row",
          plur(length(problem_conditions)),
          " ",
          paste(problem_conditions, collapse = ", "),
          " in set_definitions",
          "\n",
          sep = ""
        )
      warning(message)
      }
    }
    
    if ("type" %in% names(set_definitions))
      set_definitions$type <- toupper(set_definitions$type)
    else{
      num_rows <- nrow(set_definitions)
      set_definitions <- rbind(set_definitions, set_definitions)
      set_definitions$type <- rep(c("PFS", "OS"), each = num_rows)
      type_col <- match("type", names(set_definitions))
      set_definitions <-
        set_definitions[c(type_col, setdiff(1:ncol(set_definitions), type_col))]
    }
    set_definitions <- set_definitions[!duplicated(set_definitions),]
    dups <-
      duplicated(set_definitions[, c("treatment", "type", "set_name")])
    if (any(dups)) {
      stop("it appears you have multiple definitions of some sets: ",
           set_definitions[dups, c("treatment", "type", "set_name", "condition")])
    }
    
    set_definitions <-
      set_definitions %>%
      dplyr::mutate_(time_subtract = ~ ifelse(is.na(time_subtract),
                                              0, time_subtract))
    set_definitions
  }

#'
#' Fit survival models using different distributions
#' @export
#' @param survdata Survival data to be used.
#' @param dists Distributional forms to be considered in fitting using
#'     `flexsurvreg`.  By default, includes exponential, Weibull,
#'     lognormal, log-logistic, gamma, gompertz, and generalized gamma.  Kaplan-Meier
#'     curves will also be stored automatically, with distribution
#'     name "km".
#' @param treatment_col_name Name of the column in survdata that holds the
#'     treatment group names (for example "control", "treatment",
#'     "dose1", and so on).
#' @param time_col_name Name of the column in survdata with event times.
#' @param censor_col_name Name of the column in survdata with censorship
#'     indicators.
#' @param covariate_col_names  Not yet implemented
#' @param fit_indiv_groups Should groups be fit individually?  Default TRUE.
#' @return a matrix (with dimnames) of flexsurvreg objects,
#'     with rows corresponding to
#'     the distributions, and columns corresponding to groups.
#' @details
#'  `survdata` should be a data frame of the form required for
#'  [flexsurv::flexsurvreg()].   In the column specified
#'  by `censor_col_name` in the data itself,
#'  1 means that an event (frequently death or disease progression,
#'  depending on context) was observed, while 0
#'  means the event was censored.  So the Kaplan-Meier plot will have
#'  a drop anywhere the censor column contains a 1, and will not
#'  contain a drop when the censor column contains a 0.
#'
f_fit_survival_models <-
  function(survdata,
           dists = allowed_fit_distributions,
           time_col_name,
           censor_col_name,
           treatment_col_name,
           covariate_col_names = NULL,
           fit_indiv_groups = TRUE) {
    ##
    ## TODO:  so far, covariate_col_names is only a placeholder
    ##  what needs to be decided is whether to allow all combinations
    ##  of covariates, or have people enter the combinations they want,
    ##  or somehow allow both ways (which might be better handled by
    ##  different functions entering in)
    if (!requireNamespace("flexsurv", quietly = TRUE))
      stop("flexsurv package needed to fit survival models")
    if (!requireNamespace("survminer", quietly = TRUE))
      stop("survminer package needed for survival models")
    stopifnot(is.data.frame(survdata))
    stopifnot(nrow(survdata) > 0)
    stopifnot(time_col_name %in% names(survdata))
    stopifnot(censor_col_name %in% names(survdata))
    stopifnot(treatment_col_name %in% names(survdata))
    stopifnot(is.null(covariate_col_names) |
                all(covariate_col_names %in% names(survdata)))
    
    ## make a list of the subgroups, and add a group of all together
    unique_groups <-
      as.character(unique(survdata[, treatment_col_name]))
    if ("all" %in% unique_groups)
      stop("can't name a treatment 'all'")
    
    if (length(unique_groups) > 1) {
      groups_list <- c(list(unique_groups))
      names(groups_list) <- c("all")
    }
    else{
      groups_list <- as.list(unique_groups)
      names(groups_list) <- unique_groups
    }
    
    formula_base_string <-
      paste("survival::Surv(",
            time_col_name,
            ", ",
            censor_col_name,
            ")",
            " ~",
            sep = "")
    
    ## we're going to add "km" to the fit distributions
    ##   and put a survfit object there
    
    dists <- c(dists, "km")
    
    ## cycle through combinations of distributions and subsets,
    ##   getting survival analysis results at each step
    #    conditions <- expand.grid(dist = dists, group = names(groups_list),
    #                              stringsAsFactors = FALSE)
    conditions <- data.frame(dist = dists,
                             stringsAsFactors = FALSE)
    all_res <-
      lapply(1:nrow(conditions), function(i) {
        ##this_group_set <- groups_list[[ conditions[i, "group"] ]]
        this_data <-
          ##  survdata[survdata[[treatment_col_name]] %in% this_group_set,
          ## c(time_col_name, censor_col_name, treatment_col_name)]
          survdata[, c(time_col_name, censor_col_name, treatment_col_name)]
        
        ## only include the treatment group in the formula
        ##   if there is more than one treatment in the data set
        num_treatments <-
          length(unique(this_data[, treatment_col_name]))
        if (num_treatments > 1)
          this_formula <-
          paste(formula_base_string, treatment_col_name)
        else
          this_formula <- paste(formula_base_string, "1")
        this_formula <- stats::as.formula(this_formula)
        if (conditions[i, "dist"] == "km") {
          this_fit <- survminer::surv_fit(this_formula,
                                          data = this_data)
        }
        else{
          this_fit <- try(flexsurv::flexsurvreg(this_formula,
                                                data = this_data,
                                                dist = conditions[i, "dist"]),
                          silent = TRUE)
        }
      })
    all_res <-
      f_add_surv_fit_metrics(all_res, metrics = c("BIC", "m2LL"))
    
    tibble::tibble(
      dist = rep(dists, length(groups_list)),
      treatment = rep(names(groups_list), each = length(dists)),
      fit = all_res
    )
    
  }


#' Choose the best model out of a set based on a metric.
#'
#' @param surv_fits A list object from `f.get.survival.probs` that gives
#' 		a collection (list) of parametric survival fit object.
#' @param metric The metric to choose the model by.
#' 		Currently supports selecting the best model using one (and only one) of the
#' 		following metrics:
#' 			1. "AIC": Akaike information criterion
#' 			2. "BIC": Bayesian information criterion
#' 			3. "m2LL": -2*log likelihood
#' @details
#' Gets the best survival model, according to a particular metric.
#' 	Current implementation is limited to selecting the best model
#' 	 using one (and only one) of the following metrics:
#' 		1. Akaike information criterion (AIC)
#' 		2. Bayesian information criterion (BIC)
#' 		3. -2*log likelihood (-2LL)
#' @return The best model according to the chosen metric.

f_get_best_surv_model <-
  function(surv_fits,
           metric = c("AIC")) {
    #argument checks
    stopifnot(length(surv_fits) > 0)
    stopifnot(length(metric) == 1)
    stopifnot(metric %in% c("AIC", "BIC", "m2LL"))
    
    #order surv_fits by metric priority
    sorted = surv_fits[order(sapply(surv_fits, '[[', metric))]
    
    #return best model
    best_model <- sorted[[1]]
    attr(best_model, "dist") <- names(sorted)[1]
    best_model
  }

#'
#' Calculate additional metrics to evaluate fit of survival model.
#'
#' @param surv_fits A list object from [f_fit_survival_models()] that gives
#' 		a collection (list) of `flexsurvreg` parametric survival fit object.
#' @param metrics Metrics to calculate.
#' @return
#'   A list object of parametric survival fits, containing additional fields for
#' 		the calculated fit metrics.
#' @details Currently calculates only:
#' 		\itemize{\item Bayesian information criterion (BIC)
#' 		\item -2*log likelihood (-2LL)}  (Objects come with AIC already calculated.)

f_add_surv_fit_metrics <-
  function(surv_fits, metrics = c("BIC", "m2LL")) {
    #argument checks
    stopifnot(length(surv_fits) > 0)
    stopifnot(all(metrics %in% c("BIC", "m2LL")))
    
    #get current and previous time step in absolute (not Markov) units
    for (metric in metrics)
    {
      if (metric == "BIC") {
        surv_fits = add_BIC(surv_fits)
      }
      if (metric == "m2LL") {
        surv_fits = add_m2LL(surv_fits)
      }
    }
    
    #returns updated survival fits object
    surv_fits
  }


add_BIC <- function(surv_fits)
{
  out = surv_fits
  for (i in 1:length(out))
  {
    if (inherits(out[[i]], "flexsurvreg")) {
      this_BIC = (-2 * getElement(out[[i]], "loglik") +
                    (log(as.numeric(
                      getElement(out[[i]], "N")
                    )) * getElement(out[[i]], "npars")))
      out[[i]]$BIC <- this_BIC
    }
  }
  out
}

add_m2LL <- function(surv_fits)
{
  out = surv_fits
  for (i in 1:length(out))
  {
    if (inherits(out[[i]], "flexsurvreg")) {
      this_m2LL = -2 * getElement(out[[i]], "loglik")
      out[[i]]$m2LL <- this_m2LL
    }
  }
  out
}

#' Add fit metrics to a fit tibble
#'
#' @param fit_tib a tibble of fits, for example from
#'   [survival_fits_from_tabular()]
#' @param metric which metrics to add
#'
#' @return a tibble with the added columns
#' @export
#'
#' @examples
extract_surv_fit_metrics <-
  function(fit_tib, metric = c("AIC", "BIC", "m2LL")) {
    ## metric <- match.arg(metric)
    fit_tib <-
      fit_tib %>% dplyr::filter_( ~ dist != "km")
    fit_tib$fit <-
      lapply(fit_tib$fit, extract_fits)
    extracted <-
      fit_tib %>%
      dplyr::rowwise() %>%
      dplyr::do_( ~ data.frame(.$fit[metric]))
    tibble::as_tibble(cbind.data.frame(fit_tib, extracted))
  }



#' Recode censoring / event column
#'
#' @param this_data data set
#' @param this_censor_col the column to change
#' @param censor_code code for censoring
#' @param event_code code for events
#' @param file_name the name of the file this came from,
#'   in case it's needed for an error message
#'
#' @return the same data set, with the column recoded
#'
fix_censor_col <- function(this_data,
                           this_censor_col,
                           censor_code,
                           event_code,
                           file_name) {
  ## set up the event values: 1 for event, 0 for censored
  ## this_censor_col <- survival_specs[this_row, "censor_col"]
  if (!this_censor_col %in% names(this_data))
    stop(
      "censoring status column '",
      this_censor_col,
      "' does not exist in data file '",
      file_name,
      "'"
    )
  this_data[, this_censor_col] <-
    match(this_data[, this_censor_col],
          c(censor_code, event_code)) - 1
  if (any(is.na(this_data[, this_censor_col])))
    stop(
      "non-matching values in ",
      this_censor_col,
      "; all values should be either '",
      event_code,
      "' (for events) or '",
      censor_code,
      "' (for censoring)"
    )
  this_data
}


#' Subtract time
#'
#' @param subset_data some data
#' @param to_subtract the number to subtract from times
#' @param this_time_col the name of the time column
#' @param set_name in case needed for error message
#' @param condition in case needed for error message
#' @details This is for fitting subsets that are defined by
#'    being after a certain time, and we want to
#'   fit as if the cutoff time is the new base time.
#'
#' Moved into a separate function to get the error check
#'   and lengthy construction of an error messageout of the main flow.
#' @return the data, with `to_subtract`` subtracted from the time column
#'
subtract_times <- function(subset_data,
                           to_subtract,
                           this_time_col,
                           set_name,
                           condition) {
  ## subtract off the appropriate time, if there is one
  if (length(to_subtract) > 0 && !is.na(to_subtract)) {
    subset_data[, this_time_col] <-
      subset_data[, this_time_col] - to_subtract
  }
  ## error check:
  ##   doesn't make sense to have times = 0 or negative;
  ##   most parametric fits will fail
  
  negatives <- subset_data[, this_time_col] < 0
  zeros <- subset_data[, this_time_col] == 0
  
  
  if (any(negatives | zeros)) {
    num_zeros <- sum(zeros)
    num_negatives <- sum(negatives)
    piece1 <- piece2 <- ""
    if (num_zeros > 0) {
      piece1 <- paste(num_zeros,
                      " '0' value",
                      plur(num_zeros),
                      sep = "")
    }
    if (num_negatives > 0) {
      piece2 <- paste(num_negatives,
                      " negative value",
                      plur(num_negatives),
                      sep = "")
    }
    if (nchar(piece1) & nchar(piece2))
      message <- paste(piece1, piece2, sep = " and ")
    else
      message <- ifelse(nchar(piece1), piece1, piece2)
    message <- paste(message,
                     " in the subset ",
                     # these_sets[set_index, "set_name"],
                     set_name,
                     "'",
                     "with condition '",
                     condition,
                     "'",
                     sep = "")
    if (grepl(">=", condition))
      #these_sets[set_index, "condition"]))
      message <- paste(message,
                       "\n",
                       "perhaps you need '>' in your condition rather than '>='?")
    
    stop(message)
  }
  subset_data
}


make_filter_strings <- function(survival_specs, this_row) {
  this_treatment <-
    survival_specs[this_row, "treatment"]
  treatment_col_name <-
    survival_specs[this_row, "treatment_col"]
  choose_treatment_str <-
    paste(treatment_col_name, " == '", this_treatment, "'", sep = "")
  
  ## filter by the treatment flag column
  ## (typically intent to treat or actual treatment)
  treatment_flag_name <-
    survival_specs[this_row, "treatment_flag_col"]
  treatment_flag_str <-
    paste(
      treatment_flag_name,
      "==",
      "'",
      c("Y", "Yes", "1", "TRUE", "T"),
      "'",
      sep = "",
      collapse = " | "
    )
  list(treatment_flag_str = treatment_flag_str,
       choose_treatment_str = choose_treatment_str)
}


add_cross_treatment_fits <- function(fit_tibble, dists) {
  ##fit_groups <- fit_tibble[, c("type", "set_name", "set_def", "time_subtract")]
  
  ## the point of all this is to take only those sets that have more
  ##   than one treatment (if all subsets are defined for all treatments,
  ##   you could just use a !duplicated)
  fit_groups <-
    fit_tibble %>%
    dplyr::select(-dplyr::matches("fit")) %>%
    dplyr::select(-dplyr::matches("dist")) %>%
    dplyr::distinct() %>%
    dplyr::group_by_( ~ type, ~ set_name, ~ set_def, ~ time_subtract) %>%
    dplyr::summarize_(num = ~ n()) %>%
    dplyr::filter_( ~ num > 1) %>%
    dplyr::select(-dplyr::matches("num"))

  ## change how this is calculated from dplyr to
  ## base R so that we don't set of "no binding
  ##   for global variables" in R CMD CHECK
  # new_fits <-
  #   fit_groups %>% dplyr::rowwise() %>%
  #   dplyr::do(fit = do.call(
  #     get_cross_fit,
  #     list(
  #       fit_tibble = fit_tibble,
  #       set_name = .$set_name,
  #       type = .$type,
  #       dists = dists,
  #       time_subtract = .$time_subtract,
  #       set_def = .$set_def
  #     )
  #   ))
  
  new_fits <- lapply(1:nrow(fit_groups), function(i){
    do.call(get_cross_fit,
            list(fit_tibble = fit_tibble,
                 set_name = fit_groups[[i, "set_name"]],
                 type = fit_groups[[i, "type"]],
                 dists = dists,
                 time_subtract = fit_groups[[i, "time_subtract"]],
                 set_def = fit_groups[[i, "set_def"]])
            )
  })
  ##  dplyr::bind_rows(fit_tibble, tidyr::unnest(new_fits))
  dplyr::bind_rows(fit_tibble, dplyr::bind_rows(new_fits))
  }

get_cross_fit <- function(fit_tibble,
                          type,
                          set_name,
                          dists,
                          time_subtract,
                          set_def) {
  partial1 <-
    dplyr::filter_(
      fit_tibble,
      ~ dist == "km",
      lazyeval::interp( ~ type == var, var = type),
      lazyeval::interp( ~ set_name %in% var, var = set_name)
    )
  reassembled_data <-
    dplyr::bind_rows(lapply(partial1$fit,
                            function(x) {
                              y <- extract_data(x)
                              names(y) <-
                                c("time", "status", "treatment")
                              y
                            }))
  
  cross_treatment_fits <-
    f_fit_survival_models(
      reassembled_data,
      dists = dists,
      time_col_name = "time",
      censor_col_name = "status",
      treatment_col_name = "treatment",
      covariate_col_names = NULL,
      fit_indiv_groups = FALSE
    )
  
  cross_treatment_fits$type <- type
  cross_treatment_fits$set_name <- set_name
  cross_treatment_fits$time_subtract <- time_subtract
  cross_treatment_fits$set_def <- set_def
  cross_treatment_fits
}


## make sure that survival fit specifications are 
##   properly formatted and have reasonable data
check_survival_specs <- 
  function(surv_specs){
    ## if there are troubles with absolute file paths, 
    ##   might want to add an "absolute" column to surv_specs
    ##   to be explicit (if, possibly, a little redundant)
    if(!is.data.frame(surv_specs) || nrow(surv_specs) == 0)
      stop("surv_specs must be a data frame with at least one row")
    if(!("event_code" %in% names(surv_specs))){
      warning("event_code not defined in surv_specs; setting to 1 for all rows")
      surv_specs$event_code <- 1
    }
    if(!("censor_code" %in% names(surv_specs))){
      warning("censor_code not defined in surv_specs; setting to 0 for all rows")
      surv_specs$censor_code <- 0
    }
    if(!"treatment_flag_col" %in% names(surv_specs)){
      warning("treatment_flag_col not defined in surv_specs;",
              "setting to 'ITTFL' for PFS and OS,  'TRTFL' for ToT")
      surv_specs$treatment_flag_col <- "ITTFL"
      surv_specs[surv_specs$type == "ToT", "treatment_flag_col"] <- "TRTFL"
    }
    if(any(is.na(surv_specs$treatment_flag_col))){
      warning("treatment_flag_col is 'NA' for some rows;",
              "setting to 'ITTFL' for PFS and OS,  'TRTFL' for ToT")
      surv_specs[is.na(surv_specs$treatment_flag_col) & 
                   surv_specs$type %in% c("PFS", "OS", "pfs", "os")] <- 
        "ITTFL"
      surv_specs[is.na(surv_specs$treatment_flag_col) & 
                   surv_specs$type %in% c("ToT")] <- 
        "TRTFL"
    }
    
    surv_spec_col_names <- c("type", "treatment", "data_directory", "data_file",
                             "fit_directory", "fit_name", "fit_file",
                             "time_col", "treatment_col", "censor_col",
                             "event_code", "censor_code", "treatment_flag_col"
    )
    if(! identical(sort(names(surv_specs)), sort(surv_spec_col_names))){
      extra_names <- setdiff(names(surv_specs), surv_spec_col_names)
      missing_names <- setdiff(surv_spec_col_names, names(surv_specs))
      names_message <- paste("surv_ref must have column names:\n",
                             paste(surv_spec_col_names, collapse = ", "), 
                             sep = "")
      if(length(extra_names) > 0)
        names_message <- paste(names_message, "\n", "extra names: ",
                               paste(extra_names, collapse = ", "),
                               sep = "")
      if(length(missing_names) > 0)
        names_message <- paste(names_message, "\n", "missing names: ",
                               paste(missing_names, collapse = ", "),
                               sep = "")
      stop(names_message)
    }
    
    if(any(is.na(surv_specs) | surv_specs == ""))
      stop("all elements of surv_specs must be filled in")
    
    unknown_type <- 
      !(surv_specs$type %in% c("PFS", "pfs", "OS", "os", "ToT"))
    if(any(unknown_type))
      stop(
        "only types 'PFS', 'OS', and 'ToT' are allowed; unknown type",
        plur(sum(unknown_type)),
        ": ",
        paste(surv_specs$type[unknown_type], collapse = ", ")
      )
    
    surv_specs <- 
      surv_specs %>% dplyr::arrange_(~ treatment, ~ type)
    
    
    surv_specs$isos <- surv_specs$type %in% c("OS", "os")
    surv_specs$ispfs <- surv_specs$type %in% c("PFS", "pfs")
    surv_specs$istot <- surv_specs$type == "ToT"
    
    surv_specs_check <- 
      surv_specs %>% dplyr::group_by_(~ treatment) %>%
      dplyr::summarize_(num_os = "sum(isos)",
                        num_pfs = "sum(ispfs)",
                        num_tot = "sum(istot)"
      )
    
    if(any(surv_specs_check$num_os != 1) |
       any(surv_specs_check$num_pfs != 1)){
      stop("each treatment must have exactly one PFS and one OS entry")
    }
    
    if(any(dups <- duplicated(surv_specs[, c("type", "treatment")]))){
      print(surv_specs[dups,])
      stop("survival fit specification can only have one row for fitting ",
           "OS, PFS, or ToT for a given treatment\n")
    }
    if(any(dups <- duplicated(surv_specs[, c("fit_directory", "fit_file", 
                                             "time_col", "censor_col")]))){
      print(surv_specs[dups,])
      stop("can only specify a given data file in a given directory with ",
           "the same time column and censoring column for one fit")
    }
    surv_specs
  }
