#'
#' Find best piecewise models for different numbers of pieces.
#'
#' @param num_pieces numbers of pieces to check (for example, 2:4)
#' @param dists a list of matrices with the distributions to be considered 
#'   in the different time periods using \code{flexsurvreg} - with one row for
#' 		each set of fits and one column for each time period. The number of elements
#' 		in the list must match num_pieces.
#' @param fit_indiv_groups Should treatment groups be fit individually?  Default TRUE.
#' @param build_fits should a fit with more pieces build on the
#'   breaks from fits with fewer pieces?
#' @param ... arguments to be passed to \code{f_find_best_with_fixed_num_pieces}
#'
#' @details
#'  \code{build_fits} has meaning only when num_pieces has
#'    length > 1 (that is, we are checking models with, say, 2 pieces
#'    and also models with 3 pieces).  
#'    
#'    \code{build_fits = FALSE} means that
#'    we start from scratch for each new number of pieces.
#'    
#'    \code{build_fits = TRUE} means that the breakpoints of the best smaller
#'    model are included as required breakpoints of a larger model.  For example,
#'    if the breakpoint of a 2-piece model is at 100, then, when fitting 3-piece
#'    models, only models with one breakpoint of 100 and one other breakpoint will
#'    be considered.  A model with breakpoints at 75 and 150 would not be
#'    considered, because it does not include a breakpoint at 100.   This
#'    substantially reduces the search space for larger models, which saves time,
#'    but could in theory miss an even better model.  In practice, using
#'    \code{build_fits = TRUE} seems to work well.
#'
#' @return  A list of lists, with each sublist pertaining to a particular group.
#' 	 Each sublist will have two elements:  \code{metric}, a data frame giving
#'   the evaluation metric (either AIC, BIC, or m2LL - depending on the selection set for 
#' 	 the eval_metric parameter) for the best fit with each number of pieces, and \code{fits},
#'    a list with the best fit for each number of pieces.   Finally, each `fit` object
#'    has two elements:  `model_def` and `fit`.   `model_def` is a human-readable
#'    description of the model, and `fit` is the corresponding `flexsurvreg` object.
#' @export
#'
#' @examples
#' my_pe <- data.frame(start = c(0, 100, 200), hazard = c(0.003, 0.010, 0.006))
#' data.3piece <- 
#'   data.frame(time = r_piecewise_exponential_builder(my_pe)(1000),
#'              status = 1, group = 1)
#' try_dists <- list(
#' 					cbind("exp"),
#' 					rbind(
#' 						cbind("exp","exp"),
#' 						cbind("gamma","exp")
#' 					),
#'					rbind(
#' 						cbind("exp","exp","exp"),
#'						cbind("gamma","exp","exp")
#' 					))
#'            
#' \dontrun{
#' res <- f_find_best_piecewise_survival_models(num_pieces = c(1, 2, 3), 
#' 			dists = try_dists, build_fits = TRUE,
#'          survdata = data.3piece, time_col_name = "time", 
#'          censor_col_name = "status", treatment_col_name = "group", 
#'          fit_indiv_groups = FALSE, num_cand = 10, fixed_bp = TRUE, 
#'          min_bp_dist = 30, include_breaks = 200, eval_metric="BIC")
#' }
#' 
#' my_pe <- data.frame(start = c(0, 100, 200), hazard = c(0.003, 0.010, 0.006))
#' my_pe2 <- data.frame(start = c(0, 50, 170), hazard = c(0.006, 0.005, 0.012))
#' try_dists <- list(
#' 					cbind("exp"),
#' 					rbind(
#' 						cbind("exp","exp"),
#' 						cbind("gamma","exp")
#' 					),
#'					rbind(
#' 						cbind("exp","exp","exp"),
#'						cbind("gamma","exp","exp")
#' 					))
#' data.3piece.1 <- 
#' 		data.frame(time = r_piecewise_exponential_builder(my_pe)(1000),
#' 				status = 1, group = 1)
#' data.3piece.2 <- 
#' 		data.frame(time = r_piecewise_exponential_builder(my_pe2)(1000),
#' 				status = 1, group = 2)
#' data.3piece <- rbind.data.frame(data.3piece.1, data.3piece.2)
#' \dontrun{
#' all_res <- f_find_best_piecewise_survival_models(num_pieces = c(1, 2, 3), 
#' 		dists = try_dists, fit_indiv_groups = TRUE, build_fits = TRUE,
#' 		survdata = data.3piece, time_col_name = "time", censor_col_name = "status",
#' 		treatment_col_name = "group", num_cand = 10, fixed_bp = TRUE, 
#' 		min_bp_dist = 30, include_breaks = 200, eval_metric="BIC")
#' }

f_find_best_piecewise_survival_models <-
		function(num_pieces, dists, fit_indiv_groups = TRUE, build_fits = TRUE, ...) {
	
	num_pieces <- sort(num_pieces)
	fits <- list()
	args_for_fit <- list(...)
	
	#check that dists list has length equal to num_pieces
	if(length(dists) != length(num_pieces))
		stop("dists list must have same length as num_pieces.")
	
	## make a list of the subgroups, and add a group of all together
	survdata <- args_for_fit[["survdata"]]
	treatment_col_name <- args_for_fit[["treatment_col_name"]]
	unique_groups <- as.character(unique(survdata[,treatment_col_name]))
	groups_list <- c(list(unique_groups))
	names(groups_list) <- c("all")
	attachNamespace("flexsurv")
#	if(fit_indiv_groups)
#	{
#		groups_list <- c(as.list(unique_groups), list(unique_groups))
#		names(groups_list) <- c(unique_groups, "all")
#	}
  groups_list <- as.list(unique_groups)
  names(groups_list) <- unique_groups
  
	conditions <- names(groups_list)
	all_res <- 
			lapply(1:length(conditions), function (j) {
				
				args_for_fit$which_groups <- as.vector(groups_list[[ conditions[j] ]])
	      
				for (i in seq(along = num_pieces)) {
					args_for_fit$num_breaks <- num_pieces[i] - 1
					args_for_fit$dists <- dists[[i]]
		      args_for_fit$best_only <- TRUE
					fits[[i]] <-
							try(
							  do.call("f_find_best_with_fixed_num_pieces",
											args_for_fit),
									silent = TRUE
							)
					if (build_fits & !inherits(fits[[i]], "try-error")) {
						args_for_fit$include_breaks <-
								sort(unique(c(
										args_for_fit$include_breaks,
										fits[[i]]$multi_spec$start
								)))
					}
				}
	      wanted_metric <- args_for_fit[["eval_metric"]]
				metric <- 
	        sapply(lapply(fits, "[[", "fit"),
	               function(x){
	                 if(!inherits(x, "try-error"))
	                    x[[wanted_metric]]
	                 else
	                   NA
	               })

				fit_summary <- data.frame(num_pieces = num_pieces,
						metric = metric)
	
				list(metric = fit_summary, fits = fits)
			})

	names(all_res) <- names(groups_list)
	all_res
}

#'
#' Find the best piecewise survival models with a given number of pieces.
#'
#' @param survdata Survival data to be used.
#' @param dists A matrix with the distributions to be considered 
#'   in the different time periods using
#'     \code{flexsurvreg} - with one row for each set and one column for
#' 		each time period. By default, this function fits a piecewise 
#' 		exponential distribution with 2 pieces.
#' @param treatment_col_name Name of the column in survdata that holds the
#'     treatment group names (for example "control", "treatment", "dose1", and
#'     so on).
#' @param time_col_name Name of the column in survdata with event times.
#' @param censor_col_name Name of the column in survdata with censorship
#'     indicators.   0:  event observed; 1:  censored.
#' @param covariate_col_names  Not yet implemented
#' @param which_groups Which treatment groups should be fit separately?
#' @param num_breaks How many breakpoints should there be in the
#'   piecewise fit? By default, set as 1.
#' @param num_cand number of candidate breakpoints checked.
#' @param fixed_bp Should breakpoints be considered fixed (the default) or
#'   can they be optimized over?
#' @param min_bp_dist The minimum distance between adjacent breakpoints.
#' @param check_breaks points that should be included in the set of
#'   potential breaks (but may not end up being included).
#' @param include_breaks breaks that must be included in any model.
#' @param eval_metric the evaluation metric to use for finding the best model. 
#'   Supports AIC, BIC, and -2*log likelihood (m2LL).
#' @param best_only should only the best fit for each condition
#'   be returned?
#' @param ... additional arguments to pass to flexsurvreg.
#' @export
#' @details  Current implementation supports any specified number of breakpoints.
#'   The evaluation metric can be specified as AIC, BIC, or m2LL.
#'
#'   \code{check_breaks} and \code{include_breaks} are only active
#'   when \code{fixed_bp = TRUE}.   In addition, \code{min_bp_dist}
#'   overrides \code{include_breaks}:  if you include breaks in
#'   \code{include_breaks} that are too close according to \code{min_bp_dist},
#'   they will be removed.
#'
#'   When \code{fixed_bp = TRUE}, breakpoints are static. Different sets of
#'   breakpoints are created (subject to the restrictions in \code{min_bp_dist}),
#'   a model is fit for each set of breakpoints, and quality measures are
#'   calculated and returned.  In this case we are optimizing only  the
#'   parameters of distribution in the various pieces, and only selecting
#'   from a pre-specified list of possible sets of breakpoints.
#'   When \code{fixed_bp = FALSE}, the breakpoints themselves are also
#'   optimized over.   However, because the functions may not have continuous
#'   derivatives at the breakpoints, this often does not converge well.
#'
#'   The requested number \code{num_cand} of candidate breakpoints
#'   are evenly-spaced percentiles of times in the data).  Higher values
#'   of \code{num_cand} gives a finer mesh of potential breakpoints, and
#'   therefore increases the chance of finding a breakpoint very close
#'   to the true breakpoint.   \code{num_cand} does
#'   not include any breaks in \code{check_breaks} or \code{include_breaks},
#'   which will be added to the set of candidate breakpoints.
#'   0 is always added as a required breakpoint; it need not be specified.
#'
#'   To not have a minimum distance between breakpoints,
#'   use \code{min_bp_dist = Inf}.
#'   This may result in checking potential fits with very
#'   closely-spaced breakpoints.
#'
#' @return If \code{best_only = FALSE}, a list with two components:
#' \itemize{
  #' \item{\code{fits}, with one element for each fit.  Each element is
  #' in turn a list with components \code{model_def} and \code{fit}.}
  #' \item{\code{quality}, a data frame with one row for each fit,
  #' and the various quality metrics in columns.}
#' }
#' If \code{best_only = TRUE}, only the element of fit corresponding
#' to the lowest fit metric will be returned.
#' 
#' @examples 
#' my_pe <- data.frame(start = c(0, 100, 200), hazard = c(0.003, 0.010, 0.006))
#' data.3piece <- 
#'   data.frame(time = r_piecewise_exponential_builder(my_pe)(1000),
#'              status = 1, group = 1)           
#' \dontrun{
#' res <- f_find_best_with_fixed_num_pieces(num_breaks = 2, 
#'          survdata = data.3piece, time_col_name = "time", 
#'          censor_col_name = "status", treatment_col_name = "group", 
#'          fit_indiv_groups = FALSE, num_cand = 10, 
#'          fixed_bp = TRUE, min_bp_dist = 30, eval_metric="BIC")
#' }
#' 
f_find_best_with_fixed_num_pieces <-
		function(survdata,
				dists = cbind("exp","exp"),
				time_col_name,
				censor_col_name,
				treatment_col_name,
				covariate_col_names = NULL,
				which_groups,
				num_breaks = 1,
				num_cand = 1,
				fixed_bp = TRUE,
				min_bp_dist = 30,
				check_breaks = numeric(0),
				include_breaks = numeric(0),
				eval_metric = "AIC",
				best_only = TRUE,
				...) {
	## TODO:  covariates
	
	## error checking
	if(!requireNamespace("flexsurv", quietly = TRUE))
		stop("flexsurv packaged needed to fit find best piecewise models")
	stopifnot(time_col_name %in% names(survdata))
	stopifnot(censor_col_name %in% names(survdata))
	stopifnot(treatment_col_name %in% names(survdata))
	stopifnot(all(which_groups %in% unique(survdata[[treatment_col_name]])))
	stopifnot(is.null(covariate_col_names) |
					all(covariate_col_names %in% names(survdata)))
	stopifnot(is.numeric(survdata[, time_col_name]))
	stopifnot(is.numeric(check_breaks))
	stopifnot(is.numeric(include_breaks))
	stopifnot(is.numeric(min_bp_dist))
	stopifnot(is.numeric(num_breaks))
	stopifnot(is.numeric(num_cand))
	stopifnot(eval_metric %in% c("AIC","BIC","m2LL"))
	
	#check that dist matrix has same number of columns as number of pieces
	if(ncol(dists) != (num_breaks + 1))
		stop("dists does not have the same number of columns as there are model pieces.")
	
	survdata <- survdata[survdata[[treatment_col_name]] %in% which_groups, ]

	#get breakpoints
	these_bp <-
			get_breakpoints(
					survdata[, time_col_name],
					num_cand,
					check_breaks,
					include_breaks,
					min_bp_dist,
					num_breaks,
					fixed_bp
			)
	bp_and_dists <- expand.grid(1:nrow(these_bp),
	                            1:nrow(dists))
  
	bp_and_dists_return <- cbind(these_bp[bp_and_dists[, 1], ],
	                             dists[bp_and_dists[, 2], ])
	
	## special distribution definition for flexsurv
	this_formula <- stats::formula(paste("survival::Surv(", time_col_name,
					", ", censor_col_name, ") ~ 1",
					sep = ""))
	
	args <- list(
			formula = this_formula,
			data = survdata
	)
	
	args <- c(args, list(...))
	
	fit_quality_across_dists <-
	  data.frame(
	    set = rep(as.numeric(NA), nrow(dists)),
	    AIC = rep(as.numeric(NA), nrow(dists)),
	    BIC = rep(as.numeric(NA), nrow(dists)),
	    m2LL = rep(as.numeric(NA), nrow(dists))
	  )
	
  dist_fits <- list()
  dist_fit_breaks <- list()
  for(dist_index in 1:nrow(dists)){
    use_dists <- dists[dist_index,]
    
    best_fit_across_bp <- 
      get_best_fit_across_breakpoints(args, use_dists, these_bp, fixed_bp,
                                      survdata[, time_col_name], eval_metric)
    dist_fits[[dist_index]] <- list(model_def = best_fit_across_bp$struc,
                                    fit = best_fit_across_bp$fit)
    dist_fit_breaks[[dist_index]] <- best_fit_across_bp$breakpoints

    fit_quality_across_dists[dist_index, "set"] <- dist_index 
    fit_quality_across_dists[dist_index, "AIC"] <- best_fit_across_bp$fit$AIC
    fit_quality_across_dists[dist_index, "BIC"] <- best_fit_across_bp$fit$BIC
    fit_quality_across_dists[dist_index, "m2LL"] <- best_fit_across_bp$fit$m2LL
  }
  if(best_only){
  	best_model_ind_across_dists <- 
  	  which.min(fit_quality_across_dists[, eval_metric])
  	
  	
  	return(dist_fits[[best_model_ind_across_dists]])
  }
  else{
    return(list(fits = dist_fits,
                quality = fit_quality_across_dists))
  }
	
}

#' Return a set of breakpoints for fitting a piecewise survival model
#'
#' @param times times from survival data
#' @param num_cand number of candidate breakpoints checked (selected at
#'   evenly-spaced percentiles of times in survdata).  This number does
#'   not include any breaks in \code{check_breaks} or \code{include_breaks}.
#' @param check_breaks points that should be included in the set of
#'   potential breaks (but may not end up being included).
#' @param include_breaks breaks that must be included in any solution.
#' @param min_bp_dist The minimum distance between adjacent breakpoints.
#' @param num_breaks How many breakpoints should there be in the
#'   piecewise survival fit?
#' @param fixed_bp Should breakpoints be considered fixed (the default) or
#'   can they be optimized over?
#' @return a matrix with one row for each set of breakpoints.
#'
get_breakpoints <-
		function(times,
				num_cand,
				check_breaks,
				include_breaks,
				min_bp_dist,
				num_breaks,
				fixed_bp) {
	
	#argument checks
	stopifnot(is.numeric(times))
	stopifnot(is.numeric(check_breaks))
	stopifnot(is.numeric(include_breaks))
	stopifnot(is.numeric(min_bp_dist))
	stopifnot(is.numeric(num_breaks))
	stopifnot(is.numeric(num_cand))
	stopifnot(is.logical(fixed_bp))
	
	## get the set of breakpoints - first the percentile points
	use_percentiles <-
			seq(from = 0,
					to = 1,
					length.out = num_cand + 2)
	bp_set <- stats::quantile(times,
			                      probs = use_percentiles[-c(1, num_cand + 2)])
	bp_set <- unique(round(bp_set))
	
	if (fixed_bp) {
		if (num_breaks == 0) {
			these_bp <- matrix(0, nrow = 1, ncol = 1)
		}
		else{
			include_breaks <- setdiff(include_breaks, 0)
			## set up all combinations of breakpoints ...
			bp_set <-
					sort(unique(c(
											bp_set, check_breaks, include_breaks
									)))
			to_expand <- rep(list(bp_set), num_breaks)
			these_bp <- do.call("expand.grid", to_expand)
			
			## and filter down to the ones we want
			all_required_breaks <- apply(these_bp, 1,
					function(x) {
						all(include_breaks %in% x)
					})
			these_bp <- subset(these_bp, all_required_breaks)
			
			these_bp <- cbind(0, as.matrix(these_bp))
			
			min_distance_observed <-
					apply(these_bp, 1, function(x) {
								all(diff(x) > min_bp_dist)
							})
			these_bp <- subset(these_bp, min_distance_observed)
		}
	}
	else{
		these_bp <- matrix(c(0, bp_set), nrow = 1)
	}
	
	if (nrow(these_bp) == 0)
		stop("no set of breakpoints meets the criteria")
	
	dimnames(these_bp) <- list(NULL, NULL)
	these_bp
		}



#' Get best fit across sets of breakpoints
#'
#' @param args arguments for fitting passed from \code{f_find_best_with_fixed_num_pieces}
#' @param use_dists the distributions to be used for the pieces between breakpoints
#' @param these_bp sets of breakpoints to be compared across
#' @param fixed_bp Should breakpoints be considered fixed (the default) or
#'   can they be optimized over?
#' @param survtimes the survival times
#' @param eval_metric one of \code{AIC}, \code{BIC}, or \code{m2LL}
#'
#' @return a list with 3 components:
#'   \itemize{
#'   \item{\code{fit}; the fit object from \code{flexsurvreg};}
#'   \item{\code{breaks}; the corresponding breakpoints}; and
#'   \item{\code{params}; the parameters for the distributions in the piecewise fit}
#'   }
#'
get_best_fit_across_breakpoints <-
  function(args, use_dists, these_bp, fixed_bp, survtimes, eval_metric){

    stopifnot(eval_metric %in% c("AIC","BIC","m2LL"))
    
    bp_fit <- list()
    for(bp_index in 1:nrow(these_bp)){
      use_bp <- these_bp[bp_index, ]
      if (options()$heemod.verbose) message(use_bp)
      add_args <- get_survival_arg_list(these_bp = use_bp,
                                        dists = use_dists, 
                                        fixed_bp = fixed_bp,
                                        survtimes = survtimes)
      loadNamespace("flexsurv")
      suppressMessages(
        bp_fit[[bp_index]] <- do.call(flexsurv::flexsurvreg, c(args, add_args))
      )
   }
    bp_fit <- f_add_surv_fit_metrics(bp_fit)
    fit_quality_bp <- 
      do.call(rbind,
              lapply(bp_fit, function(x){
                      c(set = x$set, AIC = x$AIC,
                        BIC = x$BIC, m2LL = x$m2LL)
                    })
      )
    fit_quality_bp <- data.frame(set = 1:nrow(fit_quality_bp),
                              fit_quality_bp)
    best_bp_ind <- which.min(fit_quality_bp[, eval_metric])  
    
    params <- bp_fit[[best_bp_ind]]$opt$par
    optpars <- bp_fit[[best_bp_ind]]$optpars
    inv.trans <- bp_fit[[best_bp_ind]]$dlist$inv.transforms[optpars]
    inv.trans.params <- 
      sapply(seq(along = params), function(param_ind){
        inv.trans[[param_ind]](params[[param_ind]])
      })
    
    list(fit = bp_fit[[best_bp_ind]],
         breakpoints = as.vector(these_bp[best_bp_ind,]),
         struc = build_pieces_def_structure(dists = use_dists, 
                                            breaks = as.vector(these_bp[best_bp_ind,]),
                                            params = inv.trans.params)
    )
  }

#' Get a list of arguments for fitting a piecewise survival curve
#'
#' @param these_bp break points.
#' @param dists Sequence of distributions.
#' @param fixed_bp Should the break points be considered fixed?
#' @param survtimes Survival times for the analysis.
#' @return a list of arguments to be used in 
#'    \code{\link[flexsurv]{flexsurvreg}}
#'
get_survival_arg_list <-
  function(these_bp, dists, fixed_bp, survtimes) {
    num_breaks <- length(these_bp) - 1  

      if (length(dists) != (num_breaks + 1))
        stop(
          "length of dists must be one greater than number of breaks ",
          "(except for special case of piecewise exponential)"
        )
      ## piecewise distribution with multiple distributions
      hfn_args_list <- haz_fn_args(dists)
      num_params <- length(dists) + sum(hfn_args_list$num_hfn_args)
      Hfn <-
        flexsurv::unroll.function(H_piecewise_multi,
                                  pieces_def = 0:(num_params - 1))
      hfn <-
        flexsurv::unroll.function(h_piecewise_multi,
                                  pieces_def = 0:(num_params - 1))
      custom_fss <-
        get_custom_fss_piecewise_multi(dists = dists,
                                       fixed_bp = fixed_bp)
      inits <- get_multi_dist_inits(these_bp, dists = dists,
                                         t = survtimes)
      aux <- list(dists = dists)
      
      fixed_args <- 
        if(fixed_bp)
          1:length(dists)
      else 
        1

    list(dist = custom_fss,
         dfns = list(h = hfn, H = Hfn),
         fixed = fixed_args,
         inits = inits,
         aux = aux
    )
  }
