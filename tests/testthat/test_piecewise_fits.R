context("piecewise fitting functions") 
my_pe <- data.frame(start = c(0, 100, 200), hazard = c(0.003, 0.010, 0.006))

 set.seed(200)
  data.3piece <- 
   data.frame(time = r_piecewise_exponential_builder(my_pe)(1000),
              status = 1, group = 1)
 try_dists <- list(
 					rbind(
 						cbind("exp","exp"),
 						cbind("gamma","exp")
 					),
					rbind(
 						cbind("exp","exp","exp"),
						cbind("gamma","exp","exp")
 					))

 res <- f_find_best_piecewise_survival_models(num_pieces = c(2, 3), 
 			  dists = try_dists, build_fits = TRUE,
          survdata = data.3piece, time_col_name = "time", 
          censor_col_name = "status", treatment_col_name = "group", 
          fit_indiv_groups = FALSE, num_cand = 5, fixed_bp = TRUE, 
          min_bp_dist = 30, include_breaks = 200, eval_metric="BIC")[[1]]

 test_that("correct answer",
           {
             expect_equal(res$metric,
                              data.frame(num_pieces = as.numeric(2:3),
                                         metric = c(11439.2667241, 
                                                    11130.4427765))
                              )
             expect_equal(res$fits[[1]]$model_def,
                              data.frame(start = c(0, 200),
                                         dist_name = c("gamma", "exp"),
                                         shape = c(2.62499352, NA),
                                         rate = c(0.02091745542, 0.00659770225),
                                         stringsAsFactors = FALSE
                                         )
                              )
             expect_equal(res$fits[[2]]$model_def,
                              data.frame(start = c(0, 107, 200),
                                         dist_name = c("gamma", "exp", "exp"),
                                         shape = c(1.446698229, NA, NA), 
                                         rate = c(0.00693475248,
                                                  0.0211723636,
                                                  0.0065977429),
                                         stringsAsFactors = FALSE
                                         )
                              )
           })
 
 
 test_that("getting sets of breakpoints works",
           {
             times <- 1:100
             check_breaks <- 5
             include_breaks <- 50
              min_bp_dist <- 10
              num_breaks <- 2
              num_cand <- 10
              fixed_bp <- TRUE
              expect_error(get_breakpoints("times",
                               num_cand,
                               check_breaks,
                               include_breaks,
                               min_bp_dist,
                               num_breaks,
                               fixed_bp),
                           "is.numeric(times) is not TRUE",
                           fixed = TRUE
              )
              expect_error(get_breakpoints(times,
                                           "num_cand",
                                           check_breaks,
                                           include_breaks,
                                           min_bp_dist,
                                           num_breaks,
                                           fixed_bp),
                           "is.numeric(num_cand) is not TRUE",
                           fixed = TRUE
              )
              expect_error(get_breakpoints(times,
                                           num_cand,
                                           "check_breaks",
                                           include_breaks,
                                           min_bp_dist,
                                           num_breaks,
                                           fixed_bp),
                           "is.numeric(check_breaks) is not TRUE",
                           fixed = TRUE
              )
              expect_error(get_breakpoints(times,
                                           num_cand,
                                           check_breaks,
                                           "include_breaks",
                                           min_bp_dist,
                                           num_breaks,
                                           fixed_bp),
                           "is.numeric(include_breaks) is not TRUE",
                           fixed = TRUE
              )
              expect_error(get_breakpoints(times,
                                           num_cand,
                                           check_breaks,
                                           include_breaks,
                                           "min_bp_dist",
                                           num_breaks,
                                           fixed_bp),
                           "is.numeric(min_bp_dist) is not TRUE",
                           fixed = TRUE
              )
              expect_error(get_breakpoints(times,
                                           num_cand,
                                           check_breaks,
                                           include_breaks,
                                           min_bp_dist,
                                           "num_breaks",
                                           fixed_bp),
                           "is.numeric(num_breaks) is not TRUE",
                           fixed = TRUE
              )
              expect_error(get_breakpoints(times,
                                           num_cand,
                                           check_breaks,
                                           include_breaks,
                                           min_bp_dist = 250,
                                           num_breaks = 5,
                                           fixed_bp),
                           "no set of breakpoints meets the criteria",
                           fixed = TRUE
              )
              
              
           })