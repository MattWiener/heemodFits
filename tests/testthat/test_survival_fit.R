context("fitting survival models")

data <- data.frame(time = rexp(100, rate = 0.01),
                   status = rep(1, 100),
                   group = rep(1, 100),
                   ITTFL = rep("Y", 100))

test_that("input errors are caught", {
  expect_error(
  f_fit_survival_models(
    data,
    dist,
    time_col_name = "time2",
    censor_col_name = "status",
    fit_indiv_groups = F
  ),
  "time_col_name"
)
expect_error(
  f_fit_survival_models(
    data,
    dist,
    time_col_name = "time",
    censor_col_name = "status2",
    treatment_col_name = "group",
    fit_indiv_groups = F
  ),
  "censor_col_name"
)
expect_error(
  f_fit_survival_models(
    data,
    dist,
    time_col_name = "time",
    censor_col_name = "status",
    treatment_col_name = "group2",
    fit_indiv_groups = F
  ),
  "treatment_col_name"
)

expect_error(
  f_fit_survival_models(
    as.list(data),
    dist,
    time_col_name = "time",
    censor_col_name = "status",
    treatment_col_name = "group2",
    fit_indiv_groups = F
  ),
  "is.data.frame"
)

expect_error(
  f_fit_survival_models(
    data[0, ],
    dist,
    time_col_name = "time",
    censor_col_name = "status",
    treatment_col_name = "group2",
    fit_indiv_groups = F
  ),
  "nrow(survdata)", fixed = TRUE
)

})

test_that("reading set definitions works",
          {
            example_1 <- 
              get_set_definitions(system.file("surv", package = "heemodFits"),
                                  "set_definitions_1.csv")
            expect_identical(example_1,
                             data.frame(type = rep(c("PFS", "OS"), each = 4),
                                        treatment = rep(c("fake_treatment", "not_real"),4),
                                        set_name = rep(rep(c("all", "time.gt.100"), each = 2),
                                                       2),
                                        condition = rep(rep(c("TRUE", "time > 100"), each = 2),
                                                        2),
                                        time_subtract = rep(rep(c(0, 100), each = 2),
                                                            2),
                                        stringsAsFactors = FALSE
                                        )
                             )
            example_2 <-
              get_set_definitions(system.file("surv", package = "heemodFits"),
                                  "set_definitions_2.csv")
            expect_identical(example_2,
                             data.frame(type = c("PFS", "OS"), 
                                        treatment = rep("fake_treatment", 2),
                                        set_name = rep("all", 2), 
                                        condition = rep("TRUE", 2),
                                        time_subtract = c(0, 0),
                                        stringsAsFactors = FALSE)
                             )
            expect_error(get_set_definitions(system.file("surv", package = "heemodFits"),
                                             "set_definitions_error_1.csv"),
                         "set_definitions file missing column(s): treatment, set_name",
                         fixed = TRUE
                         )
            expect_error(get_set_definitions(system.file("surv", package = "heemodFits"),
                                             "set_definitions_error_2.csv"),
                         "bad offset in set_definitions",
                         fixed = TRUE
            )
            expect_warning(get_set_definitions(system.file("surv", package = "heemodFits"),
                                               "set_definitions_error_3.csv"),
                           "value of 'time_subtract' does not appear in 'condition'",
                           fixed = TRUE)
            expect_error(get_set_definitions(system.file("surv", package = "heemodFits"),
                                             "set_definitions_3.csv"),
                         "multiple definitions of some sets"
                         )
            }
          )


test_that("getting survival inputs works",
          {
           ok_surv_info <- 
             read_file(system.file("surv/survival_info.csv", 
                                   package = "heemodFits"))
             check_survival_specs(ok_surv_info)
           for(i in 1:(ncol(ok_surv_info) - 3)){
             expect_error(check_survival_specs(ok_surv_info[, -i]),
                          "missing names")
           }
            for(i in ncol(ok_surv_info) - 2:0){
              expect_warning(check_survival_specs(ok_surv_info[, -i]),
                           "not defined in surv_specs")
            }
           surv_info_extra <- ok_surv_info
           surv_info_extra$extra_col <- 5
           expect_error(check_survival_specs(surv_info_extra),
                        "extra names")
           expect_error(check_survival_specs(surv_info_extra),
                        "extra_col")
           surv_info_extra_row <- ok_surv_info[c(1, 1:nrow(ok_surv_info)),]
           expect_error(check_survival_specs(surv_info_extra_row),
                     "exactly one PFS and one OS entry")
           surv_info_dup <- ok_surv_info
           surv_info_dup[1, c("fit_directory", "fit_file", "time_col", "censor_col")] <- 
             surv_info_dup[2, c("fit_directory", "fit_file", "time_col", "censor_col")]
           expect_error(capture.output(check_survival_specs(surv_info_dup)),
                        "same time column and censoring column")
           surv_info_wrong_type <- ok_surv_info
           surv_info_wrong_type[1, "type"] <- "oops"
           expect_error(check_survival_specs(surv_info_wrong_type),
                        "only types 'PFS', 'OS', and 'ToT' are allowed; unknown type: oops",
                        fixed = TRUE
           )
           surv_info_na <- ok_surv_info
           surv_info_na[1, "time_col"] <- NA
           expect_error(check_survival_specs(surv_info_na),
                        "all elements of surv_specs must be filled in")
           
           }
          )

test_that("fitting works (including with subsets)",
          {
            location <- system.file("surv", package = "heemodFits")
            ok_surv_info <- 
              read_file(system.file("surv/survival_info.csv", 
                                    package = "heemodFits"))
            not_ok_surv_info <- 
              ok_surv_info
            not_ok_surv_info$censor_col[1] <- "status2"
            
            expect_error(
              survival_from_data(location = location,
                                          survival_specs = not_ok_surv_info,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE), 
              "censoring status column 'status2' does not exist in data file"
              
            )
              full_these_fits <- 
              survival_from_data(location = location,
                                 survival_specs = ok_surv_info,
                                 dists = c("exp", "weibull"),
                                 save_fits = FALSE) 
              these_fits <- dplyr::filter_(full_these_fits,
                                           ~treatment != "all")
            expect_identical(names(these_fits), 
                             c("type", "treatment", "set_name",
                               "dist", "fit", "set_def",
                               "time_subtract"))
            expect_identical(these_fits$dist,
                             rep(c("exp", "weibull", "km"), 10))
            expect_identical(sapply(these_fits$fit, class),
                             c(rep(c("flexsurvreg", "flexsurvreg", "survfit",
                                   "surv_shift", "surv_shift", "surv_shift"), 
                                 3),
                             c("flexsurvreg", "flexsurvreg", "survfit",
                               "flexsurvreg", "flexsurvreg", "survfit",
                               "surv_shift", "surv_shift", "surv_shift",
                               "flexsurvreg", "flexsurvreg", "survfit")))
            combos <- table(these_fits[, c("treatment", "set_def")])
            ## sorting to make sure things are in right order for tests -
            ##   otherwise sometimes had problems with different locales
            combos <- combos[c("A", "B"),]
            combos <- combos[, c("biomarker > 0.5","time > 50", "TRUE")]
            expect_identical(dimnames(combos),
                             list(treatment = c("A", "B"),
                                  set_def = c("biomarker > 0.5","time > 50", "TRUE")))
            expect_identical(as.numeric(combos),
                             c(0, 6, 6, 6, 6, 6))
            metrics <- extract_surv_fit_metrics(these_fits)
            expect_identical(names(metrics),
                             c("type", "treatment", "set_name", "dist", "fit",
                               "set_def", "time_subtract", "AIC", "BIC", "m2LL"))
            expect_equal(nrow(metrics), 20)
            expect_identical(round(metrics[1, c("AIC", "BIC", "m2LL")], 3),
                             tibble::tribble(~AIC, ~BIC, ~m2LL,
                                             345.293, 347.205, 343.293)
            
                             )
            ## now test with absolute path
            abs_path_surv_info <- ok_surv_info
            abs_path_surv_info$data_directory <-
              file.path(location, ok_surv_info$data_directory)
              abs_path_fits <- 
              survival_from_data(location = location,
                                          survival_specs = abs_path_surv_info,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE)#,
                                          # use_envir = new.env())
            abs_path_fits <- dplyr::filter_(abs_path_fits, 
                                           ~ treatment != "all")
            metrics <- extract_surv_fit_metrics(abs_path_fits)
            expect_identical(names(metrics),
                             c("type", "treatment", "set_name", "dist", "fit",
                               "set_def", "time_subtract", "AIC", "BIC", "m2LL"))
            expect_equal(nrow(metrics), 20)
            expect_identical(round(metrics[1, c("AIC", "BIC", "m2LL")], 3),
                             tibble::tribble(~AIC, ~BIC, ~m2LL,
                                             345.293, 347.205, 343.293)
                             
            )
            
            ## make sure metrics works calling with just one row
            ##metrics <- extract_surv_fit_metrics(these_fits[[1]][1,])
            metrics <- extract_surv_fit_metrics(these_fits[1,])
            expect_identical(round(metrics[, c("AIC", "BIC", "m2LL")], 3),
                             tibble::tribble(~AIC, ~BIC, ~m2LL,
                                             345.293, 347.205, 343.293)
                             )
            ## make sure we get an error if we specify incorrect event codes
            eventcode_surv_info_error <-
              read_file(system.file("surv/survival_info_eventcode_error.csv", 
                                             package = "heemodFits"))
            expect_error(
              survival_from_data(location = location,
                                          survival_specs = eventcode_surv_info_error,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE),
              "all values should be either 'event' (for events) or 'censor' (for censoring)",
              fixed = TRUE
            )
            
            ## make sure we run correctly if we specify correct event codes
            eventcode_surv_info<-
              read_file(system.file("surv/survival_info_eventcode.csv", 
                                             package = "heemodFits"))
            eventcode_fits <-
              survival_from_data(location = location,
                                          survival_specs = eventcode_surv_info,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE) 
            eventcode_fits <-
              dplyr::filter_(eventcode_fits, ~treatment != "all")
            ## this should have same results as earlier fits
            expect_identical(lapply(these_fits$fit, heemod::compute_surv, time = c(45:55)),
                             lapply(eventcode_fits$fit, heemod::compute_surv, time = c(45:55))
            )
            ## check that if we designate subsets by type, the ones
            ##    we leave out don't show up (no GT50 for OS)
            subset_fits_by_type <- 
              survival_from_data(location = location,
                                          survival_specs = ok_surv_info,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE,
                                          set_definitions = "set_def_pfs_os.csv")
            expect_equal(unique(subset_fits_by_type[, 1:3]),
                             tibble::tribble(
                               ~type, ~treatment, ~set_name,
                               "OS", "A", "all",
                               "PFS", "A", "all",
                               "PFS", "A", "GT50",
                               "OS", "B", "all",
                               "OS", "B", "B5",
                               "PFS", "B", "all",
                               "PFS", "B", "GT50",
                               "PFS", "B", "B5",
                               "OS", "all", "all",
                               "PFS", "all", "all",
                               "PFS", "all", "GT50"
                             )
            )
            expect_error(
              suppressWarnings(survival_from_data(location = location,
                                          survival_specs = ok_surv_info,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE,
                                          set_definitions = "set_def_pfs_os_error.csv")),
              "1 '0' value"
              )
            expect_error(
              suppressWarnings(survival_from_data(location = location,
                                          survival_specs = ok_surv_info,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE,
                                          set_definitions = "set_def_pfs_os_error_2.csv")),
              "negative value"
            )
            
          }
          )

test_that("we handle fitting errors",
          {
            ## some data that causes an error when
            ##  fitting with generalized gamma
            this_dat <- 
              tibble::tribble(
              ~time, ~event, ~trt,
              251,0,"A",
              250,0,"A",
              91,1,"A",
              173,1,"A",
              467,0,"A",
              382,0,"A",
              247,1,"A",
              8,1,"A",
              291,1,"A",
              49,1,"A",
              1,1,"A",
              2,1,"A",
              259,0,"A",
              288,0,"A",
              212,0,"A",
              45,1,"A",
              9,1,"A",
              257,0,"A",
              63,1,"A",
              69,1,"A",
              338,0,"A",
              292,0,"A",
              248,1,"A",
              1,1,"A",
              1,1,"A",
              73,1,"A",
              223,0,"A",
              1,1,"A",
              2,1,"A",
              378,0,"A",
              464,0,"A",
              295,0,"A",
              40,1,"A",
              2,1,"A",
              44,1,"A"
            )
            this_dat$ITTFL <- "Y"
            suppressMessages(
            fit_tib <- f_fit_survival_models(this_dat, 
                                         dist = c("exp", "weibull", "gengamma"),
                                         time_col_name = "time", 
                                         censor_col_name = "event", 
                                         treatment_col_name = "trt",
                                         fit_indiv_groups = FALSE)
            )
            # expect_equal(sapply(fit_tib$fit, class),
            #              c("flexsurvreg", "flexsurvreg", 
            #                "try-error", "survfit")
            #              )
            }
          )


