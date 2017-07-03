context("test convenience functions")

test_that("functions for AE values work",
          {
            AEs <- data.frame(
              treatment = c("A", "A", "B", "B"),
              ae = c("ae1", "ae2", "ae1", "ae2"),
              prob = c(0.1, 0.1, 0.2, 0),
              cost = c(100, 200, 100, 200),
              disutility = c(0, .1, 0, 0.1),
              stringsAsFactors = FALSE
            )
            
            res1 <- compute_vals_for_adv_ev(AEs)
            correct <- tibble::tribble( ~ treatment, ~ cost, ~ disutility,
                                        "A", 30, 1 / 100,
                                        "B", 20, 0)
            expect_equivalent(as.data.frame(res1),
                              as.data.frame(correct))
            expect_error(compute_vals_for_adv_ev(AEs[, c("treatment", "ae",
                                                         "prob")]),
                         "must have some column for a value")
            AEs2 <- AEs
            AEs2$prob <- c(0.1, 0.1, 0.2, NA)
            
            expect_warning(
              res2 <- compute_vals_for_adv_ev(AEs2),
              "some AE probabilities are missing in the table"
            )
            expect_identical(res1, res2)
            
            AEs3 <- AEs
            AEs3$cost <- c(1, 2, 3, 4)
            
            expect_error(compute_vals_for_adv_ev(AEs3),
                         "multiple values defined")
            
            AEs4 <- AEs
            names(AEs4) <- c("treatment", "AE", "prob", "cost")
            expect_error(compute_vals_for_adv_ev(AEs4),
                         "adverse events table must have columns")
            
            expect_equal(ae_val(AEs, "A", "cost"), 30)
            expect_equal(ae_val(AEs, "B", "cost"), 20)
            expect_equal(ae_val(AEs, "A", "disutility"), 0.01)
            expect_equal(ae_val(AEs, "B", "disutility"), 0)
            expect_error(ae_val(AEs, "C", "cost"),
                         "no AE information returned")
            expect_error(ae_val(AEs, "A", "novalue"),
                         "no column")
            
            AEs5 <- AEs
            AEs5$cost[c(3, 4)] <- NA
            expect_equal(ae_val(AEs5, "A", "cost"), 30)
            expect_equal(ae_val(AEs5, "B", "cost"), 20)
            expect_equal(ae_val(AEs5, "A", "disutility"), 0.01)
            expect_equal(ae_val(AEs5, "B", "disutility"), 0)
            
          })

test_that("is_dosing_period works",
          {
            expect_identical(
              is_dosing_period(
                N = 1:13,
                first = 4,
                then_every = 3,
                cap = 40
              ),
              c(
                TRUE,
                TRUE,
                TRUE,
                TRUE,
                FALSE,
                FALSE,
                TRUE,
                FALSE,
                FALSE,
                TRUE,
                FALSE,
                FALSE,
                TRUE
              )
            )
            expect_identical(
              is_dosing_period(
                N = 37:46,
                first = 4,
                then_every = 3,
                cap = 40
              ),
              c(
                TRUE,
                FALSE,
                FALSE,
                TRUE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE
              )
            )
            expect_identical(is_dosing_period(
              N = 1:8,
              first = 4,
              pattern = 0,
              cap = 40
            ),
            rep(c(TRUE, FALSE), each = 4))
            expect_error(
              is_dosing_period(
                N = 1:8,
                first = 4,
                then_every = -1,
                cap = 40
              ),
              "then_every cannot be negative"
            )
            ## demonstrating argument precedence rules
            expect_identical(
              is_dosing_period(
                N = 1:10,
                init = c(1, 0, 1),
                first = 3,
                then_every = 5
              ),
              c(TRUE, FALSE, TRUE, FALSE, FALSE,
                FALSE, FALSE, TRUE, FALSE, FALSE)
            )
            expect_identical(
              is_dosing_period(
                N = 1:10,
                init = numeric(0),
                pattern = c(1, 1, 0, 1, 0),
                then_every = 2
              ),
              c(TRUE, TRUE, FALSE, TRUE, FALSE,
                TRUE, TRUE, FALSE, TRUE, FALSE)
            )
          })

test_that("utility by time before death",
          {
            ProgFree <- round(1000 * exp(-0.4 * 0:24))
            Progressive <- round((1000 - ProgFree) * exp(-0.4 * 0:24))
            Death <- 1000 - ProgFree - Progressive
            state_names <-
              rep(c("ProgFree", "Progressive", "Death"), each = 25)
            
            counts <- data.frame(
              .strategy = rep("s1", 25),
              markov_cycle = 0:24,
              state_names = state_names,
              count = c(ProgFree, Progressive, Death)
            )
            class(counts) <- c("cycle_counts", class(counts))
            
            ## if utility is 1 in the cycle before death and 0 otherwise,
            ##   then utility should be equal to the number of people
            ##   about to die
            aa1 <- data.frame(until_lag = 1, util = 1)
            res1 <-
              utility_by_time_from_death(counts,
                                         util_before_death = aa1,
                                         util_long_before_death = 0)
            expect_identical(res1[-length(res1)],
                             diff(dplyr::filter(counts,
                                                state_names == "Death")$count))
            ## gets delayed by 1 if utility is 1 only two cycles before death
            aa2 <- data.frame(until_lag = 1:2, util = 0:1)
            res2 <-
              utility_by_time_from_death(counts,
                                         util_before_death = aa2,
                                         util_long_before_death = 0)
            expect_identical(res2[-(length(res2) - 1 + 0:1)],
                             diff(dplyr::filter(counts,
                                                state_names == "Death")$count)[-1])
            
            ## catch the error if not all the subjects die
            restricted_counts <- dplyr::filter_(counts, ~ markov_cycle <= 15)
            class(restricted_counts) <- c("cycle_counts", "data.frame")
            expect_error(
              utility_by_time_from_death(
                restricted_counts,
                util_before_death = aa2,
                util_long_before_death = 0
              ),
              "not all subjects reach the death state"
            )
            
            ## we should be able to replace util_long_before_death by
            ##  having a long lag time with the same utility.
            ## check once defining times longer than number of cycles
            
            util_df_1 <- data.frame(until_lag = c(4, 13, 26, 39, 52),
                                    util = c(.2, .3, .4, .5, .6))
            util_df_2 <- rbind(util_df_1,
                               c(until_lag = 1000, util = 0.7))
            res_1 <-
              utility_by_time_from_death(counts,
                                         util_before_death = util_df_1,
                                         util_long_before_death = 0.7)
            res_2 <-
              utility_by_time_from_death(counts,
                                         util_before_death = util_df_2,
                                         util_long_before_death = 0)
            expect_equal(res_1, res_2)
            
            util_df_1 <- data.frame(until_lag = c(2, 4, 6, 8, 10),
                                    util = c(.2, .3, .4, .5, .6))
            util_df_2 <- rbind(util_df_1,
                               c(until_lag = 1000, util = 0.7))
            res_1 <-
              utility_by_time_from_death(counts,
                                         util_before_death = util_df_1,
                                         util_long_before_death = 0.7)
            res_2 <-
              utility_by_time_from_death(counts,
                                         util_before_death = util_df_2,
                                         util_long_before_death = 0)
            expect_equal(res_1, res_2)
            
            
          })

test_that("finding least-cost combination of vials for a dose works",
          {
            units <- data.frame(
              trt = c("fake", "fake"),
              volume = c(40, 100),
              cost = c(1003.01, 2507.54)
            )
            expect_error(
              find_least_cost_partition(214.3,
                                        units,
                                        subset_col = "trt",
                                        subset_val = "fake"),
              "available_units must have columns 'size' and 'cost'"
            )
            names(units)[2] <- c("size")
            expect_equal(
              find_least_cost_partition(214.3,
                                        units,
                                        subset_col = "trt",
                                        subset_val = "fake"),
              data.frame(
                desired_dose = 214.3,
                used_dose = 220,
                waste = 5.7,
                cost = 5516.57,
                cost.no.waste = 5373.64
              )
            )
            
            expect_equal(
              find_least_cost_partition(
                c(214.3, 200),
                units,
                subset_col = "trt",
                subset_val = "fake"
              ),
              data.frame(
                desired_dose = c(214.3, 200),
                used_dose = c(220, 200),
                waste = c(5.7, 0.0),
                cost = c(5516.57, 5015.05),
                cost.no.waste = c(5373.64, 5015.05)
              )
            )
            
            expect_error(
              find_least_cost_partition(214.3,
                                        units,
                                        subset_col = "trt"),
              "subset_col and subset_val should either both be NULL"
            )
            expect_error(
              find_least_cost_partition(214.3,
                                        units,
                                        subset_val = "fake"),
              "subset_col and subset_val should either both be NULL"
            )
            expect_error(
              find_least_cost_partition(214.3,
                                        units,
                                        subset_col = "not_a_col",
                                        subset_val = "fake"),
              "is not a column of available_units"
            )
            expect_error(
              find_least_cost_partition(214.3,
                                        units,
                                        subset_col = "trt",
                                        subset_val = "not_there"),
              "does not contain the value"
            )
            expect_error(
              find_least_cost_partition(
                214.3,
                units,
                subset_col = "trt",
                subset_val = c("fake1", "fake2")
              ),
              "exactly one value"
            )
            expect_error(
              find_least_cost_partition(
                214.3,
                units,
                subset_col = c("trt1", "trt2"),
                subset_val = "fake"
              ),
              "exactly one column"
            )
            
          })

test_that("cost_at_right_time works",
          {
            expect_error(
              cost_at_right_time(
                cost = 5,
                cycle = 1:10,
                age = 50 + 1:10
              ),
              "must specify at least one of 'at_cost' and 'at_age'"
            )
            expect_error(
              cost_at_right_time(
                cost = 5,
                at_cycle = 5,
                at_age = 55
              ),
              "must specify at least one of 'cost' and 'age'"
            )
            expect_error(
              cost_at_right_time(
                cost = 5,
                cycle = 1:10,
                age = 50 + 1:11,
                at_cycle = 8
              ),
              "must have the same length"
            )
            expect_error(
              cost_at_right_time(
                cost = 5,
                cycle = 1:10,
                at_age = 55
              ),
              "at least one of both 'cycle' and 'at_cycle'"
            )
            expect_error(
              cost_at_right_time(
                cost = 5,
                age = 1:10,
                at_cycle = 5
              ),
              "at least one of both 'cycle' and 'at_cycle'"
            )
            expect_equal(cost_at_right_time(
              cost = 5,
              cycle = 1:10,
              at_cycle = 20
            ),
            rep(0, 10))
            expect_equal(cost_at_right_time(
              cost = 5,
              age = 1:10,
              at_age = 20
            ),
            rep(0, 10))
            expect_equal(cost_at_right_time(
              cost = 5,
              cycle = 1:10,
              at_cycle = c(3, 6)
            ),
            c(0, 0, 5, 0, 0, 5, 0, 0, 0, 0))
            expect_equal(cost_at_right_time(
              cost = 5,
              age = 1:10,
              at_age = c(4, 7)
            ),
            c(0, 0, 0, 5, 0, 0, 5, 0, 0, 0))
            expect_equal(
              cost_at_right_time(
                cost = 5,
                cycle = 1:10,
                age = 50 + 1:10,
                at_cycle = 3,
                at_age = 55
              ),
              c(0, 0, 5, 0, 5, 0, 0, 0, 0, 0)
            )
            
          })

test_that("intravenous dosing cost function works",
          {
            expect_equal(cost_iv_administration(0.5, 100, 20, prorate_first = TRUE),
                         50)
            expect_equal(cost_iv_administration(1.5, 100, 20, prorate_first = TRUE),
                         110)
            expect_equal(cost_iv_administration(0.5, 100, 20, prorate_first = FALSE),
                         100)
            expect_equal(cost_iv_administration(1.5, 100, 20, prorate_first = FALSE),
                         110)
          })

test_that("getting weighted costs works",
          {
            vialCost <- data.frame(treatment = "fake",
                                   size = 50,
                                   cost = 1000)
            expect_equal(
              weighted_dose_costs(
                "lnorm",
                params = c(meanlog = 4.3423559,
                           sdlog = 0.2116480),
                var_base = 50,
                dose_base = 100,
                dose_multiplier = 2,
                available_units = vialCost,
                subset_col = "treatment",
                subset_val = "fake",
                share_vials = FALSE
              ),
              weighted_lognorm_dose_costs(
                4.3423559,
                0.2116480,
                var_base = 50,
                dose_base = 100,
                dose_multiplier = 2,
                available_units = vialCost,
                subset_col = "treatment",
                subset_val = "fake",
                share_vials = FALSE
              )
            )
            expect_equal(round(
              weighted_lognorm_dose_costs(
                4.3423559,
                0.2116480,
                var_base = 50,
                dose_base = 100,
                dose_multiplier = 2,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE
              ),
              2
            ),
            3636.36)
            expect_equal(round(
              weighted_lognorm_dose_costs(
                4.3423559,
                0.2116480,
                var_base = 50,
                dose_base = 100,
                dose_multiplier = 2,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = TRUE
              ),
              2
            ),
            3140.94)
            
            ## testing that getting duplicate values for
            ##   parameters is OK
            expect_equal(
              round(
                weighted_dose_costs(
                  "lnorm",
                  params = list(meanlog = rep(4.3423559, 5),
                             sdlog = rep(0.2116480, 5)),
                  var_base = 50,
                  dose_base = 100,
                  dose_multiplier = 2,
                  available_units = vialCost,
                  subset_col = "treatment",
                  subset_val = "fake",
                  share_vials = FALSE
                ), 
                2),
              3636.36
            )
              
            
            
            vialCost <- data.frame(
              treatment = "fake",
              size = c(50, 250),
              cost = c(2521, 12297)
            )
            expect_equal(round(
              weighted_norm_dose_costs(
                1.85,
                0.25,
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE
              ),
              2
            ),
            30511.71)
            expect_equal(
              weighted_dose_costs(
                "norm",
                params = c(mean = 1.85, sd = 0.25),
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE
              ),
              weighted_norm_dose_costs(
                1.85,
                0.25,
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE
              )
            )
            
            expect_error(
              weighted_norm_dose_costs(
                1.85,
                0.25,
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE,
                qmin = 0,
                qmax = 0.99,
                by = 0.01
              ),
              "it does not make sense for 'qmin' to be <= 0",
              fixed = TRUE
            )
            expect_error(
              weighted_norm_dose_costs(
                1.85,
                0.25,
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE,
                qmin = 0.01,
                qmax = 1,
                by = 0.01
              ),
              "it does not make sense for 'qmax' to be >= 1",
              fixed = TRUE
            )
            expect_error(
              weighted_norm_dose_costs(
                1.85,
                0.25,
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE,
                qmin = 0.01,
                qmax = 0.99,
                by = -0.01
              ),
              "wrong sign in 'by' argument",
              fixed = TRUE
            )
            expect_error(
              weighted_dose_costs(
                "norm",
                params = c(MEAN = 1.85,
                           sd = 0.25),
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE
              ),
              "mismatched arguments for function 'qnorm'",
              fixed = TRUE
            )
            expect_error(
              weighted_dose_costs(
                "norm",
                params = list(mean = c(1.85, 1.87),
                           sd = 0.25),
                var_base = 0,
                dose_base = 0,
                dose_multiplier = 320,
                available_units = vialCost,
                "treatment",
                "fake",
                share_vials = FALSE
              ),
            "can only have one value for parameter 'mean'")
          })

test_that("cost_iv_compound_administration works",
          {
            exampleParams <-
              data.frame(
                compound = "X",
                param = c("cost_admin_first_hr", "cost_admin_addl_hr",
                          "iv_time_hr"),
                value = c(100, 20, 1.5)
              )
            expect_equal(cost_iv_compound_administration(exampleParams, "X"),
                         110)
            exampleParams <-
              data.frame(
                compound = "X",
                param = c("cost_admin_first_hr", "cost_admin_addl_hr",
                          "iv_time_hr"),
                value = c(100, 20, 0.75)
              )
            expect_equal(cost_iv_compound_administration(exampleParams, "X",
                                                         prorate_first = TRUE),
                         75)
            expect_equal(cost_iv_compound_administration(exampleParams, "X",
                                                         prorate_first = FALSE),
                         100)
            expect_error(
              cost_iv_compound_administration(exampleParams, "Y"),
              "Y not present in 'compound' column of data_table",
              fixed = TRUE
            )
            expect_error(
              cost_iv_compound_administration(exampleParams, "X",
                                              time_col = "iv_time_min"),
              "parameter iv_time_min not in data_table",
              fixed = TRUE
            )
          })

test_that("linear interpolation function works",
          {
            fn1 <- linear_interpolation(c(1, 10, 20), c(1, .5, .25))
            expect_equal(fn1(c(1, 5.5, 10, 11, 15, 20, 25)),
                         c(1, 0.75, 0.5, 0.475, 0.375, 0.25, 0.25))
            fn2 <-
              linear_interpolation(change_per_cycle = -0.02,
                                   max_change_cycles = 50)
            expect_equal(fn2(c(1, 20, 50, 51)),
                             c(0.98, 0.60, 0.00, 0.00))
            fn3 <-
              linear_interpolation(
                change_per_cycle = -0.02,
                max_change_cycles = 100,
                min_val = -0.8
              )
            expect_equal(fn3(
              c(1, 20, 50, 51, 69, 70, 71, 100, 150)),
              c(0.98, 0.60, 0.00,-0.02,-0.38,-0.40,-0.42,-0.80,-0.80)
            )
            fn4 <-
              linear_interpolation(
                age_start_change_per_cycle = 65,
                age_at_cycle_0 = 60,
                change_per_cycle = -0.05,
                cycles = 1:15,
                max_change_cycles = 5
              )
            expect_equal(fn4, c(rep(1, 6),
                                0.95, 0.9, 0.85, 0.8,
                                rep(0.75, 5)))
            expect_error(
              linear_interpolation(c(1, 10, 20), c(1, .5, .25, 0.1)),
              "marked_cycles and marked_values must have the same length"
            )
            expect_error(
              linear_interpolation(max_change_cycles = 5),
              "either marked_cycles and marked_values or change_per_cycle"
            )
          })
