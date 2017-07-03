context("testing survival fit output")

test_that("pad_matrix works",
          {
            orig_matrix <- matrix(1:6, nrow = 3, nc = 2)
            res_matrix <- cbind(orig_matrix, 
                                matrix(NA, nrow = 3, nc = 2))
            dimnames(res_matrix) <- list(NULL, NULL)
            expect_identical(pad_matrix(orig_matrix, 4),
                             res_matrix)
            dimnames(orig_matrix) <- dimnames(res_matrix) <- 
              list(letters[1:3], NULL)
            expect_identical(pad_matrix(orig_matrix, 4),
                             res_matrix)
          }
          )

test_that("prepare plot data from fit tibble",
          {
            location <- system.file("surv", package = "heemodFits")
            ok_surv_info <- 
              heemodFits:::read_file(system.file("surv/survival_info.csv", 
                                             package = "heemodFits"))
            full_these_fits <- 
              survival_from_data(location = location,
                                          survival_specs = ok_surv_info,
                                          dists = c("exp", "weibull"),
                                          save_fits = FALSE) 
            these_fits <- dplyr::filter_(full_these_fits,
                                         ~treatment != "all")
            aa <- prepare_plot_data_from_fit_tibble(these_fits[1:3, ], B_ci = 0)

            expect_equal(aa$est[seq(from = 50, to = 2000, by = 50)],
                         c(0.89900388, 0.80820798, 0.72658211, 0.65320014, 
                           0.58722946, 0.52792156, 0.47460354, 0.84, 0.66, 
                           0.6, 0.58, 0.56, 0.54, 0.84328836, 0.75056722, 
                           0.68373791, 0.63030209, 0.58551426, 0.54691921,
                          0.51303849, 0.08304498, 0.18951291, 0.29598084,
                            0.40244876, 0.50891669, 0.61538462, 0.72185254,
                            0.15082289, 0.32850407, 0.51082562, 0.54472718,
                            0.5798185, 0.61618614, 0.13827727, 0.26400002,
                            0.36090194, 0.44441258, 0.51957942, 0.58885222,
                            0.65365155)
                         )
            expect_equal(aa$dist,
                         rep(rep(c("exp", "km","weibull"), 
                                 c(350, 311, 350)),
                             2)
            )
            expect_equal(unique(aa$fn[1:1011]), "survival")
            expect_equal(unique(aa$fn[1011 + 1:1011]), "cumulative hazard")
            expect_error(prepare_plot_data_from_fit_tibble(these_fits$fit, B_ci = 0),
                         "fit_tib should be a tibble")
            expect_error(prepare_plot_data_from_fit_tibble(these_fits[, -1]),
                         "fit_tib is missing column with name: type",
                         fixed = TRUE)
            
          }
          )

test_that("get_component works",
          {
            expect_identical(
              get_component(list(a = 1, b = 2), "a"),
              1
            )
            expect_identical(
              get_component(list(a = 1, b = 2), "g"),
              NULL
            )
            expect_identical(
              get_component(list(a = 1, dist = list(b = 2, d = 3)),
                            "d"),
              3
            )
            expect_identical(
              get_component(list(a = 1, dist = list(b = 2, d = 3)),
                            "m"),
              NULL
            )
          })

test_that("extract_fits works",
          {
            x <- 1:5
            expect_error(extract_fits(x),
                         "unrecognized input"
                         )
            fake_flexsurvreg <- 1:5
            class(fake_flexsurvreg) <- "flexsurvreg"
            expect_identical(
              extract_fits(fake_flexsurvreg),
              fake_flexsurvreg
              )
            expect_identical(
              extract_fits(list(a = 1, dist = fake_flexsurvreg)),
              fake_flexsurvreg
              )
              
          })

test_that("fit_plot_tibble catches errors in input",
          {
            fake_fit_tib <- read_file(system.file("surv",
                                                  "fake_fit_tib.csv", 
                                                  package = "heemodFits"))
            expect_error(
              plot_fit_tibble(fake_fit_tib,
                              treatment = "wrong_treatment",
                              set_name = "all", type = "PFS"),
              "treatment wrong_treatment not present in entered fit_tibble",
              fixed = TRUE
            )
            expect_error(
              plot_fit_tibble(fake_fit_tib,
                              treatment = "A",
                              set_name = "not_a_set", type = "PFS"),
              "set_name not_a_set not present in entered fit tibble",
              fixed = TRUE
            )
            expect_error(
              plot_fit_tibble(fake_fit_tib,
                              treatment = "A",
                              set_name = "all", type = "PZZ"),
              "type PZZ not present in entered fit_tibble",
              fixed = TRUE
            )
            expect_error(
              plot_fit_tibble(fake_fit_tib[, -c(2,3)],
                              treatment = "A",
                              set_name = "all", type = "PFS"),
              "missing required columns: treatment, set_name",
              fixed = TRUE
            )
            expect_error(plot_fit_tibble(fake_fit_tib,
                                         treatment = c("A", "B"),
                                         set_name = "all", type = "PFS"),
                         "single treatment, set_name, and type",
                         fixed = TRUE
            )
            expect_error(plot_fit_tibble(fake_fit_tib,
                                         treatment = "A",
                                         set_name = "all", 
                                         type = c("PFS", "OS")),
                         "single treatment, set_name, and type",
                         fixed = TRUE
            )
            expect_error(plot_fit_tibble(fake_fit_tib,
                                         treatment = "A",
                                         set_name = "all", 
                                         set_for_km = "wrong",
                                         type = "PFS"),
                         "no Kaplan-Meier curve data using given set_for_km: 'wrong'",
                         fixed = TRUE
            )
          }
          )