
pad_matrix <- function(orig_matrix, pad_to){
  new_matrix <- matrix(NA, nrow = nrow(orig_matrix), ncol = pad_to)
  new_matrix[, 1:ncol(orig_matrix)] <- orig_matrix
  dimnames(new_matrix) <- list(rownames(orig_matrix), NULL)
  new_matrix
}

prepare_vcov <- function(fit_list){
  km_pos <- match("km", names(fit_list))
  ## drop Kaplan-Meier, if present
  km_pos <- km_pos[!is.na(km_pos)]
  if(length(km_pos) > 0)
    fit_list <- fit_list[-km_pos]
  
  max.n.pars <- max(sapply(fit_list, function(x){x$npars}))
  
  vcov_matrices <- lapply(fit_list, stats::vcov)
  vcov_matrices <- lapply(vcov_matrices, pad_matrix, 
                          pad_to = max.n.pars)
  prefixes <- names(fit_list)
  for(i in 1:length(vcov_matrices))
    rownames(vcov_matrices[[i]]) <- 
    paste(prefixes[i], rownames(vcov_matrices[[i]]), sep = "_")
  
  do.call("rbind", vcov_matrices)
}

prepare_fit_info <- function(fit_list){
  ## drop Kaplan-Meier, if present
  km_pos <- match("km", names(fit_list))
  km_pos <- km_pos[!is.na(km_pos)]
  if(length(km_pos) > 0)
    fit_list <- fit_list[-km_pos]
  par.est <- lapply(fit_list,
                    function(x) {
                      res <- x$res.t[, 1]
                      names(res) <-rownames(x$res.t)
                      res
                    })
  AIC <- sapply(fit_list, get_component, comp = "AIC") 
  BIC <- sapply(fit_list, get_component, comp = "BIC") 
  vcov <- prepare_vcov(fit_list)
  list(par.est = par.est, AIC = AIC, BIC = BIC, vcov = vcov)
}

#' Write fit information to an Excel workbook.
#'
#' @param fit_tibble a tibble of fits 
#' @param wb the workbook object, or a character string with its file path
#' @param skip_at_start how many lines to skip at the top of each page
#' @param skip_between how many lines to skip between output sections
#' @param alignment either "horizontal" or "vertical"
#' @param B_ci number of simulations to use for confidence intervals
#'   for `flexsurvreg` objects; becomes the `B` of 
#'   [flexsurv::summary.flexsurvreg()]
#' @details Called for the side effect of creating Excel output.
#' 
#' `alignment = 'vertical'` means sections are written one after
#'   the other down the page.   If `alignment = 'horizontal'`, then
#'   some parts will be written next to each other.
#' @return data suitable for plotting, invisibly.
#' @export
#'
write_fits_to_excel_from_tibble <- 
  function(fit_tibble, wb, skip_at_start = 3, skip_between = 1,
           alignment = c("horizontal", "vertical"),
           B_ci = 100){
    requireNamespace("XLConnect")
    requireNamespace("flexsurv")

    ## if necessary, get flexsurvreg fits out of surv_X objects
    if(is.character(wb)) 
      wb <- XLConnect::loadWorkbook(wb, create = TRUE)
    alignment <- match.arg(alignment)
    stopifnot(inherits(wb, "workbook"))
    stopifnot(identical(names(fit_tibble),
                c("type", "treatment", "set_name",
                  "dist", "fit", "set_def", "time_subtract"))
              )
    XLConnect::createSheet(wb, "OPCPem")

    # fit_tibble_multi_treatment <- 
    #   dplyr::filter(fit_tibble, treatment == "all")
    # fit_tibble <-
    #   dplyr::filter(fit_tibble, treatment != "all")
    # 
    fit_tibble_nest <- 
      dplyr::filter_(fit_tibble, ~ dist != "km")
    fit_tibble_nest$fit <- lapply(fit_tibble_nest$fit, extract_fits)
    fit_tibble_nest <- 
      fit_tibble_nest %>%
        dplyr::group_by_(~ type, ~ treatment, ~ set_name) %>%
          tidyr::nest()
    
    
    flexspace <- loadNamespace("flexsurv")
    ## need to know how many rows of output we're going to write
    num_rows <- 
      sapply(fit_tibble_nest$data, function(x){
        1 + sum(sapply(x$fit, function(y){1 + length(coef(y))}))
      })
    
    if(alignment == "vertical"){
    rows_per <- 
      2 * (num_rows + skip_between) + 
        sapply(fit_tibble_nest$data, nrow) + 1
    }
    if(alignment == "horizontal"){
      rows_per <- 
        num_rows + skip_between + sapply(fit_tibble_nest$data, nrow) + 1
    }
    fit_tibble_nest$start_row <- 
      skip_at_start + 1 + cumsum(rows_per) - rows_per[1]

    
    fit_tibble_nest <-
      fit_tibble_nest %>%
        dplyr::group_by_(~ type, ~ treatment, ~ set_name) %>%
          dplyr::do_(~send_info_to_workbook(., wb = wb, 
                                          skip_between = skip_between,
                                          alignment = alignment))
    use_pieces_str <- 
      "c('time', 'n.risk', 'n.event', 'surv', 'std.err', 'lower', 'upper')"
    do_str <- 
      paste("data.frame(summary(object = .$fit[[1]], summary_type = 'standard')[",
            use_pieces_str,
            "])",
            sep = ""
            )
    base_km_summaries <- 
      fit_tibble %>% 
      dplyr::filter_(~ dist == "km") %>%
      dplyr::group_by_(~ type, ~ treatment, ~ set_name) %>%
      dplyr::do_(do_str) %>%
      dplyr::ungroup()
    
    
    plot_data <- prepare_plot_data_from_fit_tibble(fit_tibble, 
                                                   B_ci = B_ci)
    plot_data_km <- plot_data %>% 
      dplyr::filter_(~ dist == "km", ~fn == "survival")
  
    XLConnect::createSheet(wb, "km")
    XLConnect::writeWorksheet(wb,
                              base_km_summaries,
                              sheet = "km",
                              startRow = skip_at_start,
                              startCol = 1)
    
    XLConnect::createSheet(wb, "km_plot")
    XLConnect::writeWorksheet(wb, 
                   plot_data_km[, c("type", "treatment", "set_name", "dist",
                                         "time", "est", "lcl", "ucl", "fn")],
                   sheet= "km_plot", 
                   startRow = skip_at_start,
                   startCol = 1)
          
    XLConnect::saveWorkbook(wb)
    invisible(plot_data)
  }


send_info_to_workbook <- 
  function(fit_tib, wb, skip_between = 1,
           alignment = c("horizontal", "vertical")){

  alignment <- match.arg(alignment)
  
  fit_list <- fit_tib[[1, "data"]][["fit"]]
  names(fit_list) <- fit_tib[[1, "data"]][["dist"]]
  start_row <- fit_tib[[1, "start_row"]]
  dist_names <- names(fit_list)
  fit_info <- prepare_fit_info(fit_list)
  
  if(alignment == "vertical"){
    start_row_2 <- start_row + length(dist_names) + skip_between + 1
    start_row_3 <- 
      start_row_2 + length(unlist(fit_info$par.est)) + skip_between
    start_col_1 <- start_col_2 <- start_col_3 <- 1
  }
  if(alignment == "horizontal"){
    start_row_2 <- start_row + length(dist_names) + skip_between + 1
    start_row_3 <- start_row_2
    start_col_1 <- 1
    start_col_2 <- 1
    start_col_3 <- 8
    
  }
  
  stat_output <- data.frame(fit_tib[[1, "type"]],
                            fit_tib[[1, "treatment"]],
                            fit_tib[[1, "set_name"]],
                            dist = dist_names, 
                           AIC = fit_info$AIC,
                           BIC = fit_info$BIC)
  names(stat_output)[1:3] <- ""
    XLConnect::writeWorksheet(wb, stat_output,
                 sheet = "OPCPem",
                  startRow = start_row, startCol = 1, header=TRUE)
 #   new_start_row <- start_row + length(dist_names) + skip_between + 1
    
    write_pars <- data.frame(
      fit_tib[[1, "type"]],
      fit_tib[[1, "treatment"]],
      fit_tib[[1, "set_name"]],
      unlist(lapply(fit_info$par.est, names)),
      rep(names(fit_info$par.est),
          sapply(fit_info$par.est, length)),
      unlist(fit_info$par.est))
    XLConnect::writeWorksheet(wb, write_pars,
                   sheet = "OPCPem",
                   startRow = start_row_2, 
                   startCol = start_col_2, 
                   header = FALSE) 
    use_names <- rownames(fit_info$vcov)
    write_vcov <- fit_info$vcov
    if(alignment == "vertical"){
    write_vcov <- data.frame(fit_tib[[1, "type"]],
                            fit_tib[[1, "treatment"]],
                            fit_tib[[1, "set_name"]],
                            do.call("rbind", strsplit(use_names, "_")),
                            write_vcov)
    write_vcov[, 4:5] <- write_vcov[, 5:4]
    }
  XLConnect::writeWorksheet(wb, write_vcov,
                 sheet = "OPCPem",
                 startRow = start_row_3,
                 startCol = start_col_3, header = FALSE)
  data.frame(numeric(0))
}


#' Prepare fit data for plotting
#'
#' @param fit_tib a tibble containing fits
#' @param B_ci number of simulations to use for confidence intervals
#'   for `flexsurvreg` objects; becomes the `B` of 
#'   [flexsurv::summary.flexsurvreg()]
#' @return a tibble with the necessary data
#' @export
#'
prepare_plot_data_from_fit_tibble <-
  function(fit_tib, B_ci = 100){
    if(!inherits(fit_tib, c("tbl_df", "tbl", "data.frame")))
      stop("fit_tib should be a tibble")
    required_names <- c("type", "treatment", "set_name", "dist", "fit")
    missing_names <- setdiff(required_names, names(fit_tib))
    if(length(missing_names) > 0)
      stop("fit_tib is missing column", 
           plur(length(missing_names)),
           " with name",
           plur(length(missing_names)),
           ": ", 
           paste(missing_names, collapse = ", ")
           )
    
    do_str_surv <- paste("heemodFits:::summary_helper(., type = 'survival', tidy = TRUE, B = ", 
                    B_ci,
                    ")", 
                    sep = "")
    survival_summaries <- 
      fit_tib %>% 
        dplyr::group_by_(~ type, ~ treatment, ~ set_name, ~ dist) %>%
      ##summary_helper(type = 'survival', tidy = TRUE, B = B_ci) %>%
          dplyr::do_(do_str_surv) %>%
            dplyr::ungroup()
    survival_summaries$fn <- "survival"
    do_str_cumhaz <- paste("heemodFits:::summary_helper(., type = 'cumhaz', tidy = TRUE, B = ", 
                    B_ci,
                    ")",
                    sep = "")
    cumhaz_summaries <- 
      fit_tib %>% 
      dplyr::group_by_(~ type, ~ treatment, ~ set_name, ~ dist) %>%
      ##summary_helper(type = 'cumhaz', tidy = TRUE, B = B_ci) %>% 
      dplyr::do_(do_str_cumhaz) %>%
      dplyr::ungroup()
    cumhaz_summaries$fn <- "cumulative hazard"
    rbind(survival_summaries, cumhaz_summaries)
  }


summary_helper <- function(fit_holder, type, ...){
  fit <- fit_holder$fit[[1]]
  stopifnot(inherits(fit, c("flexsurvreg", "survfit", "surv_shift",
                            "surv_projection")))
    res1 <- summary(fit, type = type, 
                    summary_type = "standard", ...)
    all_times <- sort(unique(c(res1[["time"]],
                               seq(from = 1, 
                                   to = max(res1[["time"]]),
                                   by = 1))))
    res1 <- summary(fit, type = type, t = all_times, 
                    summary_type = "standard", ...)

  if(inherits(res1, "summary.survfit")){
    res1 <- data.frame(res1[c("time", "surv", "lower", "upper")])
    names(res1) <- c("time", "est", "lcl", "ucl")
    if(type == "cumhaz")
      res1[, -1] <- -log(res1[, -1])
    res1
  }
  res1
}



#' Plot fit data
#'
#' @param data_to_plot a data frame from 
#'   [prepare_plot_data_from_fit_tibble()]
#' @param plot_type `survival` or `cumulative hazard`
#' @param plot_groups which groups should be included
#'   in the plot.   If NULL, all distributions in the data
#'   are included.
#' @param groups column name to use for grouping.
#' @param scale_time times are multiplied by this.  So if your
#'   fit was done on a time scale of days, and you want to plot
#'   by weeks, set this to 1/7.
#' @param time_label the `x` axis label
#' @param max_scaled_time maximum time, after scaling by `scale_time`
#' @param title for the plot
#' @param x_axis_gap distance between breaks on the x axis
#' @param legend_loc location for legend - "top", "bottom", "left", or "right"
#' @param logy should the `y` axis be on the logarithmic scale?  By default,
#'   it is TRUE for survival curves and FALSE for cumulative hazard.
#' @param km_width width for the Kaplan-Meier curve line;
#'   approximately 0.5 seems to be the standard width for lines
#' @param label_size relative multiplier for label size.
#' @param silence_warnings suppresses certain warnings
#' @return a `ggplot2` plot object
#' @export
#'
#' @examples
plot_fit_data <- function(data_to_plot, 
                          plot_type = c("survival", "cumulative hazard",
                                        "log cumulative hazard"),
                          plot_groups = NULL, 
                          groups = "dist",
                          scale_time = 1,
                          time_label = "time",
                          max_scaled_time = Inf,
                          title = NULL,
                          x_axis_gap,
                          legend_loc = "right",
                          logy = NULL, 
                          km_width = 1.25,
                          label_size = 1.25,
                          silence_warnings = FALSE){
  if(length(unique(data_to_plot$type)) > 1)
    warning("more than one type in plotting data - could cause problems")
  if(length(unique(data_to_plot$treatment)) > 1 & !silence_warnings)
    warning("more than one treatment in plotting data - could cause problems")
  if(length(unique(data_to_plot$set_name)) > 1)
    warning("more than one subset name in plotting data - could cause problems")
  if(!("km" %in% unique(data_to_plot$dist)))
    stop("no 'km' in dist column - but Kaplan-Meier curve is required")
  plot_type <- match.arg(plot_type)
  default_logy <- data.frame(
    plot_type = c("survival", "cumulative hazard", "log cumulative hazard"),
    logy = c(FALSE, FALSE, TRUE)
  )
  if(is.null(logy))
    logy <- default_logy$logy[match(plot_type, default_logy$plot_type)]
  
  stopifnot(legend_loc %in% c("right", "left", "bottom", "top"))
  data_to_plot <- dplyr::filter_(data_to_plot, 
                                 lazyeval::interp(~fn == var, 
                                                  var = plot_type))
  data_to_plot$time <- data_to_plot$time * scale_time
  data_to_plot <- dplyr::filter_(data_to_plot, ~ time <= max_scaled_time)

  unique_groups <- unique(data_to_plot[[groups]])
  use_colors <- 
    grDevices::hcl(h = seq(15, 375, length = length(unique_groups) + 1),
                   l = 65, c = 100)[1:length(unique_groups)]

  use_colors[unique_groups == "km"] <- "black"

  unique_labels <- unique_groups
  unique_labels[unique_groups == "km"] <- "Kaplan-Meier"
  line_widths <- rep(1, length(unique_groups))
  line_widths[unique_groups == "km"] <- 1.5
  if(!is.null(plot_groups)){
    not_present <- setdiff(plot_groups, data_to_plot[[groups]])
    if(length(not_present) > 0)
      warning("element",
              plur(length(not_present)),
              " of column ",
              groups, 
              " not in the data to plot: ",
              paste(not_present, collapse = ", "),
              ".\n",
              "Available elements: ",
              paste(unique(data_to_plot[[groups]]), collapse = ", ")
      )
    filter_str <- paste(groups, "%in% c(", 
                        paste("'",plot_groups, "'", sep = "", collapse = ", "), 
                        ")"
                        )
  data_to_plot <- dplyr::filter_(data_to_plot, 
                                 filter_str)
  }
  seed_frame <- 
    data.frame(unique_groups,
               lwd = line_widths,
               stringsAsFactors = FALSE)
  names(seed_frame)[1] <- groups
  data_to_plot <- 
    dplyr::left_join(data_to_plot,
                     seed_frame,
                     by = groups)
  names(use_colors) <- unique_groups
  names(unique_labels) <- unique_groups
  ## I'm sure there's a better way to do this with aes_string
  ##   or something, but haven't figured it out
  ## names(data_to_plot)[names(data_to_plot) == groups] <- "groups"
  res <- 
    ggplot2::ggplot(data_to_plot,
                    ggplot2::aes_string(x = "time", y = "est", color = groups)) + 
      ggplot2::geom_line(data = data_to_plot) +
      ggplot2::scale_color_manual(labels = unique_labels,
                                  values = use_colors) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.5)))
  if(groups == "dist")
    res <- res + 
    ggplot2::geom_step(ggplot2::aes_string(x = "time", y = "est"), 
                       data = dplyr::filter_(data_to_plot, 
                                             ~ dist == "km"), 
                       col = "black", 
                       lwd = km_width)
  if(logy) res <- res + ggplot2::scale_y_log10()
  if(plot_type == "log cumulative hazard") 
    res <- res + ggplot2::scale_x_log10()
  if(!is.infinite(max_scaled_time) & !missing(x_axis_gap)){
    breaks <- seq(from = 0, to = max_scaled_time, by = x_axis_gap)
    res <- res + ggplot2::scale_x_continuous(breaks = breaks)
  }
  res <- res + ggplot2::labs(y = plot_type, x = time_label, title = title)
  res <- 
    res + ggplot2::theme(axis.title = ggplot2::element_text(size = ggplot2::rel(label_size)))
  res <- res + ggplot2::theme(legend.position = legend_loc)
  res
}

get_component <- function(obj, comp){
  if(comp %in% names(obj)) return(obj[[comp]])
  else{
    if("dist" %in% names(obj)) return(get_component(obj[["dist"]], comp))
    else return(NULL)
  }
}

extract_fits <- function(x) {
  if (inherits(x, c("flexsurvreg", "survfit")))
    x
  else{
    if ("dist" %in% names(x) && 
        inherits(x$dist, c("flexsurvreg", "survfit")))
      x$dist
    else
      stop(
        "unrecognized input; not a flexsurvreg or survfit object and ",
        "doesn't contain a flexsurvreg or survfit object as 'dist'"
      )
  }
}

#' Assemble data from a fit tibble and plot
#'
#' @param fit_tibble the output of [survival_fits_from_tabular()].
#' @param treatment the treatment for which we want to plot.
#' @param set_name subset name for all but the Kaplan-Meier curve.
#' @param type PFS or OS
#' @param set_for_km subset name for the Kaplan-Meier curve.
#' @param B_ci number of simulations to use for confidence intervals
#'   for `flexsurvreg` objects; becomes the `B` of 
#'   [flexsurv::summary.flexsurvreg()]
#' @param scale_time times are multiplied by this.  So if your
#'   fit was done on a time scale of days, and you want to plot
#'   by weeks, set this to 1/7
#' @param ... additional parameters to pass to [plot_fit_data()]
#'
#' @return a list of two `ggplot` plots - one for survival and one
#'   for cumulative hazard
#' @export
#'
#' @examples
plot_fit_tibble <-
  function(fit_tibble, treatment, set_name, type, 
           set_for_km = "all", B_ci = 100, 
           scale_time = 1, ...){
    check_plot_tibble_input(fit_tibble, treatment, set_name, type)
    if(length(treatment) > 1 | length(set_name) > 1 | length(type) > 1)
      stop("can only enter a single treatment, set_name, and type")

    partial1 <- 
      dplyr::filter_(fit_tibble, 
                     lazyeval::interp(~type == var, var = type),
                     lazyeval::interp(~treatment == var, var = treatment)
                     )
    not_km <- 
      partial1 %>%
      dplyr::filter_(lazyeval::interp(~set_name == var, var = set_name),
                      ~ dist != "km")
    time_subtract <- unique(not_km$time_subtract)
    if(length(time_subtract) > 1)
      stop("more than one time_subtract specified - doesn't make sense")
    for_km <- partial1 %>%
      dplyr::filter_(~dist == "km", 
                     lazyeval::interp(~set_name == var, var = set_for_km))
    if(nrow(for_km) == 0)
      stop("no Kaplan-Meier curve data using given set_for_km: '",
           set_for_km,
           "'"
           )
    if(unique(not_km$set_name) != set_for_km){
    not_km$fit <- 
      lapply(1:nrow(not_km), 
             function(i){
               heemod::join(for_km$fit[[1]],
                            not_km$fit[[i]],
                            at = time_subtract) 
             }
      )
    }
    not_km_plot_data <- lapply(1:nrow(not_km),
          function(i){
            x <- not_km[i,]
            prepare_plot_data_from_fit_tibble(x, B_ci = B_ci) %>%
              dplyr::filter_(~time > x$time_subtract)
          })
    not_km_plot_data_2 <- dplyr::bind_rows(not_km_plot_data)
    km_plot_data <- prepare_plot_data_from_fit_tibble(for_km)
    fit_data <- dplyr::bind_rows(km_plot_data, not_km_plot_data_2)
    fit_data$set_name <- "for_plot"
    res_surv <- plot_fit_data(fit_data, plot_type = "survival", 
                              scale_time = scale_time, ...)
    res_cumhaz <- plot_fit_data(fit_data, 
                                plot_type = "cumulative hazard", 
                                scale_time = scale_time, 
                                ...)
    
    use_title <- paste(type, " for subset '", 
          set_name, "' of treatment '", 
          treatment, "'", sep = "")
    
    res_surv <- res_surv + ggplot2::ggtitle(use_title)
    res_cumhaz <- res_cumhaz + ggplot2::ggtitle(use_title)
    
    if(!is.na(time_subtract) && time_subtract != 0){
      res_surv <- 
        res_surv + 
        ggplot2::geom_vline(xintercept = time_subtract * scale_time, 
                            lty = 2) 
      res_cumhaz <- 
        res_cumhaz + 
        ggplot2::geom_vline(xintercept = time_subtract * scale_time, 
                            lty = 2)
    }
    list(survival = res_surv, cumhaz = res_cumhaz)
  }


#' Assemble data from a fit tibble and plot a log cumulative hazard plot
#'
#' @param fit_tibble the output of [survival_fits_from_tabular()].
#' @param treatments the treatments to include; if left `NULL`,
#'   all treatments will be included.
#' @param set_name subset name
#' @param type PFS or OS
#' @param ... additional parameters to pass to [plot_fit_data()]
#'
#' @return a `ggplot` plot
#' @export
#'

plot_cloglog_fit_tibble <- 
  function(fit_tibble, treatments = NULL, set_name, type, ...){
    check_plot_tibble_input(fit_tibble, treatments, set_name, type)
    ## pull out just Kaplan-Meier fits,
    ##   with only the type we want (PFS or OS),
    ##   and with whatever treatments we want
    loadNamespace("survminer")
    partial1 <- 
      dplyr::filter_(fit_tibble,
                     ~ dist == "km",
                     lazyeval::interp(~type == var, var = type),
                     lazyeval::interp(~set_name %in% var, var = set_name)
      )
    
    if(!is.null(treatments))
      partial1 <- partial1 %>% 
                    dplyr::filter_( lazyeval::interp(~treatment %in% var, 
                                                     var = treatments))
    
     km_cloglog_plot_data <- 
       prepare_plot_data_from_fit_tibble(partial1)  %>%
        dplyr::filter_(~fn == "cumulative hazard") %>%
         dplyr::mutate(fn = "log cumulative hazard")

    res_cloglog <- 
      plot_fit_data(km_cloglog_plot_data, 
                    plot_type = "log cumulative hazard",
                    groups = "treatment",
                    time_label = "time (log scale)",
                    silence_warnings = TRUE,
                    ...)
    res_cloglog <- 
      res_cloglog + ggplot2::geom_hline(yintercept = 1,
                                        lty = "dashed") + 
      ggplot2::geom_line(size = 1.25)
    use_title <- paste(type, " for subset '", set_name, "'", sep = "")
    res_cloglog <- res_cloglog + ggplot2::ggtitle(use_title)
   res_km <-
      survminer::ggsurvplot_combine(lapply(partial1$fit, extract_fits), legend.title = "",
                                    legend.labs = partial1$treatment,
                                    conf.int = TRUE) + 
        ggplot2::labs(title = use_title, x = "time")
    res_km$plot$theme <- res_cloglog$theme
    
    ## we're also going to get the Cox proportional hazards
    ##  information.   This is a little bit of a hack, but
    ##  lets us put together the information we need out of
    ##  what we already have.  (It is currently hard to do this
    ##  in the main fitting routine in survival_fit.R because
    ##  of an earlier decision to have subsets defined by treatment,
    ##  which led to fitting treatment by treatment.   This is 
    ##  a convenient place where we have the data together across
    ##  treatments.   May want to fix this someday.)
    
    data_for_cox <- 
      dplyr::bind_rows(lapply(partial1$fit, 
                              function(x){y <- extract_data(x)
                                          names(y) <- c("time", "status", "treatment")
                                          y
                                          }
                              )
                       )
    cox_fit <- 
      survival::coxph(survival::Surv(time, status) ~ treatment,
                    data = data_for_cox)
    cox_schoenfeld_resid_plot <-
      survminer::ggcoxdiagnostics(cox_fit,
                                  type = "schoenfeld",
                                  ox.scale = "time",
                                  hline.col = "black",
                                  hline.size = 0.5)
    
    list(cloglog = res_cloglog, km = res_km, 
         cox_fit = cox_fit, cox_resid = cox_schoenfeld_resid_plot)
      
  }


check_plot_tibble_input <- 
  function(fit_tibble, treatment, set_name, type, ...){
    required_cols <- c("treatment", "type", "set_name", "fit") 
    present_cols <- required_cols %in% names(fit_tibble)
    if(!all(present_cols))
      stop("input fit_tibble is missing required column",
           plur(sum(!present_cols)),
           ": ",
           paste(required_cols[!present_cols], collapse = ", ")
      )
    if(!all(treatment %in% fit_tibble$treatment)){
      missing_treatment <- setdiff(treatment, fit_tibble$treatment)
      stop("treatment",
            plur(length(missing_treatment)),
           " ",
           paste(missing_treatment, collapse = ", "),
           " not present in entered fit_tibble.\n",
           "Available treatments: ",
           paste(unique(fit_tibble$treatment), collapse = ", ")
      )
    }
    if(!all(type %in% fit_tibble$type)){
      missing_type <- setdiff(type, fit_tibble$type)
      stop("type",
           plur(length(missing_type)),
           " ",
           paste(missing_type, collapse = ", "),
           " not present in entered fit_tibble.\n",
           "Available types: ",
           paste(unique(fit_tibble$type), collapse = ", ")
      )
    }
    if(!all(set_name %in% fit_tibble$set_name)){
      missing_set <- setdiff(set_name, fit_tibble$set_name)
      stop("set_name",
           plur(length(missing_set)),
           " ",
           paste(missing_set, collapse = ", "),
           " not present in entered fit tibble.\n",
           "Available set_names: ",
           paste(unique(fit_tibble$set_name), collapse = ", ")
      )
    }
  }

extract_data <- 
  function(x) {
    if (inherits(x, "survfit"))
      x$call$data
    else{
      if ("dist" %in% names(x) && 
          inherits(x$dist,  "survfit"))
        x$dist$call$data
      else
        stop(
          "unrecognized input; not survfit object and ",
          "doesn't contain a urvfit object as 'dist'"
        )
    }
  }
