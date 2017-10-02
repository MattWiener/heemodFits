#' Read the accepted file formats for tabular input
#'
#' Columns starting with '.comment' are ignored.
#' 
#' @param file_name File name.
#' 
#' @return A `data.frame`.
#'   
read_file <- function(file_name) {
  
  have_xls <- is_xls(file_name)
  have_xlsx <- is_xlsx(file_name)
  have_csv <- is_csv(file_name)
  
  if(! have_csv & ! have_xls & ! have_xlsx) {
    stop("file names must be for csv, xls, or xlsx")
  }
  
  if(have_csv) {
    tab <- utils::read.csv(
      file_name, header = TRUE,
      stringsAsFactors = FALSE,
      strip.white = TRUE,
      na.strings = c(".", "NA", "")
    ) 
  } else if(have_xls | have_xlsx) {
    if (! requireNamespace("readxl", quietly = TRUE)) {
      stop("readxl packaged needed to read Excel files")
      
    } else {
      tab <- as.data.frame(
        readxl::read_excel(file_name)
      )
    }
  }
  
  ## get rid of "comment" columns, if any
  tab <- tab[! grepl("^\\.comment", names(tab))]
  
  ## get rid of NA rows that may have come from blank rows
  ## in file
  tab <- filter_blanks(tab)
  tab
}

#' Remove Blank Rows From Table
#' 
#' Remove rows were all values are `NA`.
#' 
#' Some rows can be left blanks in the input table for 
#' readability, this function ensures those rows are 
#' removed.
#' 
#' @param x A `data.frame`.
#'   
#' @return A `data.frame` without blank rows.
#'   
#' @keywords internal
filter_blanks <- function(x) {
  x[! apply(is.na(x), 1, all), , drop = FALSE]
}

#' Check File Type
#' 
#' @param x A file name.
#' @return Whether the file is (respectively)
#'  csv, xlsx, or xls.
#' @rdname file-checkers
#'   
#' @keywords internal
is_csv <- function(x) {
  tolower(tools::file_ext(x)) == "csv"
}

#' @rdname file-checkers
is_xlsx <- function(x) {
  tolower(tools::file_ext(x)) == "xlsx"
}

#' @rdname file-checkers
is_xls <- function(x) {
  tolower(tools::file_ext(x)) == "xls"
}



#' Check Wholenumbers
#' 
#' @param x numeric.
#' @param tol the smallest positive floating-point number x 
#'   such that 1 + x != 1.
#'   
#' @return A logical scalar.
#'   
#' @keywords internal
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}


#' Returns "s" if x > 1
#'
#' @param x integer.
#'
#' @return `"s"` or `""`.
#'   
#' @keywords internal
plur <- function(x) {
  if (x > 1) "s" else ""
}
#' @rdname plur
plur_y <- function(x) {
  if (x > 1) "ies" else "y"
}


to_text_dots <- function(x, name = TRUE) {
  n <- names(x)
  ex <- if (is.atomic(x)) {
    format(x)
  } else {
    unlist(lapply(
      x,
      function(y) if (any(is.na(y))) NA else
        deparse(y$expr, width.cutoff = 500L)
    ))
  }
  
  if (name) {
    stopifnot(
      length(n) == length(ex)
    )
    paste(n, ex, sep = " = ")
  } else {
    ex
  }
}


appendEnv = function(e1, e2) {
       listE1 = ls(e1)
       listE2 = ls(e2)
       for(v in listE2) {
             if(v %in% listE1) stop(sprintf("Variable %s is in e1, too!", v))
             e1[[v]] = e2[[v]]
         }
   }
