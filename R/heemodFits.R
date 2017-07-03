#' heemodFits: Survival Fits Organized for Heemod
#' 
#' @docType package
#' @name heemodFits
#'   
#'   
#' @importFrom heemod compute_surv
#'    
#' @importFrom dplyr filter_
#' @importFrom dplyr mutate_
#' @importFrom dplyr do_
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise_
#' @importFrom dplyr as.tbl
#' @importFrom dplyr data_frame
#' @importFrom dplyr as_data_frame
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom dplyr "%>%"
#' @importFrom dplyr desc
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate_if
#' @importFrom dplyr funs
#' 
#' @importFrom heemod discount
#' @importFrom heemod compute_surv
#' 
#' @importFrom plyr ldply
#' @importFrom plyr ddply
#'   
#' @importFrom lazyeval lazy
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval as.lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom lazyeval interp
#'   
#' @importFrom pryr standardise_call
#'   
#' @importFrom utils head
#' @importFrom utils modifyList
#' @importFrom utils globalVariables
#' @importFrom utils as.roman
#'   
#' @importFrom stats pnorm
#' @importFrom stats qbeta
#' @importFrom stats qbinom
#' @importFrom stats qgamma
#' @importFrom stats qlnorm
#' @importFrom stats qnorm
#' @importFrom stats terms
#' @importFrom stats setNames
#' @importFrom stats reorder
#' @importFrom stats na.omit
#' @importFrom stats update
#' @importFrom stats as.formula
#' @importFrom stats var
#' @importFrom stats coef
#' @importFrom stats model.matrix
#' @importFrom stats formula
#' @importFrom stats stepfun
#'   
#' @importFrom flexsurv flexsurvreg
#'    
#' @importFrom graphics par
#' 
#' @importFrom lpSolve lp
#'   
#' @importFrom mvnfast rmvn
#'   
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 "%+replace%"
#'   
#' @importFrom grDevices dev.off
#' @importFrom grDevices cairo_pdf
#' @importFrom grDevices png
#' 
#' @importFrom graphics plot
#'   
#' @importFrom memoise memoise
#' @importFrom memoise timeout
#'   
#' @importFrom rmarkdown render
#' 
#' @importFrom survival Surv
#' @importFrom survival coxph
#' 
#' @importFrom survminer ggsurvplot_combine
#' @importFrom survminer surv_fit
#' @importFrom survminer ggcoxdiagnostics
#'    
#' @importFrom tibble tibble
#' @importFrom tibble tibble_

#' @importFrom tidyr gather_
#'   
#' @importFrom tools file_ext
#'   
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @importFrom utils packageVersion
#' 
#' @importFrom XLConnect loadWorkbook
#' @importFrom XLConnect createSheet
#' @importFrom XLConnect writeWorksheet
#' @importFrom XLConnect saveWorkbook
#' 
#' @importFrom readxl read_excel

NULL

