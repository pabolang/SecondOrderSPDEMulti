#' Spatial grid filter
#'
#' Filters the spatial grid by 'delta', where each entry has to be greater or equal \code{spatialDelta} and less or equal 1-\code{spatialDelta}
#' @param SG a data frame with the multi-dimensional spatial grid coordinates.
#' @param spatialDelta  filter variable. Each entry has to  be greater or equal \code{spatialDelta} and less or equal 1-\code{spatialDelta}
#' @keywords Filter spatial grid
#' @export
#' @seealso [SecondOrderSPDEMulti::SecondOrderSPDEMulti]
#' @return a data frame containing the filtered spatial grid.

#' @examples
#' library(data.table)
#' M <-10
#' d <- 3
#' SG <- do.call(CJ,replicate(d,seq(0,1,1/M),F))
#' delta <- 0.05
#'
#' get_filtered_SP(SG,delta)


get_filtered_SP <- function(SG,spatialDelta){
  if (!require("pacman")) {
    install.packages("pacman")
    require(pacman)
  } else {
    require(pacman)
  }
  pacman::p_load(dplyr)
  delta <- spatialDelta
  y <- filter_all(SG,all_vars(.>=delta & .<=1-delta))
  return(y)
}
