#' Spatial grid filter index
#'
#' Filters the spatial grid by 'delta', where each entry has to be greater or equal \code{spatialDelta} and less or equal 1-\code{spatialDelta} and returns index of the original spatial grid data frame.
#' @param SG a data frame with the multi-dimensional spatial grid coordinates.
#' @param spatialDelta  filter variable. Each entry has to be greater or equal \code{spatialDelta} and less or equal 1-\code{spatialDelta}.
#' @keywords index for filtered spatial grid
#' @export
#' @seealso [SecondOrderSPDEMulti::SecondOrderSPDEMulti]
#' @return a vector containing the indexset of the spatial grid coordinates which satisfies the 'delta' condition.

#' @examples
#' library(data.table)
#' M <-10
#' d <- 3
#' SG <- do.call(CJ,replicate(d,seq(0,1,1/M),F))
#' delta <- 0.05
#'
#' get_filtered_SP_index(SG,delta)






get_filtered_SP_index <- function(SG,spatialDelta){

  d <- dim(SG)[2]
  delta <- spatialDelta
  l <- lapply(1:dim(SG)[1],function(i){
    if(sum(SG[i,]>= delta & SG[i,]<= 1-delta)==d){i}
  })
  return(unlist(l))
}
