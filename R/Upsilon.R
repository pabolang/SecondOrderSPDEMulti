#' Calculate Upsilon
#'
#' Approximation of the constant Upsilon which results by non-negligible covariances of squared increments in the SPDE model.
#' @param n natural number for the approximation of the infinite sum. Default is 10000.
#' @param d  natural number which denotes the spacial dimension. 'd' must be greater 1.
#' @param alphaDsh a real number in \code{(0,1)} controlling the roughness of the sample paths.
#' @keywords Approximation of Upsilon
#' @references PhD thesis Bossert, P.
#' @export
#' @seealso [SecondOrderSPDEMulti::SecondOrderSPDEMulti].
#' @return A real number of the approximated constant.

#' @examples
#' Upsilon <- Upsilon(d=2, alphaDash = 1/3)


Upsilon <- function(alphaDash,d,n=10000){
  l <- lapply(0:n, function(r){
    (2*(r+1)^alphaDash-r^alphaDash-(r+2)^alphaDash)^2
  })
  erg <- (sum(unlist(l))+2)
  return(erg)
}





