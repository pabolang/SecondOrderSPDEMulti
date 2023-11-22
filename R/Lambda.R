#' Calculate Lambda
#'
#' Approximation of the constant Lambda which results by non-negligible covariances of squared increments of two temporal resolutions in the SPDE model.
#' @param n natural number for the approximation of the infinite sum. Default is 10000.
#' @param alphaDsh a real number in \code{(0,1)} controlling the roughness of the temporal marginal process.
#' @keywords Approximation of Lambda
#' @references Bossert, P. (2023), 'Parameter estimation for second-order SPDEs in multiple space dimensions'

#' @export
#' @seealso [SecondOrderSPDEMulti::SecondOrderSPDEMulti].
#' @return A real number of the approximated constant.

#' @examples
#' Lambda <- Lambda(alphaDash = 1/3)



Lambda <- function(alphaDash,cutoff=10000){
  res1 <- (2^alphaDash-1)^2-(2^alphaDash-2)^2-1
  l <- lapply(0:cutoff, function(r){
    (-(r+1)^alphaDash+2*(r+2)^alphaDash-(r+3)^alphaDash)*(-(r)^alphaDash+2*(r+1)^alphaDash-(r+2)^alphaDash)
  })
  res2 <- sum(unlist(l))
  return(res1+res2)
}
