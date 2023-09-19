#' Asymptotic variacne function
#'
#' Generates the asymptotic variance of the parameter estimators in the multi-dimensional SPDE model.
#' @param estimationMethod If \code{"OracleSigma"} is chosen, the funciton calculates the asymptotic variance for the oracle-volatility estimator. Therefore, provide the true parameters \code{alphaDah,sigma} of the model.
#' If \code{"SigmaAndKappa"} is chosen, please provide the true parameters \code{alphaDash, sigma_0_squared,indexset,SG}. Note that in this case, the function returns the approximated asymptotic variance, dependent on the provided spatial coordinates.
#' For additional information, see [SecondOrderSPDEMulti::estimateParametersMulti].
#' @param d  Spatial dimension of the SPDE model.
#'
#' @keywords asymptotic variance
#' @export
#' @seealso [SecondOrderSPDEMulti::SecondOrderSPDEMulti],[SecondOrderSPDEMulti::estimateParametersMulti].
#' @return a value/matrix containing the asymptotic variance of the provided estimation method.

#' @examples
#' library(data.table)
#' d <- 2
#' alphaDash <- 0.5
#' sigma_0_squared <- 1
#' indexset <- c(15,47,83)
#' SG <- do.call(CJ,replicate(d,seq(0,1,1/10),F))
#' av <- asymp_var("SigmaAndKappa",d,alphaDash=alphaDash,sigma_0_squared=sigma_0_squared,indexset = indexset,SG = SG)




asymp_var <- function(estimationMethod,d,alphaDash=NA,sigma=NA,sigma_0_squared=NA,SG=NA,indexset=NA){
  if(!(estimationMethod %in% c("OracleSigma","SigmaAndKappa","alphaDash","all")) ){
    stop("'estimationMethod' must be either 'OracleSigma', 'SigmaAndKappa', 'alphaDash', or 'all'.")
  }
  if (!require("pacman")) {
    install.packages("pacman")
    require(pacman)
  } else {
    require(pacman)
  }
  pacman::p_load(dplyr, pbmcapply,data.table,MASS, matlib)
  numCores <- detectCores() - 1

  if(estimationMethod == "OracleSigma"){
    res <- Upsilon(alphaDash,d)*sigma^4
    return(res)
  }
  if(estimationMethod == "SigmaAndKappa"){
    U <- Upsilon(alphaDash,d)
    I <- -diag(d+1)
    I[1,1] <- sigma_0_squared

    S <- SG[indexset,]
    M <- length(indexset)
    ones <- rep(1,M)
    X <- as.matrix(cbind(ones,S))
    Sigma <- crossprod(X)/M
    res <- I %*% inv(Sigma) %*% I

    return(res*U)
  }
}
