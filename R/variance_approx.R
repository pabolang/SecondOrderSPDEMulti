#' Calculates the approximated variance of the replacement random normals
#'
#' Calculates the approximated variance of the replacement random normals within the replacement method. Of use, when considering a Monte carlo sudy. For details, see references or [SecondOrderSPDEMulti::simulateSPDEmodelMulti].
#' @param d a natural number, which denotes the spacial dimension. 'd' must be greater 1.
#' @param theta0 a real number, which controls the drift of the solution field.
#' @param nu a real vector of dimension 'd', which controls the curvature of the solution field on each axis, respectively.
#' @param eta a real number greater than zero, which reduces the noise level of parameter \code{sigma} and the curvature effect of parameter \code{nu}.
#' @param sigma a real number greater than zero, which controls the overall noise level of the solution field.
#' @param alphaDash a real number in \code{(0,1)}, controlling the roughness of the temporal marginal process.
#' @param numberOfSpatialPoints number of equidistant spatial points M on each axis, where we consider the domain \code{[0,1]}.
#' @param L a natural number indicating the replacement bound 'LM', dependent on multiples of the d-dimensional vector with entries \code{M}. The default is \code{L=10}.
#' @param K a natural number greater 'L' denoting the cut-off if the variance calculation for the replacement random variables (only of use when method is 'replacement'). 'K' is set to be 'L+10' by default. A too small choice of 'K' produces a bias of the data.
#' @param numCores specify the cores to be used for simulation. Default (NA) detects the available cores of the PC and uses all cores minus one.
#' @keywords Approximatex variance for replacement method
#' @references Bossert, P. (2023), 'Parameter estimation for second-order SPDEs in multiple space dimensions'
#' @export
#' @seealso [SecondOrderSPDEMulti::simulateSPDEmodelMulti]
#' @return A vector containing the approximated variance
#' @examples
#' d <- 2
#' M <- 10
#' theta0 <- 0
#' eta <- 1
#' nu <- c(2,1)
#' sigma <- 1
#' alphaDash <- 0.5
#' L <- 20
#' K <- 50
#'
#' va <- variance_approx(d,theta0,nu,eta,sigma,alphaDash,M,L,K)


variance_approx <- function(d,theta0,nu,eta,sigma,alphaDash,numberOfSpatialPoints,L,K,numCores=NA){

  if (!require("pacman")) {
    install.packages("pacman")
    require(pacman)
  } else {
    require(pacman)
  }

  is.naturalnumber <-
    function(x, tol = .Machine$double.eps^0.5)  x > tol & abs(x - round(x)) < tol

  if(!is.naturalnumber(d) || d <= 1){
    stop("Invalid input for 'd'. 'd' must be greter 1 and a natural number.")  }

  if(length(nu) != d){
    stop("'nu' must be a numeric vector of length 'd'.")
  }
  if(length(L) != 1 || !is.naturalnumber(L)){
    stop("Invalid input for 'L' or 'K'. Both must be a natural number of length 1.")
  }
  if(!is.na(K)){
    if(length(K) != 1 || !is.naturalnumber(K)){
      stop("Invalid input for 'L' or 'K'. Both must be a natural number of length 1.")
    }
    if(K<=L){
      stop("'K' must be greater than 'L'")
    }
  }

  alpha <- d/2-1+alphaDash
  if(length(alpha) != 1 || alpha <= d/2-1 || alpha >= d/2){
    stop("alphaDash must be a numeric of length 1 with: 0 < 'alphaDash' < 1.")
  }
  if(length(eta) != 1 || eta <= 0){
    stop("'eta' must be a positive numeric of length 1.")
  }
  if(length(theta0) != 1 ){
    stop("'theta0' must be a numeric of length 1.")
  }
  if(length(sigma)!= 1 || sigma <= 0){
    stop("'sigma' must be a positive numeric of length 1.")
  }

  pacman::p_load(dplyr, pbmcapply,data.table)


  M <- numberOfSpatialPoints
  M_df <- do.call(CJ,replicate(d,seq(1,(M-1)),F))
  if(is.na(numCores)){
    numCores <- detectCores() - 1
  }

  lambdall <- function(mm,theta0,nu,eta){
    return(-theta0+sum(nu^2)/(4*eta)+pi^2*eta*sum(mm^2))
  }

  ImPlus <- function(m,K){
    l <- lapply(1:d, function(i){
      l2 <- lapply(0:K, function(l){
        m[i]+2*l*M
      })
      unlist(l2)
    })
    return(l)
  }

  ImMinus <- function(m,K){
    l <- lapply(1:d, function(i){
      l2 <- lapply(0:K, function(l){
        2*M-m[i]+2*l*M
      })
      unlist(l2)
    })
    return(l)
  }


  Immm <- function(m,L,K){
    plus <- ImPlus(m,K)
    minus <- ImMinus(m,K)
    Im_list <- list(plus,minus)
    MM <- rep(M,d)


    combinations <- do.call(CJ,replicate(d,1:2,F))
    l <- lapply(1:dim(combinations)[1],function(i){
      l2 <- lapply(1:d,function(l){
        comb <- as.matrix(combinations[i,])[1,]
        Im_list[[comb[l]]][[l]]
      })
      do.call(expand.grid,l2)
    })
    erg <- do.call(rbind,l)


    Im_L <- filter_all(erg,all_vars(.<L*M))
    Im_K <- dplyr::setdiff(erg, Im_L)
    return(Im_K)
  }

  varianceApprox <- function(Im_K){

    l <- lapply(1:dim(Im_K)[1], function(i){
      ll <- lambdall(mm = Im_K[i,],theta0 = theta0, nu = nu, eta = eta)
      sigma^2/(2*ll^(1+alpha))
    })
    return(sum(unlist(l)))
  }

    l2 <- pbmclapply(1:dim(M_df)[1],function(i){
    Im_K <- Immm(as.vector(as.matrix(M_df[i,])[1,]),L,K)
    varianceApprox(Im_K)
  },mc.cores = numCores)
  variance_approx_data <- unlist(l2)

  return(variance_approx_data)

}
