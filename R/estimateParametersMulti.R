#' Estimate the parameters of a secon-order SPDE model in multiple space dimensions
#'
#' Oracle and non-oracle estimations for SPDE models. All estimators are consistent and estimate the natural parameters of the d-dimensional SPDE model. Non-oracle estimations for second-order SPDEs are based on a log-linear model.
#' @param data_list return of the function 'simulateSPDEmodelMulti' or a list containing multiple returns of the function 'simulateSPDEmodelMulti'.
#' @param estimationMethod a string indicating the parameter/parameters to be estimated. If only sigma is unknown choose \code{"OracleSigma"} and provide \code{delta,alphaDash,kappa,eta} respectively.
#' If alphaDash is unknown, choose 'alphaDash'. No additional information is needed.
#' If \code{sigma and kappa} are known choose \code{"SigmaAndKappa"} and provide the known parameter \code{alphaDash} as well as a \code{indexset} satisfying the respective Assumptions.
#' If none of the parameters are known, choose \code{"alphaDash"} or \code{"all"}. For \code{"alphaDash"} provide 'spatialDelta' and for \code{"all"} provide 'spatialDelta' and 'indexset'.
#' @param spatialDelta a real number greater than zero and less than 1/2 for selecting only the data points which are delta away from the Dirichlet boundary condition. The default is 0.05.
#' @param ... further arguments depending on the chosen estimation method. See \code{estiamtionmethod} or the examples below.
#' @keywords Parameter Estimation for SPDEs.
#' @references PhD thesis  Bossert, P.
#' @export
#' @return a named numeric vector denoting the respective estimate. For method \code{"OracleSigma"} the returned value denotes the estimation of 'sigma^2'. For method \code{"SigmaAndKappa"} the returned value denotes the estimated quotient theta1/theta2 of the respective axis.
#' For the method \code{"alphaDash"} returns the estimate for 'alphaDash'. \code{"all"} first estimates 'alphaDash' and based on this estimate 'sigma0^2' and each 'kappa'.
#' See references for details on estimation methods.
#' @seealso [SecondOrderSPDEMulti::simulateSPDEmodelMulti], [SecondOrderSPDEMulti::MCSPDESamplesMulti],[SecondOrderSPDEMulti::SecondOrderSPDEMulti], 
#' @examples
#' d <- 2
#' N <- 10000
#' M <- 10
#' theta0 <- 0
#' eta <- 1
#' nu <- c(2,1)
#' sigma <- 1
#' alphaDash <- 0.5
#' L <- 20
#' K <- 50
#'
#' res <- simulateSPDEmodelMulti(d=d,theta0=theta0,nu=nu,eta=eta,sigma=sigma,
#' alphaDash=alphaDash,numberOfSpatialPoints=M,
#' numberOfTemporalPoints=N,L=L,K=K)
#'
#' estimateParametersSPDEMulti(res,estimationMethod="OracleSigma",spatialDelta=0.05,alphaDash=alphaDash,kappa=nu/eta,eta=eta)
#' estimateParametersSPDEMulti(res,estimationMethod="SigmaAndKappa",alphaDash=alphaDash, indexset = c(16,53,57))
#' estimateParametersSPDEMulti(res,estimationMethod="alphaDash",spatialDelta=0.05)
#' estimateParametersSPDEMulti(res,estimationMethod="all",spatialDelta=0.05,indexset = c(16,53,57))





estimateParametersSPDEMulti <- function(data_list,estimationMethod,kappa=NA,eta=NA,alphaDash=NA,spatialDelta=0.05,indexset=NA,ignoreWarnings=F){

  if(!(estimationMethod %in% c("OracleSigma","SigmaAndKappa","alphaDash","all")) ){
    stop("'estimationMethod' must be either 'OracleSigma', 'SigmaAndKappa', 'alphaDash', or 'all'.")
  }
  if (!require("pacman")) {
    install.packages("pacman")
    require(pacman)
  } else {
    require(pacman)
  }
  pacman::p_load(dplyr, pbmcapply,data.table,MASS)
  numCores <- detectCores() - 1




  dat <- data_list$data
  SG <- data_list$SG
  d <- dim(SG)[2]
  delta <- spatialDelta


  y_WB_df <- function(SG,delta){
    y <- filter_all(SG,all_vars(.>=delta & .<=1-delta))
    return(y)
  }

  y_WB_index <- function(SG,delta){
    d <- dim(SG)[2]
    l <- lapply(1:dim(SG)[1],function(i){
      if(sum(SG[i,]>= delta & SG[i,]<= 1-delta)==d){i}
    })
    return(unlist(l))
  }


  RV <- function(yIndex,dat){
    datPath <- dat[yIndex,]
    n <- length(datPath)
    datPath1 <- datPath[-1]
    datPath2 <- datPath[-n]

    erg <- sum((datPath1-datPath2)^2)
    return(erg)
  }

  RV_thinGrid <- function(yIndex,dat){
    if((dim(dat)[2]-1)%%2 != 0 ){
      dat <- dat[,1:dim(dat)[2]-1]
    }
    datPath <- dat[yIndex,]
    n <- length(datPath)
    selection <- seq(1,n,2)
    datPath2 <- datPath[selection]
    n2 <- length(datPath2)
    datPath_1 <- datPath2[-1]
    datPath_2 <- datPath2[-n2]

    erg <- sum((datPath_1-datPath_2)^2)
    return(erg)
  }

  RV_rescaled <- function(yIndex,dat,alphaDash){
    rv <- RV(yIndex,dat)
    N <- dim(dat)[2]-1
    return(rv/(N*(1/N)^(alphaDash)))
  }
  # RV_rescaled(30,erg$data,0.5)


  V <- function(dat,alphaDash,SG,kappa,delta){
    N <- dim(dat)[2]-1
    index <- y_WB_index(SG,delta)
    l <- lapply(index,function(i){
      RV(i,dat)*exp(sum(kappa*SG[i,]))
    })
    erg <- sum(unlist(l))
    rescale <- N*(1/N)^alphaDash
    return(erg/rescale)
  }

  hat_sigma_squared <- function(data,SG,delta,alphaDash,kappa,eta){
    N <- dim(data)[2]-1
    d <- dim(SG)[2]
    M <- length(y_WB_index(SG,delta))
    sp <- V(data,alphaDash,SG,kappa,delta)
    fac <- gamma(d/2)*alphaDash*2^d*(pi*eta)^(d/2)/(M*gamma(1-alphaDash))
    return(sp*fac)
  }

  hat_sigma_squared_y <- function(data,yIndex,SG,alphaDash,kappa,eta){
    N <- dim(data)[2]-1
    d <- dim(SG)[2]
    sp <- RV(yIndex,dat)*exp(sum(kappa*SG[yIndex,]))/(N*(1/N)^alphaDash)
    fac <- gamma(d/2)*alphaDash*2^d*(pi*eta)^(d/2)/(gamma(1-alphaDash))
    return(sp*fac)
  }

  hat_psi <- function(data,SG,indexset,alphaDash){
    X <- as.matrix(cbind(1,SG[indexset,]))
    A <- t(X) %*% X
    if(round(det(A),12) == 0) {
      stop("Input for 'indexMatrix' is not a spanning set of (d+1)-dimensional real numbers.")
    }
    n <- dim(data)[2]-1
    m <- length(indexset)
    d <- dim(SG)[2]

    if(!ignoreWarnings){
      if(n<= m^((d+2)/alphaDash)){
        stop("Statistical Assumptions are not satisfied! Length of 'indexset' with a  power of ((d/2)/'alphaDash') must  less than N")
      }
    }

    l <- lapply(indexset, function(i){
      log(RV_rescaled(i,data,alphaDash))
    })
    Y <- as.matrix(unlist(l),ncol = 1)
    return(MASS::ginv(A)%*% t(X) %*% Y)
  }

  hat_beta <- function(data,SG,indexset,alphaDash){
    d <- dim(SG)[2]
    K <- gamma(1-alphaDash)/(alphaDash*2^d*pi^(d/2)*gamma(d/2))
    erg <- hat_psi(data,SG,indexset,alphaDash)
    erg1 <- exp(erg[1])/K
    erg2 <- -erg[2:(d+1)]
    return(c(erg1,erg2))
  }

  hat_alpha <- function(data,SG,delta){
    index <- y_WB_index(SG,delta)
    l <- lapply(index, function(i){
      1+log(RV_thinGrid(i,data)/RV(i,data))/(log(2))
    })
    return(sum(unlist(l))/length(index))

  }



  if(estimationMethod == 'OracleSigma'){
    se <- hat_sigma_squared(dat,SG,delta,alphaDash,kappa,eta)
    names(se) <- "sigma_0Squared Est"
    return(list(estimate = se, filteredSGIndices = y_WB_index(SG,delta)))
  }
  if(estimationMethod == 'SigmaAndKappa'){
    be <- hat_beta(dat,SG,indexset,alphaDash)
    names(be) <- c("sigma_0Squared Est",paste("kappa",1:d, "Est",sep=" "))
    return(be)
  }
  if(estimationMethod == 'alphaDash'){
    ae <- hat_alpha(dat,SG,delta)
    names(ae) <- "alphaDash Est"
    return(list(estimate = ae, filteredSGIndices = y_WB_index(SG,delta)))
  }
  if(estimationMethod == 'all'){
    ad <- hat_alpha(dat,SG,delta)
    beta <- hat_beta(dat,SG,indexset,ad)
    all_est <- c(ad,beta)
    names(all_est) <- c("alphaDash Est","sigma_0Squared Est",paste("kappa",1:d, "Est",sep=" "))
    return(list(estimate = all_est, filteredSGIndicesForAlpha = y_WB_index(SG,delta)))
  }
}
