#' Plot realized volatility against the theoretical result
#'
#' Provides a plot for the realized volatility statistics against the theoretical results. In addition it provides a plot for the deviation of theoretical and empirical results.
#' @param data_list return of the function 'simulateSPDEmodelMulti' or a list containing the returns of the function 'simulateSPDEmodelMulti'.
#' @param spatialDelta  a real number greater than zero and less than 1/2 for selecting only the data points which are delta away from the Dirichlet boundary condition. The default is 0.05.
#' @param newColors if True, it changes the colors of the returned plots. Default is F.
#' @keywords Plotting realized volatilities
#' @references PhD thesis Bossert, P.
#' @export
#' @seealso [SecondOrderSPDEMulti::simulateSPDEmodelMulti],[SecondOrderSPDEMulti::SecondOrderSPDEMulti]
#' @return a list with 2 plots, the main plot displays the theoretical and empirical realized volatilities. The other plot contains its deviation. Also both data sets for creating these plots are provided.

#' @examples
#' d <- 2
#' N <- 1000
#' M <- 10
#' theta0 <- 0
#' eta <- 1
#' nu <- c(2,1)
#' sigma <- 1
#' alphaDash <- 0.5
#' L <- 20
#' K <- 50
#' va <- variance_approx(d,theta0,nu,eta,sigma,alphaDash,M,L,K)
#'
#' res <- simulateSPDEmodelMulti(d=d,theta0=theta0,nu=nu,eta=eta,sigma=sigma,
#' alphaDash=alphaDash,numberOfSpatialPoints=M,numberOfTemporalPoints=N,
#' L=L,approx_var=va)
#'
#' p <- RV_plot(res,0.05,nu,eta,sigma,alphaDash)
#' p$plot






RV_plot <- function(data_list,spatialDelta=0.05,newColors = F){

  erg <- data_list
  d <- dim(erg$SG)[2]
  delta <- spatialDelta
  param <- data_list$param
  nu <- param[paste("nu",1:d,sep="")]
  eta <- param["eta"]
  sigma <- param["sigma"]
  alphaDash <- param["alphaDash"]
  kappa <- nu/eta

  rgb2 <- function(r,g,b){
    return(rgb(r/255,g/255,b/255))
  }
  darkmint <- c( rgb2(123, 188, 176), rgb2(85, 156, 158), rgb2(58, 124, 137), rgb2(35, 93, 114), rgb2(18, 63, 90))
  adobeColorsDiscrete <- c("#323E40","#F2AB27","#BF6B04","#732002","#D95323","#254021","#2E5902","#024873")
  
  is.naturalnumber <-
    function(x, tol = .Machine$double.eps^0.5)  x > tol & abs(x - round(x)) < tol
  
  if(length(nu) != d){
    stop("'nu' must be a numeric vector of length 'd'.")
  }
  alpha <- d/2-1+alphaDash
  if(length(alpha) != 1 || alpha <= d/2-1 || alpha >= d/2){
    stop("alphaDash must be a numeric of length 1 with: 0 < 'alphaDash' < 1.")
  }
  if(length(eta) != 1 || eta <= 0){
    stop("'eta' must be a positive numeric of length 1.")
  }
  if(length(sigma)!= 1 || sigma <= 0){
    stop("'sigma' must be a positive numeric of length 1.")
  }

  y_WB_index <- function(SG,delta){
    d <- dim(SG)[2]
    l <- lapply(1:dim(SG)[1],function(i){
      if(sum(SG[i,]>= delta & SG[i,]<= 1-delta)==d){i}
    })
    return(unlist(l))
  }

  if (!require("pacman")) {
    install.packages("pacman")
    require(pacman)
  } else {
    require(pacman)
  }
  pacman::p_load(pbmcapply,ggplot2)
  numCores <- detectCores() - 1


  SGWB <- erg$SG
  SGWB2 <- y_WB_index(erg$SG,delta)

  RVVV <- function(data,SG){
    N <- dim(data)[2]-1
    l <- pbmclapply(SG, function(j){
      erg <- c()
      for (i in 2:dim(data)[2]) {
        erg <- c(erg,(data[j,i]-data[j,i-1])^2)
      }
      sum(erg)/(N*(1/(N+1))^(alphaDash))
    },mc.cores = numCores)
    return(unlist(l))

  }
  x1 <- RVVV(erg$data,SGWB2)



  theoretical <- function(SG,SGIndex){
    l <- lapply(SGIndex, function(j){
      sigma^2*exp(-sum(kappa*SG[j,]))*gamma(1-alphaDash)/(2^d*(pi*eta)^(d/2)*alphaDash*gamma(d/2))
    })
    return(unlist(l))
  }
  x2 <- theoretical(erg$SG,SGWB2)

  y <- 1:length(SGWB2)
  dat <- data.frame(x=y,y=c(x1,x2),group = as.factor(rep(c("Empirical","Theoretical"),each=length(y))))
  if(newColors){
    p <- ggplot(dat,aes(x=x,y=y,group=group,color = group))+geom_line()+
      theme_minimal()+labs(x="Index of filtered spatial grid",y="Rescaled realized volatility",color=NULL)+
      theme(legend.position = "bottom")+
      scale_fill_manual(values=adobeColorsDiscrete)
  } else {
    p <- ggplot(dat,aes(x=x,y=y,group=group,color = group))+geom_line()+
      theme_minimal()+labs(x="Index of filtered spatial grid",y="Rescaled realized volatility",color=NULL)+
      theme(legend.position = "bottom")
  }


  dat$index <- SGWB2
  dat2 <- data.frame(deviation = abs(x1-x2),index = SGWB2)

  
  if(newColors){
    p2 <- ggplot(dat2,aes(x=index,y=deviation,color=deviation))+geom_line()+
      theme_minimal()+labs(x="Index of filtered spatial grid",y="Deviation of theoretical and empirical results")+
      theme(legend.position = "none")+
      scale_colour_gradientn(colours = darkmint)
  } else {
    p2 <- ggplot(dat2,aes(x=index,y=deviation,color=deviation))+geom_line()+
      theme_minimal()+labs(x="Index of filtered spatial grid",y="Deviation of theoretical and empirical results")+
      theme(legend.position = "none")
  }
  return(list(SG=SGWB , data = dat,deviation = dat2,plot = p,deviation_plot = p2))
}
