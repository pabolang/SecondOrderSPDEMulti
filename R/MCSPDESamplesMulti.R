#' Simulates multiple second-order stochastic partial differential equations in multiple space dimension using Monte-Carlo repetitions
#'
#' Simulate Monte-Carlo samples of a SPDE model in multiple space dimension on a discrete grid with \code{N} temporal and \code{M} spatial points, where the grid points are equidistant within the unit square. Each repetition uses the "replacement"-method by default, see [SecondOrderSPDEMulti::simulateSPDEmodelMulti].
#' @param repetitions A natural number giving the number of the Monte-Carlo repetitions.
#' @param ... for further arguments see [SecondOrderSPDEMulti::simulateSPDEmodelMulti].
#' @keywords Monte-Calro samples SPDE.
#' @references PhD thesis Bossert, P.
#' @export
#' @seealso [SecondOrderSPDEMulti::simulateSPDEmodelMulti],[SecondOrderSPDEMulti::SecondOrderSPDEMulti].
#' @return A list of the returns of [SecondOrderSPDEMulti::simulateSPDEmodelMulti].
#' @examples
#' reps <- 10
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
#' MCSPDESamplesMulti(repetitions=reps,d=d,theta0=theta0,nu=nu,eta=eta,sigma=sigma,alphaDash=alphaDash,numberOfSpatialPoints=M,numberOfTemporalPoints=N,L=L,approx_var=va)

MCSPDESamplesMulti <- function(repetitions,saveData=F,start=1,save_path="",d,theta0,nu,eta,sigma,alphaDash,numberOfSpatialPoints,numberOfTemporalPoints,L,K=NA,em_data=NA,
                          approx_var=NA,
                          save_approx_var=F,path_approx_var="",
                          save_em_data=F,path_em_data="",
                          mehtod="replacement",
                          numCores = NA){

  if(saveData){
    end <- repetitions + start -1
    l <- lapply(start:end, function(i){
      res <- simulateSPDEmodelMulti( d,theta0,nu,eta,sigma,alphaDash,numberOfSpatialPoints,numberOfTemporalPoints,L,K,em_data,
                              approx_var,
                              save_approx_var,path_approx_var,
                              save_em_data,path_em_data,
                              mehtod,numCores)
      saveRDS(res,paste(save_path,"/dat",i,".RDS",sep=""))
    })
  } else {
    MC_list <- lapply(1:repetitions, function(i){
      simulateSPDEmodelMulti( d,theta0,nu,eta,sigma,alphaDash,numberOfSpatialPoints,numberOfTemporalPoints,L,K,em_data,
                              approx_var,
                              save_approx_var,path_approx_var,
                              save_em_data,path_em_data,
                              mehtod,numCores)
    })
    return(MC_list)
  }


}


