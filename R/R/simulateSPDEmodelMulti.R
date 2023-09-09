#' Simulate a second-order stochastic partial differential equations in multiple space dimension
#'
#' Simulate a sample of a SPDE model in d space dimension on a discrete grid with \code{N}-temporal and \code{M}-spatial (each axis) points, where the grid points are equidistant within the unit square. The initial condition is set to be zero. For simulating SPDE samples with a general initial condition use the build-in cut-off method. Note that this method results in a systematic bias and dramatically increases computational costs. Furthermore, the SPDE model is using the Dirichlet boundary condition. The multi-dimensional replacement method only approximates the replacement variance. Therefore the optional parameter 'K' appears. A too small choice of 'K' results in an bias of the simulated data.
#' @param d a natural number which denotes the spacial dimension. 'd' must be greater 1.
#' @param theta0 a real number which controls the drift of the solution field.
#' @param nu a real vector of dimension 'd' which controls the curvature of the solution field on each axis respectively.
#' @param eta a real number greater than zero which reduces the noise level of parameter sigma and the curvature effect of parameter \code{nu}.
#' @param sigma a real number greater than zero which controls the overall noise level of the solution field.
#' @param alphaDash a real number in \code{(0,1)} controlling the roughness of the sample paths.
#' @param numberOfSpatialPoints number of equidistant spatial points M on each axis, where we consider the domain \code{[0,1]}.
#' @param numberOfTemporalPoints number of equidistant temporal points N in \code{[0,1]}.
#' @param L a natural number indicating the replacement bound LM dependent on multiples of the d-dimensional vector with entries \code{M}. The default is \code{L=10}.
#' @param xi initial condition. The default is \code{xi = 0}. For general initial condition choose \code{method = "cutoff"}.
#' @param cutoff a natural number for the cut-off frequency of the Fourier series. Only used when \code{method = "cutoff"}. The default is \code{cutoff=10000}.
#' @param method either \code{"replacement"} (Default) or \code{"cutoff"}. Note, that the multi-dimensional replacement method only allows for the initial condition to be zero. For general initial conditions choose \code{method="cutoff"}.
#' @param K a natural number greater 'L' denoting the cut-off if the variance calculation for the replacement random variables (only of use when method is 'replacement'). 'K' is set to be 'L+10' by default. A too small choice of 'K' produces a bias of the data.
#' @param em_data a data frame containing every combination of the orthonormal system which is used for the Fourier transform. Alternatively a path of type String can be used.
#' @param approx_var a vector especially useful when simulating more than one data set. The approximated variance for the replacement random variables can be computed via the function [SecondOrderSPDEMulti::variance_approx]. Alternatively a path of type String to the approx_var type vector can be used. Only of use when method is 'replacement'.
#' @param save_approx_var a Boolean variable being 'F' by default. If being 'T', the generated approximated variance or tghe replacement random varaibles is going to be saved. Therefore the parameter 'name_approx_var' needs to be specified as a string. Only of use when method is 'replacement'.
#' @param path_approx_var a string for locating the folder of the approximated variance to be saved. Only being used, if 'save_approx_var' is true. Only of use when method is 'replacement'.
#' @param save_em_data a Boolean variable being 'F' by default. If being 'T', the generated orthonormal system is going to be saved. Therefore the parameter 'name_em_data' needs to be specified as a string. Only of use when method is 'replacement'.
#' @param path_em_data a string for locating the folder of the data to be saved. Only being used, if 'save_em_data' is true. Only of use when method is 'replacement'.
#' @keywords Sample of one SPDE in multiple space dimensions
#' @references PhD thesis Bossert, P.
#' @export
#' @seealso [SecondOrderSPDEMulti::MCSPDESamplesMulti],[SecondOrderSPDEMulti::plot_SPDEMulti],[SecondOrderSPDEMulti::RV_plot],[SecondOrderSPDEMulti::SecondOrderSPDEMulti].
#' @return A list containing 'data', 'SG', 'TG' and 'param'. 'data' is a (M+1)^d x (N+1) data-frame of the simulated SPDE model. 'SG' contains the generated spatial grid, where 'TG' the generated temporal grid. 'param' contains the parameters used for the simulation.
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


simulateSPDEmodelMulti <- function(d,theta0,nu,eta,sigma,alphaDash,numberOfSpatialPoints,numberOfTemporalPoints,L,
                                   K=NA,em_data=NA,
                                   approx_var=NA,
                                   save_approx_var=F,path_approx_var="",
                                   save_em_data=F,path_em_data="",
                                   method="replacement"){
  if (!require("pacman")) {
    install.packages("pacman")
    require(pacman)
  } else {
    require(pacman)
  }
  pacman::p_load(dplyr, pbmcapply,data.table)
  numCores <- detectCores() - 1

  is.naturalnumber <-
    function(x, tol = .Machine$double.eps^0.5)  x > tol & abs(x - round(x)) < tol

  read_em_data <- function(name){
    erg <- readRDS(paste(name,".RDS",sep=""))
    return(erg)
  }

  read_approx_var <- function(name){
    erg <- readRDS(paste(name,".RDS",sep=""))
    return(erg)
  }


  M <- numberOfSpatialPoints
  N <- numberOfTemporalPoints


  # check input
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




  if(method == "replacement"){
    print("Using method: 'replacement'")




    kappa <- nu/eta
    sigma0 <- sigma/eta
    Delta <- 1/N
    Mbound <- rep(M,d)

    y1Dim <- seq(0,1,1/M)
    SG <- do.call(CJ,replicate(d,y1Dim,F))
    M_df <- do.call(CJ,replicate(d,seq(1,(M-1)),F))
    TG <- seq(0,1,1/N)

    dimSG <- dim(SG)[1]
    dimM_df <- dim(M_df)[1]
    dimTG <- length(TG)

    em <- function(d,kappa,y,m){
      erg <- 2^(d/2)*prod(exp(-(kappa*y)/2))*prod(sin(pi*m*y))
      return(erg)
    }

    # create data
    if(sum(is.na(em_data))==1){
      print("Simulate Fourier Basis:")
      l1 <- pbmclapply(1:dimM_df,function(m){
        l2 <- lapply(1:dimSG, function(j){
          em(d,kappa,as.matrix(SG[j,])[1,],as.matrix(M_df[m,])[1,])
        })
        unlist(l2)
      },mc.cores = numCores)

      # 1 Zeile = M_df[1,] und alle y etc.
      em_data <- do.call(rbind,l1)
      if(save_em_data){
        saveRDS(em_data,paste(path_em_data,".RDS",sep=""))
      }
    }
    if(is.character(em_data)){
      em_data <- read_em_data(em_data)
    }

    # coordinate process

    lambdal <- function(l,theta0,nu,eta){
      return(-theta0+sum(nu^2)/(4*eta)+pi^2*eta*sum(l^2))
    }


    xl <- function(N,l,sigma,Delta,alpha,theta0,nu,eta){
      erg <- rep(NULL,N+1)
      erg[1] <- 0
      r <- rnorm(N)
      ll <- lambdal(l,theta0,nu,eta)
      b <- sigma*sqrt((1-exp(-2*ll*Delta))/(2*ll^(1+alpha)))*r

      for (i in 2:(N+1)) {
        a <- erg[i-1]*exp(-ll*Delta)+b[i-1]
        erg <- c(erg,a)
      }
      return(erg)
    }


    # Im
    # m <- as.matrix(M_df[123,])[1,]
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

    Imm <- function(m,L){
      plus <- ImPlus(m,L)
      minus <- ImMinus(m,L)
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

      Im_L <- dplyr::filter_all(erg,all_vars(.<L*M))
      return(Im_L)
    }







    # calc approx variance for all m
    if(sum(is.na(approx_var))==1){
      print("Calculate Approximated Variance:")
      if(is.na(K)){K=L+10}
      variance_approx_data <- variance_approx(d,theta0,nu,eta,sigma,alphaDash,numberOfSpatialPoints,L,K)

      if(save_approx_var){
        saveRDS(variance_approx_data,paste(path_approx_var,".RDS",sep=""))
      }
    } else {
      if(is.character(approx_var)){
        variance_approx_data <- read_approx_var(approx_var)
      } else {
        variance_approx_data <- approx_var
      }
    }




    xl_Im <- function(m,Im_L){
      l <- lapply(1:dim(Im_L)[1], function(i){
        xl(N,Im_L[i,],sigma,Delta,alpha,theta0,nu,eta)
      })
      erg <- do.call(rbind,l)
      return(colSums(erg))
    }



    print("Generate coordinate processes 1:")
    l <- pbmclapply(1:dimM_df,function(m){
      Im_L <- Imm(as.matrix(M_df[m,])[1,],L)
      xl_Im(M_df[m,],Im_L)
    },mc.cores = numCores)

    xm_data <- do.call(rbind,l)


    Um <- function(dimM_df){
      l <- pbmclapply(1:dimM_df, function(i){
        r <- c(0,rnorm(N,0,sqrt(variance_approx_data[i])))
        xm_data[i,]+r
      },mc.cores = numCores)
      return(do.call(rbind,l))
    }
    print("Generade coordinate process 2:")
    Um_data <- Um(dimM_df)
    #plot(xm_data[2,],type="l")
    #plot(Um_data[2,],type="l")


    print("Merge Components:")
    l <- pbmclapply(1:dimTG, function(i){
      l2 <- lapply(1:dimSG, function(j){
        l3 <- lapply(1:dimM_df, function(m){
          Um_data[m,i]*em_data[m,j]
        })
        sum(unlist(l3))
      })
      unlist(l2)
    },mc.cores = numCores)
    result <- do.call(cbind,l)

    param <- c(d,theta0,nu,eta,sigma,alphaDash,numberOfSpatialPoints,numberOfTemporalPoints,L,
               K)
    names(param) <- c("d","theta0",paste("nu",1:d,sep=""),"eta","sigma","alphaDash","numberSP","numberTP","L","K")

    return(list(data = result, SG = SG, TG = TG, param = param))

  } else {
    print("Using method: 'Cut-off'")
    print("Warning: implementation must be varified! Please use 'replacement' method!")
    pacman::p_load(parallel,pbapply,cubature, pbmcapply,data.table)
    numCores <- detectCores() - 1

    y <- seq(0,1,1/M)
    t <- seq(0,1,1/N)
    xiZero=T
    k_df <- do.call(CJ, replicate(d, 1:K, FALSE))
    names(k_df) <- paste("Dim",1:d)
    y_df <- do.call(CJ, replicate(d, y, FALSE)) # spatial coords not index
    names(y_df) <- paste("Dim",1:d)
    k_numberComb <- dim(k_df)[1]
    y_numberComb <- dim(y_df)[1]
    kappa <- nu/eta

    # x_kMat: [i,] gives all values for fixed time t_i and all k {1,...,K}^d.
    l_x <- pbmclapply(1:k_numberComb,function(kk){
      k1 <- as.vector(t(k_df[kk,]))
      x_k(k1,t,d,alpha,theta0,eta,nu,sigma,xiZero,xi)
    },mc.cores = numCores)
    x_kMat <- do.call(rbind,l_x)

    l_e <- pbmclapply(1:y_numberComb, function(yy){
      y1 <- as.vector(t(y_df[yy,]))
      l1 <- lapply(1:k_numberComb, function(kk){
        k1 <- as.vector(t(k_df[kk,]))
        e_k(y=y1,kappa=kappa,k=k1,d=d)
      })
      unlist(l1)
    },mc.cores = numCores)
    e_kMat <- do.call(rbind,l_e)


    l_y <- pbmclapply(1:y_numberComb,function(yy){
      l_t <- lapply(1:(N+1),function(tt){
        l_k <- lapply(1:k_numberComb, function(kk){
          x_kMat[kk,tt]*e_kMat[yy,kk]
        })
        sum(unlist(l_k))
      })
      unlist(l_t)
    },mc.cores = numCores)

    dat <- do.call(cbind,l_y)
    # return dat with [i,] fixed time t_i and all space coords
    #                 [,j] fixed spatial coord y_j and all time coord



    return_list <- list(data = dat, time_coords = t, space_coords = y_df, k_indices = k_df )
    return(return_list)
  }



}



