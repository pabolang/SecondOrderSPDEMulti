# SecondOrderSPDEMulti

Simulate and plot multi-dimensional linear second-order SPDE models. This R-package also provides estimation methods for estimating and the natural parameters in this model.

## Installation
```r
require("devtools")
devtools::install_github('pabolang/ParabolicSPDEsMulti')
```

## Usage
We consider the following second-order stochastic partial differential equation 
$$\text{d} X_t(\textbf{y}) = A_\vartheta X_t(\textbf{y})\text{d} t+\sigma\text{d} B_t(\textbf{y}),$$

where $\vartheta_0,\nu_1,\ldots,\nu_d\in\mathbb{R}$, $\eta,\sigma>0$ 
and a cylindrical Brownian motion $B_t(y)$ within a Sobolev space 
on the spatial domain $\[0,1\]^d$. The corresponding differential operator is given by
$$A_\vartheta = \eta \sum_{l=1}^d \frac{\partial}{\partial y_l^2}+\sum_{l=1}^d \nu_l\frac{\partial}{\partial y_l}+\vartheta_0. $$
Further, we consider a Dirichlet boundary condition and an initial condition $\xi\equiv 0$.


By using this package, we can simply simulate a SPDE model on an equidistant discrete $N\times M^d$ grid, 
where $N$ denotes the temporal and $M$ the spatial resolution on each of the $d$ spatial axes, using the function `simulateSPDEmodel`:
```r
library(ParabolicSPDEs;ulti)
d <- 2
N <- 1000
M <- 10
theta0 <- 0
eta <- 1
nu <- c(2,1)
sigma <- 1
alphaDash <- 0.5
L <- 20
K <- 50
va <- variance_approx(d,theta0,nu,eta,sigma,alphaDash,M,L,K)

res <- simulateSPDEmodelMulti(d=d,theta0=theta0,nu=nu,eta=eta,sigma=sigma,
alphaDash=alphaDash,numberOfSpatialPoints=M,numberOfTemporalPoints=N,
L=L,approx_var=va)
res
```


The function `simulateSPDEmodel` returns a list, containg among others a matrix of the simultated model, which we can be plotted using the following function:
```r
SPDE_plot(data_list = res, coord_plot = 2, spatialCoordsRemainingAxes = c(0.3,0.7))
```

<img width="983" alt="Bildschirmfoto 2023-07-22 um 17 21 01" src="https://github.com/pabolang/ParabolicSPDEsMulti/assets/78961989/1ad819a1-7a67-46ca-be5a-869b0cc82a52">



For creating multiple SPDE samples, use the function `MCSPDESamplesMulti`. 
This package also includes the function `estimateParametersSPDEMulti` for estimating the parameters of a SPDE model. For the respective statistical assumptions, see the documentation.
