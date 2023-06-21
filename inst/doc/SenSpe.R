## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE, message=FALSE, warning=FALSE------------------------
#  install.packages("SenSpe")

## ----snsp1m, eval=TRUE, message=FALSE, warning=FALSE--------------------------
library("SenSpe")
## simulate biomarkers of 100 cases and 100 controls
n1 <- 100
n0 <- 100
mk <- c(rnorm(n1,1,1),rnorm(n0,0,1))
## estimate specificity at controlled 0.95 sensitivity
snsp1m(mk, n1=n1, s0=0.95)

## ----snsp2mp, eval=TRUE, message=FALSE, warning=FALSE-------------------------
## simulate paired biomarkers X and Y, with correlation 0.5, 100 cases and 100 controls
n1 <- 100
n0 <- 100
rho <- 0.5
mkx <- rnorm(n1+n0,0,1)
mky <- rho*mkx + sqrt(1-rho^2)*rnorm(n1+n0,0,1)
mkx <- mkx + c(rep(2,n1),rep(0,n0))
mky <- mky + c(rep(1,n1),rep(0,n0))
mk <- rbind(mkx,mky)
## compare specificity at controlled 0.95 sensitivity
snsp2mp(mk, 100, 0.95)

## ----snsp2mup, eval=TRUE, message=FALSE, warning=FALSE------------------------
## simulate biomarker X with 100 cases and 100 controls
mkx <- c(rnorm(100,2,1),rnorm(100,0,1))
## simulate biomarker Y with 100 cases and 100 controls
mky <- c(rnorm(100,1,1),rnorm(100,0,1))
## compare specificity at controlled 0.95 sensitivity
snsp2mup(mkx, 100, mky, 100, 0.95)

