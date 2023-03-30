\name{snsp2mp}
\alias{snsp2mp}
\title{Two-biomarker paired comparison in specificity (or sensitivity) at a
  controlled sensitivity (or specificity) level}
\description{Point estimation and exact bootstrap-based inference}
\usage{snsp2mp(mk, n1, s0, covp=0.95, fixsens=TRUE, lbmdis=TRUE)}
\arguments{
  \item{mk}{Each of two rows corresponds to a biomarker, cases followed by
  controls.}
  \item{n1}{case size.}
  \item{s0}{controlled level of sensitivity or specificity.}
  \item{covp}{norminal level of confidence intervals.}
  \item{fixsens}{fixing sensitivity if True, and specificity otherwise.}
  \item{lbmdis}{larger value of a biomarker is more associated with cases if
  True, and controls otherwise.}
}
\value{
  \item{diff}{diff[1]: difference of empirical point estimates;
  hss[2]: difference of oscillating bias-corrected estimates.}
  \item{btmn}{bootstrap mean of the empirical difference.}
  \item{btva}{exact bootstrap variance estimate for diff[1].}
  \item{btdist}{exact bootstrap probability mass function at (-n0:n0)/n0 with n0
  being the size of controls if sensitivity is controlled, or at (-n1:n1)/n1
  otherwise.}
  \item{wald_ci}{wald_ci[1,]: Wald confidence interval using diff[1];
  wald_ci[2,]: Wald confidence interval using diff[2].}
  \item{pct_ci}{percentile confidence interval.}
  \item{scr_ci}{scr_ci[1,]: score confidence interval using diff[1];
  scr_ci[2,]: score confidence interval using diff[2].}
  \item{zq_ci}{extension of the BTII in Zhou and Qin (2005,
  Statistics in Medicine 24, pp 465--477).}
}
\examples{
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
est <- snsp2mp(mk, 100, 0.95)
}
\references{Huang, Y., Parakati, I., Patil, D. H.,and Sanda, M. G. (2023).
Interval estimation for operating characteristic of continuous biomarkers with
controlled sensitivity or specificity, \emph{Statistica Sinica} 33, 193--214.}
\author{Yijian Huang}
