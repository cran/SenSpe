\name{snsp2mup}
\alias{snsp2mup}
\title{Two-biomarker unpaired comparison in specificity (or sensitivity) at a
  controlled sensitivity (or specificity) level}
\description{Point estimation and exact bootstrap-based inference}
\usage{snsp2mup(mkx, n1x, mky, n1y, s0, covp=0.95, fixsens=TRUE, lbmdisx=TRUE, lbmdisy=TRUE)}
\arguments{
  \item{mkx}{values of biomarker X, cases followed by controls.}
  \item{n1x}{case size of biomarker X.}
  \item{mky}{values of biomarker Y, cases followed by controls.}
  \item{n1y}{case size of biomarker Y.}
  \item{s0}{controlled level of sensitivity or specificity.}
  \item{covp}{norminal level of confidence intervals.}
  \item{fixsens}{fixing sensitivity if True, and specificity otherwise.}
  \item{lbmdisx}{larger value of biomarker X is more associated with cases if
  True, and controls otherwise.}
  \item{lbmdisy}{larger value of biomarker Y is more associated with cases if
  True, and controls otherwise.}
}
\value{
  \item{diff}{diff[1]: difference of empirical point estimates;
  diff[2]: difference of oscillating bias-corrected estimates.}
  \item{hvar}{exact bootstrap variance estimate for diff[1].}
  \item{wald_ci}{wald_ci[1,]: Wald confidence interval using diff[1];
  wald_ci[2,]: Wald confidence interval using diff[2].}
  \item{pct_ci}{percentile confidence interval.}
  \item{scr_ci}{scr_ci[1,]: score confidence interval using diff[1];
  scr_ci[2,]: score confidence interval using diff[2].}
  \item{zq_ci}{extension of the BTII in Zhou and Qin (2005,
  Statistics in Medicine 24, pp 465--477).}
}
\examples{
## simulate biomarker X with 100 cases and 100 controls
mkx <- c(rnorm(100,2,1),rnorm(100,0,1))
## simulate biomarker Y with 100 cases and 100 controls
mky <- c(rnorm(100,1,1),rnorm(100,0,1))

## compare specificity at controlled 0.95 sensitivity
est <- snsp2mup(mkx, 100, mky, 100, 0.95)
}
\references{Huang, Y., Parakati, I., Patil, D. H.,and Sanda, M. G. (2023).
Interval estimation for operating characteristic of continuous biomarkers with
controlled sensitivity or specificity, \emph{Statistica Sinica} 33, 193--214.}
\author{Yijian Huang}
