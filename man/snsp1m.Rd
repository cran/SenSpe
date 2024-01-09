\name{snsp1m}
\alias{snsp1m}
\title{Estimating specificity (or sensitivity) at a controlled sensitivity
  (or specificity) level}
\description{Point estimation and exact bootstrap-based inference}
\usage{snsp1m(mk, n1, s0, covp=0.95, fixsens=TRUE, lbmdis=TRUE)}
\arguments{
  \item{mk}{biomarker values of cases followed by controls.}
  \item{n1}{size of cases.}
  \item{s0}{controlled level of sensitivity or specificity.}
  \item{covp}{norminal level of confidence intervals.}
  \item{fixsens}{fixing sensitivity if True, and specificity otherwise.}
  \item{lbmdis}{larger biomarker value is more associated with cases if True,
  and controls otherwise.}
}
\value{
  \item{threshold}{estimated threshold, at and beyond which the empirical
  sensitivity or specificity is the smallest no less than the controlled level
  s0.}
  \item{hss}{hss[1]: empirical point estimate of specificity at controlled
  sensitivity, or vice versa; hss[2]: oscillating bias-corrected estimate.}
  \item{hvar1}{estimated variance component from cases if specificity at
  controlled sensitivity is estimated, or from controls otherwise.}
  \item{hvar2}{estimated variance component from controls if specificity at
  controlled sensitivity is estimated, or from cases otherwise.}
  \item{hvar}{exact bootstrap variance estimate, =hvar1+hvar2.}
  \item{btpdf}{exact bootstrap probability mass function at (0:n0)/n0 with n0
  being the size of controls if sensitivity is controlled, or at (0:n1)/n1
  otherwise.}
  \item{wald_ci}{wald_ci[1,]: Wald confidence interval using hss[1];
  wald_ci[2,]: Wald confidence interval using hss[2].}
  \item{pct_ci}{percentile confidence interval.}
  \item{scr_ci}{scr_ci[1,]: score confidence interval using hss[1];
  scr_ci[2,]: score confidence interval using hss[2].}
  \item{zq_ci}{exact bootstrap version of the BTII in Zhou and Qin (2005,
  Statistics in Medicine 24, pp 465--477).}
}
\examples{
## simulate biomarkers of 100 cases and 100 controls
set.seed(1234)
mk <- c(rnorm(100,1,1),rnorm(100,0,1))

## estimate specificity at controlled 0.95 sensitivity
est <- snsp1m(mk, 100, 0.95)
}
\references{Huang, Y., Parakati, I., Patil, D. H.,and Sanda, M. G. (2023).
Interval estimation for operating characteristic of continuous biomarkers with
controlled sensitivity or specificity, \emph{Statistica Sinica} 33, 193--214.}
\author{Yijian Huang}
