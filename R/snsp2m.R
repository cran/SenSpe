snsp2mup <- function(mkx, n1x, mky, n1y, s0, covp=0.95, fixsens=TRUE,
                     lbmdisx=TRUE, lbmdisy=TRUE)
{  estx <- snsp1m(mkx, n1x, s0, covp, fixsens, lbmdisx)
   esty <- snsp1m(mky, n1y, s0, covp, fixsens, lbmdisy)

   n0x <- length(mkx)-n1x
   n0y <- length(mky)-n1y
   zq <- qnorm(1-(1-covp)/2)

   diff <- estx$hss-esty$hss
   hvar <- estx$hvar+esty$hvar
   wald_ci <- NULL
   for(i in 1:2) {
     tmp <- diff[i]+c(-1,1)*zq*sqrt(hvar)
     wald_ci <- rbind(wald_ci,c(max(tmp[1],-1),min(tmp[2],1)))}
   btdf <- as.vector(outer((0:n0x)/n0x,(0:n0y)/n0y,FUN="-"))
   btpdf <- as.vector(estx$btpdf %o% esty$btpdf)
   ord <- order(btdf)
   btdf <- btdf[ord]
   btpdf <- btpdf[ord]
   i <- 0
   cp <- 0
   while(i < (n0x+1)*(n0y+1) && cp < (1-covp)/2)
   {i <- i+1
    cp <- cp + btpdf[i]}
   pct_ci <- btdf[i]
   i <- (n0x+1)*(n0y+1)+1
   cp <- 0
   while(i > 1 && cp < (1-covp)/2)
   {i <- i-1
    cp <- cp + btpdf[i]}
   pct_ci <- c(pct_ci,btdf[i])
   ax <- sum((0:n0x + zq^2/2)/(n0x+zq^2) * estx$btpdf)
   bx <- sum(((0:n0x + zq^2/2)/(n0x+zq^2))^2 * estx$btpdf)-ax^2
   ay <- sum((0:n0y + zq^2/2)/(n0y+zq^2) * esty$btpdf)
   by <- sum(((0:n0y + zq^2/2)/(n0y+zq^2))^2 * esty$btpdf)-ay^2
   zq_ci <- ax-ay+c(-1,1)*zq*sqrt(bx+by)
   zq_ci <- c(max(zq_ci[1],-1),min(zq_ci[2],1))
   b1 <- ifelse(estx$hvar2 != 0, estx$hvar/estx$hvar2, 1)/n0x
   b2 <- ifelse(esty$hvar2 != 0, esty$hvar/esty$hvar2, 1)/n0y
   aa <- 1+zq^2*(b1+b2)
   bb <- -2*diff-zq^2*(b1*(1-2*esty$hss)-b2*(1-2*estx$hss))
   cc <- diff^2-zq^2*(b1*esty$hss*(1-esty$hss)+b2*estx$hss*(1-estx$hss))
   scr_ci <- NULL
   for(i in 1:2) {
     tmp <- (-bb[i]+c(-1,1)*sqrt(bb[i]^2-4*aa*cc[i]))/(2*aa)
     scr_ci <- rbind(scr_ci,c(max(tmp[1],-1),min(tmp[2],1)))}
   
   list(diff=diff, hvar=hvar,
        wald_ci=wald_ci,
        pct_ci=pct_ci,
        zq_ci=zq_ci,
        scr_ci=scr_ci)
}

snsp2mp <- function(mk, n1, s0, covp=0.95, fixsens=TRUE, lbmdis=TRUE)
{if(!lbmdis) mk <- -mk
  if(!fixsens) {
    mk <- -cbind(mk[,-(1:n1)],mk[,1:n1])
    n1 <- dim(mk)[2]-n1
  }

   n0 <- dim(mk)[2]-n1
   est <- .Fortran("eb2mp",as.double(mk),as.integer(n1),as.integer(n0),
                 as.double(s0),diff=double(1),btmn=double(1),btva=double(1),
                 btdist=double(2*n0+1),
                 integer(2*(n1+n0)),integer(2*n1),double((2*n1-1)*2))
                 
   zq <- qnorm(1-(1-covp)/2)
   
   m1est <- snsp1m(mk[1,],n1,s0)
   m2est <- snsp1m(mk[2,],n1,s0)

   diff <- m1est$hss-m2est$hss
   wald_ci <- NULL
   for(i in 1:2) {
     tmp <- diff[i]+c(-1,1)*zq*sqrt(est$btva)
     wald_ci <- rbind(wald_ci,c(max(tmp[1],-1),min(tmp[2],1)))}
   i <- -n0-1
   cp <- 0
   while(i < n0 && cp < (1-covp)/2)
   {i <- i+1
    cp <- cp + est$btdist[n0+1+i]}
   pct_ci <- i/n0
   i <- n0+1
   cp <- 0
   while(i > -n0 && cp < (1-covp)/2)
   {i <- i-1
    cp <- cp + est$btdist[n0+1+i]}
   pct_ci <- c(pct_ci,i/n0)
   zq_ci <- (est$btmn+c(-1,1)*zq*sqrt(est$btva))/(1+zq^2/n0)
   zq_ci <- c(max(zq_ci[1],-1),min(zq_ci[2],1))
   a <- m1est$hvar+m2est$hvar
   a <- ifelse(a != 0, est$btva/a, 1)*zq^2/n0
   b1 <- ifelse(m1est$hvar2 != 0, m1est$hvar/m1est$hvar2, 1)
   b2 <- ifelse(m2est$hvar2 != 0, m2est$hvar/m2est$hvar2, 1)
   aa <- 1+a*(b1+b2)
   bb <- -2*diff-a*(b1*(1-2*m2est$hss)-b2*(1-2*m1est$hss))
   cc <- diff^2-a*(b1*m2est$hss*(1-m2est$hss)+b2*m1est$hss*(1-m1est$hss))
   scr_ci <- NULL
   for(i in 1:2) {
     tmp <- (-bb[i]+c(-1,1)*sqrt(bb[i]^2-4*aa*cc[i]))/(2*aa)
     scr_ci <- rbind(scr_ci,c(max(tmp[1],-1),min(tmp[2],1)))}
   
   list(diff=diff, btmn=est$btmn, btva=est$btva, btdist=est$btdist,
        wald_ci=wald_ci,
        pct_ci=pct_ci,
        zq_ci=zq_ci,
        scr_ci=scr_ci)
}

