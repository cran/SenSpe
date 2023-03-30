snsp1m <- function(mk, n1, s0, covp=0.95, fixsens=TRUE, lbmdis=TRUE)
{if(!lbmdis) mk <- -mk
  if(!fixsens) {
    mk <- -c(mk[-(1:n1)],mk[1:n1])
    n1 <- length(mk)-n1
  }

 if(s0 <= 0 || s0 >= 1) stop("Control level out of permissible bound (0,1)")
 if(covp <= 0 || covp >= 1) stop("CI level out of permissible bound (0,1)")

 n <- length(mk)
 n0 <- n-n1
 mk <- c(sort(mk[1:n1]),sort(mk[-(1:n1)]))
 zq <- qnorm(1-(1-covp)/2)

 r <- floor(n1*(1-s0))+1
 btp1 <- rep(0,n1)
 for(i in 1:n1) btp1[i] <- pbeta(i/n1,r,n1-r+1)
 btp1[2:n1] <- btp1[2:n1] - btp1[1:(n1-1)]

  ind <- c(rep(1,n1),rep(0,n0))
  ord <- order(mk,1-ind)
  ind <- ind[ord]
  crtind <- ind[1]
  str <- NULL
  cnt <- 0
  for(i in 1:n)
  {if(ind[i] == crtind) cnt <- cnt+1
   else {str <- c(str,cnt)
         cnt <- 1
         crtind <- ind[i]}
   if(i == n) str <- c(str,cnt)}
  btpdf <- rep(0,n0+1)  
  m <- length(str)
  nsp <- 0
  nsn <- 0
  crtind <- ind[1]
  bt_mn_spec <- NULL
  for(i in 1:m)
  {if(crtind == 0) nsp <- nsp+str[i]
   else {pcs <- sum(btp1[nsn+(1:str[i])])
         bt_mn_spec <- rbind(bt_mn_spec,c(pcs,nsp/n0))
         if(nsp == n0 || nsp == 0) btpdf[nsp+1] <- btpdf[nsp+1]+pcs
         else btpdf <- btpdf+dbinom(0:n0,n0,nsp/n0)*pcs
         nsn <- nsn+str[i]}
   crtind <- 1-crtind}

   threshold <- mk[r]
   hss <- sum(mk[(n1+1):n]<threshold)/n0
   bmn <- sum(bt_mn_spec[,1]*bt_mn_spec[,2])
   hvar1 <- sum(bt_mn_spec[,1]*bt_mn_spec[,2]^2)-bmn^2
   hvar2 <- sum(bt_mn_spec[,1]*bt_mn_spec[,2]*(1-bt_mn_spec[,2]))/n0
   hvar <- hvar1+hvar2
   scale <- ifelse(hvar2 != 0, hvar/hvar2, 1)*zq^2/n0
   cis <- waldscr(hss,hvar,zq,scale)
   wald_ci <- cis$wald_ci
   scr_ci <- cis$scr_ci
   i <- 0
   cp <- 0
   while(i < n0+1 && cp < (1-covp)/2)
   {i <- i+1
    cp <- cp + btpdf[i]}
   pct_ci <- (i-1)/n0
   i <- n0+2
   cp <- 0
   while(i > 1 && cp < (1-covp)/2)
   {i <- i-1
    cp <- cp + btpdf[i]}
   pct_ci <- c(pct_ci,(i-1)/n0)
   a <- sum((0:n0 + zq^2/2)/(n0+zq^2) * btpdf)
   b <- sqrt(sum(((0:n0 + zq^2/2)/(n0+zq^2))^2 * btpdf)-a^2)
   zq_ci <- a+c(-1,1)*zq*b
   zq_ci <- c(max(zq_ci[1],0),min(zq_ci[2],1))

   rd <- (n1+1)*(1-s0)
   r <- c(max(floor(rd),1),min(ceiling(rd),n1))
   if(r[1] == r[2])
   {hss <- hss * c(1,1)
    wald_ci <- rbind(wald_ci,wald_ci)
    scr_ci <- rbind(scr_ci,scr_ci)} else
   {hss <- c(hss,(sum(mk[(n1+1):n]<mk[r[1]]) * (r[2]-rd)
            +sum(mk[(n1+1):n]<mk[r[2]]) * (rd-r[1]))/n0)
    cis <- waldscr(hss[2],hvar,zq,scale)
    wald_ci <- rbind(wald_ci,cis$wald_ci)
    scr_ci <- rbind(scr_ci,cis$scr_ci)}

   rownames(wald_ci) <- NULL
   rownames(scr_ci) <- NULL
   list(threshold=threshold,
        hss=hss,
        hvar1=hvar1,
        hvar2=hvar2,
        hvar=hvar,
        btpdf=btpdf,
        wald_ci=wald_ci,
        pct_ci=pct_ci,
        scr_ci=scr_ci,
        zq_ci=zq_ci)
}

waldscr <- function(hss,hvar,zq,scale)
{wald_ci <- hss+c(-1,1)*zq*sqrt(hvar)
 wald_ci <- c(max(wald_ci[1],0),min(wald_ci[2],1))
 tmp <- sqrt((hss*(1-hss)+scale/4)*scale)
 tmp <- hss+scale/2+c(-1,1)*tmp
 scr_ci <- tmp/(1+scale)

 list(wald_ci=wald_ci,scr_ci=scr_ci)
}
