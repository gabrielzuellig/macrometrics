plotirf1 <- function(irf, irfbands, printvars){
  
  H <- nrow(irf)-1
  n <- ncol(irf)
  lo <- irfbands[,,1]
  up <- irfbands[,,2]
  numrows <- ceiling(n^(1/2))

  par(mfrow=c(numrows,ceiling(n/numrows)), mai = c(0.4, 0.4, 0.4, 0.4))
  for (vv in 1:n){
    maxim <- max(c(irf[,vv]), lo[,vv], up[,vv])
    maxim <- maxim + 0.1*abs(maxim)
    minim <- min(c(irf[,vv]), lo[,vv], up[,vv])
    minim <- minim - 0.1*abs(minim)
    plot(0:H, irf[,vv], 'l', lwd=2, xaxt='n', xaxs='i',
         xlab='', ylab='', main=printvars[vv], ylim=c(minim, maxim))
    grid()
    abline(h=0, lwd=0.5)
    lines(0:H, lo[,vv], lty=2)
    lines(0:H, up[,vv], lty=2)
    if (H <= 30){
      axis(1, at=seq(0, H, 4), labels=seq(0, H, 4))
    } else {
      axis(1, at=seq(0, H, 12), labels=seq(0, H, 12)) 
    }
    box()
  }
  par(mfrow=c(1,1), mar=c(4.6, 4.1, 2.1, 2.1))
}


