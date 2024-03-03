plotirf3 <- function(irf1, irfbands1, irf2, irfbands2, irf3, irfbands3, printvars, printstate, poslegend='topright'){
  
  H <- nrow(irf1)-1
  n <- ncol(irf1)
  lo1 <- irfbands1[,,1]
  up1 <- irfbands1[,,2]
  lo2 <- irfbands2[,,1]
  up2 <- irfbands2[,,2]
  lo3 <- irfbands3[,,1]
  up3 <- irfbands3[,,2]
  numrows <- ceiling(n^(1/2))
  
  par(mfrow=c(numrows,ceiling(n/numrows)), mai = c(0.4, 0.4, 0.4, 0.4))
  for (vv in 1:n){
    maxim <- max(c(irf1[,vv], irf2[,vv], irf3[,vv], lo1[,vv], lo2[,vv], lo3[,vv], up1[,vv], up2[,vv], up3[,vv]))
    maxim <- maxim + 0.1*abs(maxim)
    minim <- min(c(irf1[,vv], irf2[,vv], irf3[,vv], lo1[,vv], lo2[,vv], lo3[,vv], up1[,vv], up2[,vv], up3[,vv]))
    minim <- minim - 0.1*abs(minim)
    plot(0:H, irf1[,vv], 'l', xaxt='n', xaxs='i',
         xlab='', ylab='', main=printvars[vv], ylim=c(minim, maxim))
    grid()
    polygon(c(0:H, H:0), c(lo1[,vv], rev(up1[,vv])), col='gray90', border='gray90')
    abline(h=0, lwd=0.5)
    lines(0:H, irf1[,vv], lwd=2)
    lines(0:H, irf2[,vv], col='blue', lwd=2, lty=2)
    lines(0:H, lo2[,vv], col='blue', lty=2)
    lines(0:H, up2[,vv], col='blue', lty=2)
    lines(0:H, irf3[,vv], col='red', lwd=2, lty=3)
    lines(0:H, lo3[,vv], col='red', lty=3)
    lines(0:H, up3[,vv], col='red', lty=3)
    if (H <= 30){
      axis(1, at=seq(0, H, 4), labels=seq(0, H, 4))
    } else {
      axis(1, at=seq(0, H, 12), labels=seq(0, H, 12)) 
    }
    if (vv == n){
      legend(poslegend, legend=printstate, col=c('black', 'blue','red'), lwd=2, lty=c(1,2,3), box.lty=0)
    }
    box()
  }
  par(mfrow=c(1,1), mar=c(4.6, 4.1, 2.1, 2.1))
}