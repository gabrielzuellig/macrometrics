bootstrapVAR <- function(VAR, nboot, alpha=90, method){
  
  # unpacking
  t <- nrow(VAR$u)
  n <- VAR$n
  u <- VAR$u
  A <- VAR$A
  c_case <- VAR$c_case
  X <- VAR$X
  Xex <- VAR$Xex
  if (!is.null(Xex)){
    exdata <- VAR$exdata[(p+1):nrow(VAR$exdata),]
  } else {
    exdata <- NULL
  }
  p <- VAR$p
  IRFb <- array(NA, dim=c(nrow(VAR$IRF), ncol(VAR$IRF), nboot))
  h <- nrow(VAR$IRF)
  
  # define upper and lower thresholds
  lo <- (100-alpha)/2
  up <- 100 - lo
  
  # combine constant, endogenous and exogenous variables
  if (c_case == 0){
    c <- NULL
  } else if (c_case == 1){
    c <- matrix(1, nrow=nrow(X), ncol=1)
  } else if (c_case == 2){
    c <- cbind(matrix(1, nrow=nrow(X), ncol=1), matrix(1:nrow(X), nrow=nrow(X), ncol=1))
  }
  endogvars <- c_case + seq(1,ncol(X),1)
  X <- cbind(c, X, Xex)
  if ((!is.null(Xex)) & (method == 'residual')){
    method <- 'wild'
    print('Changing method to Wild bootstrapping (due to exogenous variables).')
  }
  if ((VAR$ident == 'proxy') & (method == 'residual')){
    method <- 'wild'
    print('Changing method to Wild bootstrapping (due to external instrument).')
  }
  
  print('Progress...')
  pb <- txtProgressBar(min = 0, max = nboot, initial = 0, style=3) 
  for (bb in 1:nboot){  # 'nboot' bootstraps
    
    # generate pseudo-disturbance 'ub'
    if (method == 'residual'){
      # residual bootstrapping assumes E(y|Xb) = X*b and u are iid (i.e.
      # homoskedasticity => Draw from the empirical distribution of u
      segment <- c(1:t)/t
      ub <- matrix(0, nrow=ceiling(1.25*t), ncol=n)
      for (i in 1:nrow(ub)){
        draw <- runif(1,0,1); # uniform distribution
        ub[i,] <- u[min(which(segment>=draw)), ]
      }
    } else if (method == 'wild'){
      # allows for heteroskedasticity
      fu <- matrix(1 - 2*(runif(t,0,1)>0.5), ncol=1)
      ub <- matrix(0, nrow=nrow(VAR$u), ncol=n)
      for (i in 1:n){
        ub[,i] <- u[,i]*fu # flip sign of randomly selected 50%
      }
    }
    
    # generate pseudo-sample based on drawn u's
    Yb <- matrix(0, nrow=nrow(ub), ncol=n)
    r <- X[1,]
    for (i in 1:nrow(ub)){
      Yb[i,] <- r %*% A + ub[i,]
      r <- c(Yb[i,], r[endogvars[1:(length(endogvars)-n)]])
      if (c_case == 1){
        r <- c(1, r)
      } else if (c_case == 2){
        r <- c(1, i, r)
      }
      if (!is.null(Xex)){
        r <- c(r, Xex[i,])
      }
    }
    data <- Yb[(nrow(Yb)-nrow(VAR$u)+1):(nrow(Yb)),]
    if (VAR$ident == 'proxy'){
      zb <- VAR$z * fu
      zb <- as.matrix(zb[(p+1):nrow(zb), ], nrow=nrow(data)-p)
    }
    
    # estimate VAR and IRF
    VARb <- estimateVAR(data, p, c_case, exdata)
    VARb$C <- dyn_multipliers(VARb, h)
    
    # identify shock
    if (VAR$ident == 'chol'){
      for (hh in 1:h){
        IRFb[hh,,bb] <- t(VARb$C[,,hh] %*% VARb$S %*% VAR$shock)
      }
    } else if (VAR$ident == 'proxy'){
      b <- seq(1,n,1) == VAR$shockpos
      u_p <- as.matrix(VARb$u[, b], nrow=nrow(data)) # reduced-form residuals to be instrumented
      u_q <- VARb$u[, !b] # all other residuals 
      nona <- !is.na(zb)
      u_p <- as.matrix(u_p[nona, ], nrow=sum(nona))
      u_q <- u_q[nona,]
      zb <- cbind(1, as.matrix(zb[nona, ], nrow=sum(nona)))
      # 1st-stage of IV: u_p on z
      beta1 <- solve(t(zb) %*% zb) %*% (t(zb) %*% u_p) 
      u_p_hat <- zb %*% beta1      
      # 2nd-stage
      sb <- matrix(NA, nrow=1, ncol=n)
      sb[!b] <- solve(t(u_p_hat) %*% u_p_hat) %*% (t(u_p_hat) %*% u_q) 
      sb[b] <- 1  # normalize s at shockpos to 1
      # impulse responses
      for (hh in 1:h){
        IRFb[hh,,bb] <- t( VARb$C[,,hh] %*% t(sb))
      }  
    }
    
    # progress bar
    setTxtProgressBar(pb,bb)
  }
  
  # retrieve intervals
  IRFbands <- array(NA, dim=c(h,n,2))
  for (i in 1:n){
    IRFbands[,i,] <- t(apply(IRFb[,i,], 1, FUN=function(x) quantile(x, probs=c(lo/100,up/100))))  
  }
  
  return(IRFbands)
 
}