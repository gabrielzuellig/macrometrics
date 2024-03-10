
## Preamble
# This script replicates a small VAR in which a monetary policy shock is 
# identified using sign restrictions, similar to Uhlig (2005), but with 
# notable differences (e.g. restriction on GDP)

rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2024/macrometrics/')
source('_tbx/supportfct/compileRfct.R');


## Settings
p <- 12
h <- 48
c_case <- 1
vars <- c('GDPgap','Unemp','CoreCPIGr12','FFR')
exvars <- NULL
n <- length(vars)
ident <- 'sign'
shockpos <- 1 # position of the shock in the Rest matrix columns
shocksize <- 0 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear='no')
nboot <- 1000
alpha <- 90 # confidence level

## Load data, generate lags, subset relevant subsample, etc.
load(file='data/data_m_ready.RData')
source('_tbx/supportfct/subset_data.R', echo=T);

## Estimate matrices A, Omega, S, dynamic multipliers
VAR <- estimateVAR(data, p, c_case, exdata)
VAR$C <- dyn_multipliers(VAR, h) # not identified

## Identification: Sign restrictions
VAR$ident <- ident
# define restrictions // each column is one of q identified shocks 
# and each of n rows contains -1 and 1 for negative/positive restrictions
# of the row-variable; 
Rest <- matrix(c(c(-1, NA, -1, 1), # monetary policy shock
               c(1, NA, 1, 1), # demand shock
               c(-1, NA, 1, 1), # supply shock
               c(NA,NA,NA,NA)), nrow=n, byrow=F)
print(Rest)
if (ncol(Rest) != n){
  print('Error. Number of rows in Rest matrix has to be equal to n.')
}
Rest <- Rest[, colSums(!is.na(Rest)) > 0]  # select only columns that actually contain identifying restrictions
Start <- matrix(c(c(5, NA, 5, 0),
                  c(5, NA, 5, 0),
                  c(5, NA, 5, 0),
                  c(NA, NA, NA, NA)), nrow=n, byrow=F) # first period for binding restriction
Start <- Start[1:nrow(Rest), 1:ncol(Rest)]
End <- 5*matrix(1, nrow=nrow(Rest), ncol=ncol(Rest)) # final period for binding restriction
End <- End[1:nrow(Rest), 1:ncol(Rest)]
nshocks <- ncol(Rest)
Rest_check <- rep(0, nshocks)
perc_accepted <- rep(NA, nboot)
Qb <- array(NA, dim=c(n,n,nboot))

# while-loop to look for nboot different Q's that satisfy the restrictions 
accepteddraws <- 0
attempts <- 1
print('Progress...')
pb <- txtProgressBar(min = 0, max = nboot, initial = 0, style=3) 
while (accepteddraws < nboot){
  
  # draw a rotation of Q
  draw <- matrix(rnorm(n^2), ncol=n)
  QR <- qr(draw)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  for (nn in 1:n){
    if (R[nn,nn]<0){
      Q[,nn] = -Q[,nn]
    }
  }
  
  ii <- 0
  while (ii<nshocks){ # loop over each of the shocks (columns in 'Rest')
    ii <- ii+1
    shock <- matrix(0, n, 1)
    shock[ii,1] <- 1
    IRFcandidate <- matrix(NA, nrow=max(End)+1, ncol=n)
    for (hh in 1:(max(End)+1)){
      IRFcandidate[hh,] <- t( VAR$C[,,hh] %*% VAR$S %*% Q %*% shock )
    }
    
    junk <- matrix(seq(1,n,1), nrow=n)
    Which_Variab <- matrix(junk[!is.na(Rest[,ii])], ncol=1)
    Which_Sign <- matrix(Rest[!is.na(Rest[,ii]), ii], ncol=1)
    Which_Start <- matrix(Start[!is.na(Rest[,ii]), ii], ncol=1)
    Which_End <- matrix(End[!is.na(Rest[,ii]), ii], ncol=1)
    
    jj <- 1
    while (jj <= nrow(Which_Variab)){
      a <- 1
      rev_cond <- 1; # check all restrictions the way the are written
      a_temp <- testIRFsign(IRFcandidate[(Which_Start[jj,1]+1):(Which_End[jj,1]+1), Which_Variab[jj, 1]], Which_Sign[jj, 1], rev_cond)
      a <- a*a_temp;
      if (a_temp == 1){
        jj <- jj+1
      } else {
        jj <- nrow(Which_Variab)+100
      }
    }
    
    if (a == 0){ # try to invert all the signs 
      jj <- 1
      while (jj <= nrow(Which_Variab)){
        a <- 1
        rev_cond <- 0; # check all restrictions the way the are written 
        a_temp <- testIRFsign(IRFcandidate[(Which_Start[jj,1]+1):(Which_End[jj,1]+1), Which_Variab[jj, 1]], Which_Sign[jj, 1], rev_cond)
        a <- a * a_temp 
        if (a == 1){
          jj <- jj+1
        } else {
          jj <- nrow(Which_Variab)+100
        }
      }
      if (a == 1){
        Q[,ii] <- -Q[,ii]
      }
    }
    
    if (a == 0){
      Rest_check[ii] <- 0
      ii <- nshocks+100;   #The restrictions are not satisfied.  Go the beginning to redraw.
    } else {
      Rest_check[ii] <- 1
    }

  }
  
  if (mean(Rest_check) == 1){
    # save accepted Q's and show progress bar
    setTxtProgressBar(pb,accepteddraws)
    accepteddraws <- accepteddraws + 1;
    perc_accepted[accepteddraws] = (accepteddraws / attempts)*100;
    Qb[,,accepteddraws] = Q;
  }
  attempts <- attempts + 1
  
}

# figure: share of accepted draws (in %)
plot(1:nboot, perc_accepted, 'l', lwd=2, main='Percent of draws accepted', xlab='Draw', ylab='')
grid()
rm(accepteddraws, attempts, Rest_check, junk, Which_End, Which_Sign,
   Which_Start, Which_Variab, ii, jj, a, a_temp, draw, R, rev_cond,
   perc_accepted, nn, IRFcandidate)

# Produce IRFs with accepted Q's for each shock/column in 'Rest'
for (s in 1:nshocks){
  
  ## Impulse response functions
  # Constructions
  IRFb <- array(NA, dim=c(h,n,nboot))
  shock <- matrix(0, n, 1)
  shock[s,1] <- 1
  for (bb in 1:nboot){
    P <- VAR$S %*% Qb[,,bb]
    for (hh in 1:h){
      IRFb[hh,,bb] <- t( VAR$C[,,hh] %*% P %*% shock )
    }
  }
  IRF <- apply(IRFb, c(1,2), mean)
  IRFbands <- array(NA, dim=c(h,n,2))
  lo <- (100-alpha)/2
  up <- 100 - lo
  for (i in 1:n){
    IRFbands[,i,] <- t(apply(IRFb[,i,], 1, FUN=function(x) quantile(x, probs=c(lo/100,up/100))))  
  }
  
  if (s == shockpos){
    VAR$shock <- shock
    VAR$Rest <- Rest
    VAR$Start <- Start
    VAR$End <- End
    VAR$Q <- apply(Qb, 3, mean)
    VAR$IRF <- IRF
    VAR$IRFbands <- IRFbands
  }
  
  # Plot
  plotirf1(IRF, IRFbands, printvars)
  
}

# Organize
rm(s, bb, hh, Rest, Start, End, shock, shockpos, Qb, IRFb, lo, up, P, IRF, Q, pb, QR, IRFbands, nshocks)


