
## Preamble
# This script estimates a nonlinear model of monetary policy whose effects
# depend on whether or not the interest rate in the previous 12 months has
# increased or decreased, inspired by Berger et al (2021)

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
ident <- 'chol'
shockpos <- 4 # position of the shock in the Cholesky ordering
shocksize <- 1 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear='yes',
              logistic='yes',
              interacted='no',
              statevar='dFFR',
              gamma=-200,  # because gamma enters negatively, a double negative indicates positive values = regime 2
              cq=50)
alpha <- 90 # confidence level

## Load data, generate lags, subset relevant subsample, etc.
load(file='data/data_m_ready.RData')
source('_tbx/supportfct/subset_data.R', echo=T);

## Overwrite indicator function with interest level change of past 12 months
scum = state$s;
for (pp in 1:11){
  scum <- scum + c(rep(NA, pp), state$s[1:(length(state$s)-pp)])
}
state$s <- matrix(scum[12:length(scum)], ncol=1) # shorten by first 11 periods (NA's)
state$Fs = matrix(as.numeric(state$s > 0), ncol=1);
data <- data[12:nrow(data),]

## Estimate matrices A, Omega, S, dynamic multipliers
LP <- estimateLPnonlin(data, state$Fs, h, p, c_case, exdata, alpha, NWSE=T)

## Identification: Estimate VAR to get S
VAR <- estimateVAR(data, p, c_case, exdata);
LP$ident = ident;
LP$shock = matrix(0, n, 1)
LP$shock[shockpos] <- 1;
if (shocksize != 0){  # absolute values, e.g. 25bp = 0.25
  LP$shock1 <- LP$shock / VAR$S[shockpos, shockpos] * shocksize
  LP$shock2 <- LP$shock / VAR$S[shockpos, shockpos] * shocksize
}
# structural impulse responses
LP$IRF1 <- matrix(0, h, n)
LP$IRF1bands <- array(0, dim=c(h, n, 2))
LP$IRF2 <- matrix(0, h, n)
LP$IRF2bands <- array(0, dim=c(h, n, 2))
for (hh in 1:h){
  LP$IRF1[hh,] <- t( LP$Gamma1[,,hh] %*% VAR$S %*% LP$shock1 )
  LP$IRF1bands[hh,,1] <- t( LP$Gamma1lo[,,hh] %*% VAR$S %*% LP$shock1 )
  LP$IRF1bands[hh,,2] <- t( LP$Gamma1up[,,hh] %*% VAR$S %*% LP$shock1 )
  LP$IRF2[hh,] <- t( LP$Gamma2[,,hh] %*% VAR$S %*% LP$shock2 )
  LP$IRF2bands[hh,,1] <- t( LP$Gamma2lo[,,hh] %*% VAR$S %*% LP$shock2 )
  LP$IRF2bands[hh,,2] <- t( LP$Gamma2up[,,hh] %*% VAR$S %*% LP$shock2 )
}

## Plot impulse responses
plotirf2(LP$IRF1, LP$IRF1bands, LP$IRF2, LP$IRF2bands, printvars, c('Easing cycle','Tightening cycle'))
