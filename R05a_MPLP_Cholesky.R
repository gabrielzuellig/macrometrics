
## Preamble
# This script replicates our small VAR where the IRFs are estimated with
# local projections à la Jordà (2005) instead.

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
shocksize <- 0 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear='no')
alpha <- 90 # confidence level

## Load data, generate lags, subset relevant subsample, etc.
load(file='data/data_m_ready.RData')
source('_tbx/supportfct/subset_data.R', echo=T);

## Estimate matrices A, Omega, S, dynamic multipliers
LP <- estimateLP(data, h, p, c_case, exdata, alpha)

## Identification: Estimate VAR to get S
VAR <- estimateVAR(data, p, c_case, exdata);
LP$ident = ident;
LP$shock = matrix(0, n, 1)
LP$shock[shockpos] <- 1;
if (shocksize != 0){  # absolute values, e.g. 25bp = 0.25
  LP$shock <- LP$shock / VAR$S[shockpos, shockpos] * shocksize
}
# structural impulse responses
LP$IRF <- matrix(0, h, n)
LP$IRFbands <- array(0, dim=c(h, n, 2))
for (hh in 1:h){
  LP$IRF[hh,] <- t( LP$Gamma[,,hh] %*% VAR$S %*% LP$shock )
  LP$IRFbands[hh,,1] <- t( LP$Gammalo[,,hh] %*% VAR$S %*% LP$shock )
  LP$IRFbands[hh,,2] <- t( LP$Gammaup[,,hh] %*% VAR$S %*% LP$shock )
}

## Plot impulse responses
plotirf1(LP$IRF, LP$IRFbands, printvars)

