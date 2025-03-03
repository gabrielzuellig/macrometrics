
## Preamble
# This script replicates a small VAR in which a monetary policy shock is 
# identified using the Cholesky decomposition, similar to Christiano, 
# Eichenbaum & Evans (1999), but with notable differences (e.g. monthly frequency)


rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2025/macrometrics/')
source('_tbx/supportfct/compileRfct.R');


## Settings
p <- 12
h <- 48
c_case <- 1
vars <- c('GDPgap','Unemp','CoreCPIGr12','FFR')
exvars <- NULL
n <- length(vars)
ident <- 'chol'
shockpos <- 4 # position of the shock in the Rest matrix columns
shocksize <- 0 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear='no')
nboot <- 1000
alpha <- 90 # confidence level
ignoreyears <- c(2020, 2021)  # years to ignore in regression (e.g. covid)

## Load data, generate lags, subset relevant subsample, etc.
load(file='data/data_m_ready.RData')
source('_tbx/supportfct/subset_data.R', echo=T);

## Estimate matrices A, Omega, S, dynamic multipliers
VAR <- estimateVAR(data, p, c_case, exdata, timeinreg)
VAR$C <- dyn_multipliers(VAR, h) # not identified

## Identification: Cholesky
# organization
VAR$ident = ident
VAR$shock = matrix(0, n, 1)
VAR$shock[shockpos] <- 1
if (shocksize != 0){  # absolute values, e.g. 25bp = 0.25
  VAR$shock <- VAR$shock / VAR$S[shockpos, shockpos] * shocksize
}
# Cholesky decomposition (already done in estimation function)
print(VAR$S)
eps = VAR$u %*% solve(VAR$S)
VAR$eps = eps[, shockpos]
rm(eps)
# impulse responses
VAR$IRF = matrix(0, nrow=h, ncol=n)
for (hh in 1:h){ 
  VAR$IRF[hh,] = t( VAR$C[,,hh] %*% VAR$S %*% VAR$shock )
}

## Bootstrap
VAR$IRFbands = bootstrapVAR(VAR, nboot, alpha, 'residual')

## Plot impulse responses and save VAR structure
plotirf1(VAR$IRF, VAR$IRFbands, printvars)
