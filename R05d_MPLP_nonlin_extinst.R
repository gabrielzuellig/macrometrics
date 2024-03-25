
## Preamble
# This script estimates local projections with Gertler-Karadi shocks in
# expansions vs. recessions.

rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2024/macrometrics/')
source('_tbx/supportfct/compileRfct.R');


## Settings
p <- 6  # easily overfitted
h <- 48
c_case <- 1
vars <- c('GDPgap','Unemp','CoreCPIGr12','FFR')
exvars <- NULL
n <- length(vars)
ident <- 'proxy'
proxyvar <- c('mpsprGK4')
shockpos <- 4 # position of the shock in the Cholesky ordering
shocksize <- 0 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear='yes',
              logistic='yes',
              interacted='no',
              statevar='Unemp',
              gamma=-5, # because gamma enters negatively, a double negative indicates positive values = regime 2
              cq=60)
alpha <- 90 # confidence level

## Load data, generate lags, subset relevant subsample, etc.
load(file='data/data_m_ready.RData')
source('_tbx/supportfct/subset_data.R', echo=T);

## Estimate matrices A, Omega, S, dynamic multipliers
LP <- estimateLPznonlin(data, z, state$Fs, h, p, c_case, exdata, alpha, NWSE=F)

## Plot impulse responses
plotirf2(LP$gamma1, LP$gamma1bands, LP$gamma2, LP$gamma2bands, printvars, c('Expansion','Recession'))
