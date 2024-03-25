
## Preamble
# This script replicates a small VAR in which a monetary policy shock is 
# identified using Gertler & Karadi high-frequency shocks as an external
# instrument

rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2024/macrometrics/')
source('_tbx/supportfct/compileRfct.R');


## Settings
p <- 12
h <- 48
c_case <- 1
vars <- c('GDPgap','Unemp','CoreCPIGr12','EBP','FFR')
exvars <- NULL
n <- length(vars)
ident <- 'proxy'
proxyvar = ('mpsprGK4');
shockpos <- 5 # position of the variable to instrument
shocksize <- 1 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear='no')
nboot <- 1000
alpha <- 90 # confidence level

## Load data, generate lags, subset relevant subsample, etc.
load(file='data/data_m_ready.RData')
source('_tbx/supportfct/subset_data.R', echo=T);

## Estimate matrices A, Omega, S, dynamic multipliers
VAR <- estimateVAR(data, p, c_case, exdata)
VAR$C <- dyn_multipliers(VAR, h) # not identified

## Identification: Instrumental variable (z)
# prepare
VAR$ident = ident;
VAR$shock = matrix(0, n, 1)
VAR$shock[shockpos] <- 1;
VAR$z <- as.matrix(z[(p+1):nrow(z), ], nrow=nrow(data))
VAR$shockpos <- shockpos
s <- matrix(NA, 1, n)
b <- seq(1,n,1) == shockpos
# save residuals
u_p <- as.matrix(VAR$u[, b], nrow=nrow(data)) # reduced-form residuals to be instrumented
u_q <- VAR$u[, !b] # all other residuals 
# only use periods for which z is not NA 
nona <- !is.na(VAR$z)
u_p <- as.matrix(u_p[nona, ], nrow=sum(nona))
u_q <- u_q[nona,]
z <- cbind(1, as.matrix(VAR$z[nona, ], nrow=sum(nona)))
# 1st stage of IV: u_p on z 
beta1 <- solve(t(z) %*% z) %*% (t(z) %*% u_p) 
u_p_hat <- z %*% beta1  # fitted values of residuals
xi <- u_p - u_p_hat
print('First-stage F-statistic is: ')
Fstat <- ((t(u_p) %*% u_p - t(xi) %*% xi) / (length(beta1))) / ((t(xi) %*% xi) / (length(z)-length(beta1)))
print(Fstat)
# 2nd stage:
s[!b] <- solve(t(u_p_hat) %*% u_p_hat) %*% (t(u_p_hat) %*% u_q) 
s[b] <- 1
# impulse responses
VAR$IRF <- matrix(0, h, n)
for (hh in 1:h){
  VAR$IRF[hh,] <- t( VAR$C[,,hh] %*% t(s))
}                                            
                                             
## Bootstrap
VAR$IRFbands = bootstrapVAR(VAR, nboot, alpha, 'wild')

## Scaling and cleaning up
# so far, s gives us the vector of relative responses. We want absolute values for a 1sd
# shock. To get that, follow Piffer (2020):
if (shocksize == 0){
  sigma11 <- VAR$Omega[b,b]
  sigma12 <- matrix(VAR$Omega[b,!b], nrow=1)
  sigma21 <- matrix(VAR$Omega[!b,b], ncol=1)
  sigma22 <- VAR$Omega[!b,!b]
  mu <- matrix(s[!b], ncol=1)
  Gamma <- sigma22 + mu %*% sigma11 %*% t(mu) - sigma21 %*% t(mu) - mu %*% t(sigma21)
  b11b11prime <- sigma11 - t(sigma21 - mu %*% sigma11) %*% solve(Gamma) %*% (sigma21 - mu %*% sigma11)
  b11 <- t(chol(b11b11prime))
  VAR$shock[shockpos] <- b11
  rm(sigma11, sigma12, sigma21, sigma22, mu, Gamma, b11b11prime)
} else {
  b11 <- shocksize
}
s <- b11 %*% s
VAR$IRF <- c(b11) * VAR$IRF
VAR$IRFbands <- c(b11) * VAR$IRFbands
VAR$s <- c(b11) * VAR$s
VAR$Fstat <- Fstat 
rm(b, xi, u_p, u_q, nona, beta1, u_p_hat, Fstat, s, b11)

## Plot impulse responses
plotirf1(VAR$IRF, VAR$IRFbands, printvars)
