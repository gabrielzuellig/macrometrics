
# select T x N matrix of endogenous variables
data <- data_lib[, vars]
printvars <- NULL
for (i in 1:n){ printvars <- c(printvars, printlabels_lib[which(labels_lib %in% vars[i])])}

# select T x N matrix of exogenous variables
exdata <- data_lib[, exvars]
if (ncol(exdata) == 0){ exdata <- NULL }

# select instrument (if necessary)
z = NULL;
if (ident == 'proxy'){
  z <- as.matrix(data_lib[, proxyvar], nrow=nrow(data_lib))
}

# select state (if necessary)
s = NULL;
if (state$nonlinear == 'yes'){
  if (state$logistic == 'yes'){
    s <- as.matrix(data_lib[, state$statevar], nrow=nrow(data_lib))
  }
}

# select interaction (if necessary)
if (state$nonlinear == 'yes'){
  if (state$interacted == 'yes'){
    if (sum(labels_lib %in% c(state$shockvar, state$statevar)) < 2){
      print('Both shockvar and statevar need to be part of the vector of endogenous variables.')
    } else {
      exdata <- matrix(data[, state$shockvar] * data[, state$statevar], ncol=1)
      state$shockpos <- which(vars %in% state$shockvar)
      state$statepos <- which(vars %in% state$statevar)
    }
  }
}

# ignore non-NA columns
nona <- rowSums(is.na(data)) == 0
if (!is.null(exdata)){
  nona[rowSums(is.na(exdata)) > 0] <- FALSE
}
if (!is.null(s)){
  nona[rowSums(is.na(s)) > 0] <- FALSE
}
data <- data[nona, ]
time <- time[nona]
if (!is.null(exdata)){
  exdata <- as.matrix(exdata[nona, ], nrow=sum(nona))
}
if (!is.null(s)){
  s <- as.matrix(s[nona, ], nrow=sum(nona))
}
if (!is.null(z)){
  z <- as.matrix(z[nona, ], nrow=sum(nona))
}

# ignore data before 1965
data <- data[time >= 1965, ]
if (!is.null(exdata)){
  exdata <- as.matrix(exdata[time >= 1965, ], nrow=nrow(data))
}
if (!is.null(s)){
  s <- as.matrix(s[time >= 1965, ], nrow=nrow(data))
}
if (!is.null(z)){
  z <- as.matrix(z[time >= 1965, ], nrow=nrow(data))
}
time <- time[time >= 1965]

# ignore 2020 data
data <- data[time < 2020, ]
if (!is.null(exdata)){
  exdata <- as.matrix(exdata[time < 2020, ], nrow=nrow(data))
}
if (!is.null(s)){
  s <- as.matrix(s[time < 2020, ], nrow=nrow(data))
}
if (!is.null(z)){
  z <- as.matrix(z[time < 2020, ], nrow=nrow(data))
}
time <- time[time < 2020]

# Logistic transformation of s-variable 
if (state$nonlinear == 'yes'){
  if (state$logistic == 'yes'){
    state$s_orig <- s
    state$s_orig_cons <- unname(quantile(s, probs=state$cq/100))
    state$s_orig_std <- sd(s, na.rm=T)
    state$s <- (state$s_orig - state$s_orig_cons) / state$s_orig_std
    state$Fs <- exp(-state$gamma*state$s)/(1+exp(-state$gamma*state$s))
  }
}

# dummy of state
if (state$nonlinear == 'yes'){
  if (state$interacted == 'yes'){
    state$absval <- unname(quantile(data[, state$statepos], probs=state$cq/100))
    state$s <- as.numeric(data[, state$statepos] <= state$absval)
  }
}

# generate plot of actual data that goes into model
frequency <- sum(floor(time) == 2000) # 4 if quarterly, 12 if monthly
xt00 <- which(time == 2000)
if (frequency == 4){
  xt <- c(sort(seq(xt00, 0, -40)), seq(xt00+40, length(time), 40))
} else if (frequency == 12){
  xt <- c(sort(seq(xt00, 0, -120)), seq(xt00+120, length(time), 120))
}
  
numrows <- ceiling(n^(1/2))
par(mfrow=c(numrows,ceiling(n/numrows)), mai = c(0.4, 0.4, 0.4, 0.4))
for (vv in 1:n){
  plot(time, data[,vv], 'l', lwd=2, xlab='', ylab='', main=printvars[vv])
  grid()
  abline(h=0, lwd=0.5)
}
par(mfrow=c(1,1), mar=c(4.6, 4.1, 2.1, 2.1))

if (state$nonlinear == 'yes'){
  if (state$logistic == 'yes'){
    par(mfrow=c(2,2), mai = c(0.4, 0.4, 0.4, 0.4))
    # subplot 1: s vs. time
    plot(time, state$s_orig, 'l', lwd=2, xaxs='i', main='s(t)')
    grid()
    abline(h=0)
    # subplot 2: F(shat) vs. s
    temps <- unname(quantile(state$s_orig, probs=c(seq(0.05, 0.95, 0.01))))
    tempshat <- (temps - state$s_orig_cons) / state$s_orig_std
    tempFs <- exp(-state$gamma * tempshat) / (1 + exp(-state$gamma * tempshat))
    plot(temps, tempFs, 'l', lwd=2, ylim=c(0, 1), main='F(s)')
    grid()
    points(state$s_orig, state$Fs, pch=4, cex=0.5, col='blue')
    rm(temps, tempshat, tempFs)
    # subplot 3: F(s) vs. time
    plot(time, state$Fs, 'l', lwd=2, xaxs='i', main='F(t)', ylim=c(0,1))
    grid()
  }
  if (state$interacted == 'yes'){
    par(mfrow=c(2,1), mai = c(0.4, 0.4, 0.4, 0.4))
    # subplot 1:
    plot(time, data[,state$shockpos], 'l', lwd=2, main='Interacted variable 1: Shock', xaxs='i')
    grid()
    abline(h=0)
    # subplot 2:
    plot(time, data[,state$statepos], 'l', lwd=2, main='Interacted variable 2: State', xaxs='i')
    polygon(c(time, rev(time)), 
            c(1.1*min(data[, state$statepos])*state$s, rev(1.1*max(data[, state$statepos])*state$s)), 
            col='gray90', border='gray90')
    lines(time, data[,state$statepos], 'l', lwd=2)
    grid()
    abline(h=0)
    legend('topleft',legend='Regime 2', fill='gray90', cex=0.75, box.lty=0)
    box()
  }
}
par(mfrow=c(1,1), mar=c(4.6, 4.1, 2.1, 2.1))

rm(vv, nona, xt00, xt, numrows, s)
