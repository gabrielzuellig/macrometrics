
## Preamble
# This script reads in monthly time series data from an .xlsx, computes some
# transformations and saves all series in a database called data_m_ready.RData, 
# later to be used for the estimation of time series models.

rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2024/macrometrics/')
source('_tbx/supportfct/compileRfct.R');


## Import
xlsdata <- read.xls('data/data_m_in.xlsx')
data <- xlsdata[,3:ncol(xlsdata)]
labels <- colnames(data)

## Timing and labelling
time <- xlsdata[,1]+(1/12)*(xlsdata[,2]-1)
labelmat <- read.xls('data/data_m_in.xlsx', sheet=2, header=F)
labels_print <- labels
for (vv in 1:length(labels_print)){
  col <- which(labelmat[,1] %in% labels[vv])
  if (length(col) > 0) {
    labels_print[vv] <- as.character(labelmat[col, 2])
  }
}
rm(xlsdata, labelmat, vv, col)

## Data treatment
# 100*log of variables in levels (replace series)
var <- c('IP','CPI','CoreCPI','GDPIHS','PotGDP','PPI')
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  if (length(col) > 0) {
    data[, col] <- 100*log(data[, col])
  }
}

# y/y growth rates (generate new series)
var <- c('IP','CPI','CoreCPI','PPI');
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  if (length(col) > 0) {
    data <- cbind(data, NA)
    data[13:nrow(data), ncol(data)] <- data[13:nrow(data), col] - data[1:(nrow(data)-12), col]
    colnames(data)[ncol(data)] <- paste0(var[vv], 'Gr12')
    labels <- colnames(data)
    labels_print <- c(labels_print, paste0(var[vv],' growth'))
  }
}
labels_print[which(labels %in% 'CPIGr12')] <- 'Inflation'
labels_print[which(labels %in% 'CoreCPIGr12')] <- 'Inflation'
labels_print[which(labels %in% 'PPIGr12')] <- 'Commodity price inflation'

# first difference
var <- c('FFR','FFRshadow','CPI','CoreCPI','GDPgap')
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  if (length(col) > 0) {
    data <- cbind(data, NA)
    data[2:nrow(data), ncol(data)] <- data[2:nrow(data), col] - data[1:(nrow(data)-1), col]
    colnames(data)[ncol(data)] <- paste0('d', var[vv])
    labels <- colnames(data)
    labels_print <- c(labels_print, paste0(labels_print[col],' growth'))
  }
}
labels_print[which(labels %in% 'dFFR')] <- 'Interest rate change'
labels_print[which(labels %in% 'dCPI')] <- 'Inflation'
labels_print[which(labels %in% 'dCoreCPI')] <- 'Inflation'

# moving average
var <- c('dGDPgap')
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  if (length(col) > 0) {
    data <- cbind(data, NA)
    data[, ncol(data)] <- c(stats::filter(data[, col], rep(1/20, 20), sides=1))
    colnames(data)[ncol(data)] <- paste0(var[vv], 'MA20')
    labels <- colnames(data)
    labels_print <- c(labels_print, paste0(labels_print[col], ' Trend'))
  }
}

# change the volatility of the monthly Chicago index implies the same std of
# its quarterly aggregate and the real GDP growth rate
col <- which(labels %in% 'oriChicago')
chicagom <- data[, col]
year <- floor(time)
month <- (time - year)*12+1
quarter <- ifelse(month <= 3.1, 1,
                  ifelse((month >= 3.9 & month <= 6.1), 2,
                  ifelse((month >= 6.9 & month <= 9.1), 3, 4)))
quarter <- year + (quarter/4-0.25)
qvalues <- unique(quarter[year < 1992]) # only need period prior to 1992
chicagoq <- rep(NA, length(qvalues))
for (qq in 1:length(qvalues)){
  chicagoq[qq] <- sum(chicagom[quarter == qvalues[qq]], na.rm=T)
}
dataq <- read.xls('data/data_m_in.xlsx', sheet=4, header=T)
GDPgrowth <- dataq$GDPgrowth[1:length(chicagoq)]
plot(1:length(chicagoq), chicagoq, 'l')
lines(1:length(chicagoq), GDPgrowth, 'l', col='blue')
chicagoq <- chicagoq[2:length(chicagoq)]
GDPgrowth <- GDPgrowth[2:length(chicagoq)]
lm(GDPgrowth ~ 0 + chicagoq)
rm(chicagom, chicagoq, GDPgrowth, month, year, quarter, qq, qvalues)

# monthly output gap: compare own measure to quarterly CBO measure and
# hp-filtered version
dataq <- dataq[1:215, c('GDP','GDPgap')] # real GDP level; CBO output gap
dataq$GDP = 100*log(dataq$GDP);  # 100*log(GDP)
hp <- hpfilter(dataq$GDP, type='lambda', freq=1600)
dataq$GDPhppot <- hp$trend #  hp-filtered potential
dataq$GDPhpgap <- hp$cycle
dataq$idx <- seq(1, nrow(dataq), 1)
dataq3 <- rbind(dataq, dataq, dataq)
dataq3 <- dataq3[order(dataq3$idx), ]
while (nrow(dataq3) < nrow(data)){
  dataq3 <- rbind(dataq3, NA)
}

# Plot descriptive time series
xt00 <- which(time == 2000)
xt <- c(sort(seq(xt00, 0, -60)), seq(xt00+60, length(time), 60))
          
# output gap measures
plot(time, data$GDPgap, 'l', lwd=1.5, ylab='Output gap', xlab='', ylim=c(-12.5, 5))
grid()
lines(time, dataq3$GDPgap, col='blue', lwd=2)
lines(time, dataq3$GDPhpgap, col='red', lty=2, lwd=2)
abline(h=0)        
legend('bottomleft', legend=c('Monthly own measure', 'Quarterly CBO estimate',
                              'Quarterly hp-filteres series (lambda=1600)'),
       col=c('black','blue','red'), lwd=c(1.5,2,2), lty=c(1,1,2), box.lty=0)          
rm(dataq, dataq3, hp)
         
# 'state' graph
s <- data$GDPgap
s2 <- (s - mean(s, na.rm=T)) / sd(s, na.rm=T)
s3 <- (s - quantile(s, 0.05, na.rm=T)) / sd(s, na.rm=T)
Fs <- exp(-1.5*s2) / (1+exp(-1.5*s2))
plot(time, Fs, 'l', lwd=2, ylim=c(0,1.3), xlab='')
grid()  
Fs <- exp(-5*s2) / (1+exp(-5*s2))
lines(time, Fs, col='blue', lwd=2, lty=2)
Fs <- exp(-1.5*s3) / (1+exp(-1.5*s3))
lines(time, Fs, col='red', lwd=2, lty=3)        
legend('topleft', legend=c('gamma = 1.5', 'gamma = 5', 'gamma = 1.5, mu = q10'),
       col=c('black','blue','red'), lwd=2, lty=c(1,2,3), box.lty=0)
rm(s, s2, s3, Fs)
      
# plot all time series
for (vv in 1:ncol(data)){
  plot(time, data[,vv], 'l', lwd=2, xlab='', ylab=labels_print[vv])
  grid()
}

# Export, housekeeping
data_lib <- data
labels_lib <- labels          
printlabels_lib <- labels_print
rm(col, vv, var, data, labels, labels_print, xt, xt00)
save(file='data/data_m_ready.RData', list=c('data_lib','labels_lib','printlabels_lib','time'))

          