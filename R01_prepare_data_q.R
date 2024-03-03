
## Preamble
# This script reads in quarterly time series data from an .xlsx, computes some
# transformations and saves all series in a database called data_q_ready.RData, 
# later to be used for the estimation of time series models.

rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2024/macrometrics/')
source('_tbx/supportfct/compileRfct.R');


## Import
xlsdata <- read.xls('data/data_q_in.xlsx')
data <- xlsdata[,3:ncol(xlsdata)]
labels <- colnames(data)

## Timing and labelling
time <- xlsdata[,1]+(1/4)*(xlsdata[,2]-1)
labelmat <- read.xls('data/data_q_in.xlsx', sheet=2, header=F)
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
var <- c('GDP','PotGP','CPI','CoreCPI','NCONS','GDPDEF','LF','CONS',
         'NINV','INV','GDPPC','HOUR','EMP','THOURS','WAGE','Mortgages','FedAss','Credit');
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  if (length(col) > 0) {
    data[, col] <- 100*log(data[, col])
  }
}

# y/y growth rates (generate new series)
var <- c('GDP','CPI','CoreCPI','MortgGDP');
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  if (length(col) > 0) {
    data <- cbind(data, NA)
    data[5:nrow(data), ncol(data)] <- data[5:nrow(data), col] - data[1:(nrow(data)-4), col]
    colnames(data)[ncol(data)] <- paste0(var[vv], 'Gr4')
    labels <- colnames(data)
    labels_print <- c(labels_print, paste0(var[vv],' growth'))
  }
}
labels_print[which(labels %in% 'CPIGr4')] <- 'Inflation'
labels_print[which(labels %in% 'CoreCPIGr4')] <- 'Inflation'

# first difference
var <- c('GDP','CPI','CoreCPI','FFR','FFRshadow','GDPDEF','GDPPC','CONS','INV','WAGE')
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
labels_print[which(labels %in% 'dFFRshadow')] <- 'Interest rate change'
labels_print[which(labels %in% 'dCPI')] <- 'Inflation'
labels_print[which(labels %in% 'dCoreCPI')] <- 'Inflation'
labels_print[which(labels %in% 'dGDPDEF')] <- 'Inflation'

# hp filtering
var <- c('GDP')
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  if (length(col) > 0) {
    data <- cbind(data, NA)
    notna <- !is.na(data[,col])
    hp <- hpfilter(data[notna, col], type='lambda', freq=1600)
    data[notna, ncol(data)] <- hp$cycle
    colnames(data)[ncol(data)] <- paste0(var[vv], 'hpgap')
    labels <- colnames(data)
    labels_print <- c(labels_print, labels_print[col])
  }
}
rm(hp, notna)

# interest rate in deviation from natural rate
var = c('FFR','FFRshadow');
for (vv in 1:length(var)){
  col <- which(labels %in% var[vv])
  col2 <- which(labels %in% 'Rstar')
  if (length(col) > 0) {
    data <- cbind(data, NA)
    data[, ncol(data)] <- data[, col] - data[, col2]
    colnames(data)[ncol(data)] <- paste0(var[vv], 'dev')
    labels <- colnames(data)
    labels_print <- c(labels_print, labels_print[col])
  }
}
rm(col2)

# Plot descriptive time series
xt00 <- which(time == 2000)
xt <- c(sort(seq(xt00, 0, -20)), seq(xt00+20, length(time), 20))
          
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
save(file='data/data_q_ready.RData', list=c('data_lib','labels_lib','printlabels_lib','time'))

          