
## Preamble
# This script load and cleans the Michigan Survey of Consumers dataset and runs
# a few discriptive regressions.

rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2025/macrometrics/')


## Load raw data 
msc <- read.csv('./data/data_msc_in.csv')
colnames(msc)

## Timeset the data
library(data.table)
# Get numeric variable of date
msc$date <- as.Date(paste0(msc$YYYYMM, '01'), format='%Y%m%d')
msc$dateqtr <- year(msc$date) + (quarter(msc$date)-1)/4  # numeric variable used throughout
msc$date <- year(msc$date) + (month(msc$date)-1)/12
msc$prevdate <- as.Date(paste0(msc$DATEPR, '01'), format='%Y%m%d')
msc$prevdate <- year(msc$prevdate) + (month(msc$prevdate)-1)/12
msc$interval <- msc$date - msc$prevdate   # time intervals between 2 interviews, if any
table(msc$interval, useNA='always')   # almost exclusively 6 months
# Mark re-interviews
msc$hasprevious <- !is.na(msc$IDPREV)
table(msc$hasprevious)  # around 1/3 are re-interviews
# Check for duplicates
msc$dupl <- duplicated(msc[, c('CASEID','YYYY','YYYYQ')]) | duplicated(msc[, c('CASEID','YYYY','YYYYQ')], fromLast=TRUE)
summary(msc$dupl)  #none


## Show number of interviews by month
ninterviews <- data.table(msc)
ninterviews <- ninterviews[, j=list(nint=length(CASEID),
                                    reint=sum(hasprevious)), by=list(date)]
plot(ninterviews$date, ninterviews$nint, 'l', lwd=2, ylim=c(0, 1500),
     xlab='', ylab='Number of interviews') 
grid() 
lines(ninterviews$date, ninterviews$nint, lwd=2, lty=1)
lines(ninterviews$date, ninterviews$reint, lwd=2, lty=2)
rm(ninterviews)


## Prepare variables 
# Consumption: index for consumer sentiment (based on several other variables)
summary(msc$ICS)
msc$consumption <- msc$ICS
# Inflexp: inflation expectations (numeric value)
summary(msc$PX1)
msc$inflexp <- ifelse(abs(msc$PX1) < 95, msc$PX1, NA) #Abs > 95 are values like don't know or invalid answers 
quantile(msc$inflexp, probs=c(0, 0.01, 0.02, seq(0.05, 0.95, 0.05), 0.98, 0.99, 1), na.rm=T)
hist(msc$inflexp, breaks=40, freq=F,
     xlab='Inflation expectations', main='')
grid()
box()
# Define function to recode (not needed, just to show how to treat categorical variables)
myrecode <- function(fromvar, tovar, fromval, toval, data){
  
  oldnames <- colnames(data)
  old <- data[, fromvar]
  new <- rep(NA, length(old))
  for (s in 1:length(fromval)){
    new[old == fromval[s]] <- toval[s] 
  }
  
  data <- cbind(data, new)
  colnames(data) <- c(oldnames, tovar)
  return(data)
  
}
# Unempexp: unemployment outlook (categorical value made to numeric)
msc <- myrecode('UNEMP', 'unempexp', c(1, 3, 5, 8, 9), c(1, 0, -1, NA, NA), data=msc)
hist(as.numeric(msc$unempexp), freq=F,
     xlab='Unemployment outlook', main='')
grid()
box()
# Keep only the variables relevant (for now)
msc <- msc[, c('CASEID', 'ID', 'date', 'hasprevious', 'prevdate', 'IDPREV', 'consumption', 'inflexp')]


## Average expected inflation over time
avg_pi_p4 <- data.table(msc)
avg_pi_p4 <- avg_pi_p4[, j=list(inflexp=mean(inflexp, na.rm=T)),
                       by=list(date)]
plot(avg_pi_p4$date, avg_pi_p4$inflexp, 'l', lwd=2, ylim=c(0, 13),
     xlab='', ylab='Average expected inflation') 
grid() 


## Regress consumption attitude on inflation expectations
# Plain
reg1 <- lm(consumption ~ inflexp, data=msc)
summary(reg1)  # negative sign: high expected inflation <=> low appetite for consumption
# Include time fixed effect 
library(lfe)
reg2 <- felm(consumption ~ inflexp | date, data=msc)
library(stargazer)
stargazer(reg1, reg2, type='text')  # still negative

