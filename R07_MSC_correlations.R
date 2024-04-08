
## Preamble
# This script load and cleans the Michigan Survey of Consumers dataset and runs
# a few discriptive regressions.

rm(list = ls())

# before running, set working directory to source file location and comment out the following line
setwd('~/Dropbox/Teaching/Basel2024/macrometrics/')


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


## Recode survey variables
# Define function
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
# u_p4: unemployment outlook (categorical value made to numeric)
msc <- myrecode('UNEMP', 'u_p4', c(1, 3, 5, 8, 9), c(1, 0, -1, NA, NA), data=msc)
hist(as.numeric(msc$u_p4), freq=F,
     xlab='Unemployment outlook', main='')
grid()
box()
# ec_p4: economic outlook (categorical value made to numeric)
msc <- myrecode('BUS12', 'ec_p4', c(1, 2, 3, 4, 5, 8, 9), c(1, .5, 0, -.5, -1, NA, NA), data=msc)
hist(as.numeric(msc$ec_p4), freq=F,
     xlab='Economic outlook', main='')
grid()
box()
# pi_p4: inflation expectations (numeric value)
summary(msc$PX1)
msc$pi_p4 <- ifelse(abs(msc$PX1) <= 95, msc$PX1, NA) #Abs > 95 are values like don't know or invalid answers 
quantile(msc$pi_p4, probs=c(0, 0.01, 0.02, seq(0.05, 0.95, 0.05), 0.98, 0.99, 1), na.rm=T)
hist(msc$pi_p4, breaks=40, freq=F,
     xlab='Inflation expectations', main='')
grid()
box()


## Keep only the variables relevant (for now)
msc <- msc[, c('CASEID', 'ID', 'date', 'hasprevious', 'prevdate', 'IDPREV', 'u_p4', 'ec_p4', 'pi_p4')]


## Average expected inflation over time
avg_pi_p4 <- data.table(msc)
avg_pi_p4 <- avg_pi_p4[, j=list(pi_p4=mean(pi_p4, na.rm=T)),
                       by=list(date)]
plot(avg_pi_p4$date, avg_pi_p4$pi_p4, 'l', lwd=2, ylim=c(0, 13),
     xlab='', ylab='Average expected inflation') 
grid() 


## Regress inflation outlook on unemployment
# Plain
reg1a <- lm(pi_p4 ~ u_p4, data=msc)
summary(reg1a)  # positive sign: high unemployment <= higher inflation expectations
# Include time fixed effect 
library(lfe)
reg2a <- felm(pi_p4 ~ u_p4 | date, data=msc)
library(stargazer)
stargazer(reg1a, reg2a, type='text')


## Match previous interview with same person and calculate revisions
mscprev <- msc
msc$by <- paste0(round(msc$prevdate, 2), '-', msc$IDPREV)
msc$by[msc$by == 'NANA'] <- NA
mscprev$by = paste0(round(mscprev$date, 2), '-', mscprev$ID)
msc <- merge(msc, mscprev[, c('by','u_p4','ec_p4','pi_p4')], by='by',
             all.x=TRUE, suffixes = c('','.l1'))
msc <- subset(msc, select=-c(by, CASEID, ID, IDPREV))
msc <- msc[order(msc$date), ]
rm(mscprev)
# Transform wide to long
msc2 <- msc[msc$hasprevious == TRUE, ]
msc2$hhid = c(1:nrow(msc2))
colnames(msc2)[which(colnames(msc2) == 'prevdate')] <- 'date.l1'
msc2 <- reshape(msc2, direction='long', 
                varying=c('ec_p4', 'u_p4', 'pi_p4', 'ec_p4.l1', 'u_p4.l1', 'pi_p4.l1'), 
                timevar='intid',
                v.names=c('ec_p4', 'u_p4', 'pi_p4'),
                idvar='hhid')
msc2 <- msc2[order(msc2$hhid, msc2$date), ]


## Regress inflation outlook on unemployment with time and person fixed effect
reg3a <- felm(pi_p4 ~ u_p4 | date + hhid, data=msc2)
stargazer(reg1a, reg2a, reg3a, type='text')


## Repeat everything with economic outlook
reg1b <- lm(pi_p4 ~ ec_p4, data=msc)
reg2b <- felm(pi_p4 ~ ec_p4 | date, data=msc)
reg3b <- felm(pi_p4 ~ ec_p4 | date + hhid, data=msc2)
stargazer(reg1b, reg2b, reg3b, type='text')
