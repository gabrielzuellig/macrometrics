
/*
# This script load and cleans the Michigan Survey of Consumers dataset and runs
# a few discriptive regressions.
*/

*** Program setup
capture log close
clear all
macro drop _all
set more off
set segmentsize 2g
program drop _all
global path = "/Users/gabrielzullig/Dropbox/Teaching/Basel2025/macrometrics/"


** Load raw data 
import delim using $path/data/data_msc_in.csv, varn(1) clear
describe

** Timeset the data 
gen mon = yyyymm - yyyy*100
gen date = ym(yyyy, mon)
format date %tm
gen dateqtr = qofd(dofm(date))
gen yyyy_pr = substr(datepr, 1, 4)
destring yyyy_pr, replace
gen mon_pr = substr(datepr, 5, 6)
destring mon_pr, replace
gen prevdate = ym(yyyy_pr, mon_pr)
format prevdate %tm
gen interval = date - prevdate   // time intervals between 2 interviews, if any
tab interval, m   // almost exclusively 6 months
// mark re-interviews 
destring idprev, replace
gen hasprevious = idprev != .
tab hasprevious  // around 1/3 are re-interviews 


** Show number of interviews by month 
preserve // everything between 'preserve' and 'restore' should be run combined so original data is not erased from workspace

collapse (count) nint=caseid (sum) reint=hasprevious, by(date)
graph twoway (line nint date) ///
	(line reint date, lp(dash)), ///
	legend(pos(6) order(1 "Interviews" 2 "of which re-interviews"))

restore 


** Prepare variables
* consumption: index for consumer sentiment (based on several other variables)
sum ics   // numeric values
rename ics consumption
* inflexp: inflation expectations
destring px1, replace
sum px1, det
gen inflexp = px1 if abs(px1) < 95 // abs > 95 are value like don't know or invalid answers
sum inflexp, det
hist inflexp
* Keep only the relevant variables
keep caseid id date hasprevious prevdate idprev consumption inflexp


** Average expected inflation over time 
preserve 

collapse (mean) inflexp, by(date)
graph twoway (line inflexp date), xtitle("") ytitle("Average expected inflation")

restore 


** Regress consumption attitude on inflation expectations
* Download packages if necessary (only required once)
ssc install estout
ssc install ftools
ssc install reghdfe
* Plain
eststo reg1: reg consumption inflexp  // negative sign: high expected inflation <=> low appetite for consumption
* Include time fixed effect
eststo reg2: reg consumption inflexp, absorb(date)  // still negative
estout reg1 reg2, cells(b(star fmt(3)) se) stats(N r2)
