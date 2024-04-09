
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
global path = "/Users/gabrielzullig/Dropbox/Teaching/Basel2024/macrometrics/"


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


** Recode survey variables
* u_p4: unemployment outlook (numeric value with label)
recode unemp (1=1) (3=0) (5=-1) (8=.) (9=.), gen(u_p4)
tab u_p4, m
label define u_label 1 "higher" 0 "same" -1 "lower"
label values u_p4 u_label
tab u_p4, m
graph bar, over(u_p4) b1title("Unemployment outlook")
* ec_p4: economic outlook (numeric value without label, just to see)
recode bus12 (1=1) (2=0.5) (3=0) (4=-0.5) (5=-1) (8=.) (9=.), gen(ec_p4)
tab ec_p4, m
graph bar, over(ec_p4) b1title("Economic outlook")
* pi_p4: inflation expectations (numeric value)
sum px1, det
gen pi_p4 = px1 if abs(px1) < 95 // abs > 95 are value like don't know or invalid answers
sum pi_p4, det
hist pi_p4


** Keep only the variables relevant (for now)
keep caseid id date hasprevious prevdate idprev u_p4 ec_p4 pi_p4


** Average expected inflation over time 
preserve 

collapse (mean) pi_p4, by(date)
graph twoway (line pi_p4 date), xtitle("") ytitle("Average expected inflation")

restore 


** Regress inflation outlook on unemployment
* Download packages if necessary (only required once)
ssc install estout
ssc install ftools
ssc install reghdfe
* Plain
eststo reg1a: reg pi_p4 u_p4  // positive sign: high unemployment <=> higher inflation expectations
* Include time fixed effect 
eststo reg2a: reghdfe pi_p4 u_p4, absorb(date)
esttab reg1a reg2a, se ar2
* repeat everything with economic outlook (needs to be done now because data will be lost)
eststo reg1b: reg pi_p4 ec_p4 
eststo reg2b: reghdfe pi_p4 ec_p4, absorb(date)


** Match previous interview with same person and calculate revisions 
* First, tempsave id, date and survey variables but rename them as lags
preserve 
keep id date *_p4
rename (id date *_p4) (idprev prevdate *_p41)
save $path/data/tempdata, replace
restore 
* Second, merge tempsaved data onto original dataset 
merge m:1 idprev prevdate using $path/data/tempdata, keep(match master)
drop caseid id idprev _merge 
order date
* Transform wide to long 
keep if hasprevious == 1
gen hhid = _n
rename (date *_p4) (date2 *_p42)  // second interview
rename (prevdate) (date1)   // previous, i.e. first interview
reshape long date u_p4 ec_p4 pi_p4, i(hhid) j(intid)
* xtset data 
xtset hhid date


** Regress inflation outlook on unemployment with time and person fixed effect 
eststo reg3a: reghdfe pi_p4 u_p4, absorb(date hhid)
esttab reg1a reg2a reg3a, se ar2
* repeat everything with economic outlook 
eststo reg3b: reghdfe pi_p4 ec_p4, absorb(date hhid)
esttab reg1b reg2b reg3b, se ar2

