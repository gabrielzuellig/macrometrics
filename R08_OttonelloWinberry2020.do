
/* 
This script replicated Ottonello & Winberry (2020) with firm-level balance
sheet data and high-frequency monetary policy shocks for the US und the euro area.
All firm-level variables have been anonomyzed and contain random noise to mask
the identity of the firm.
*/

*** Program setup
capture log close
clear all
macro drop _all 
set more off 
set segmentsize 2g
program drop _all  
global path "/Users/gabrielzullig/Dropbox/Teaching/Basel2024/macrometrics/"
cd $workpath


** Load raw data

use $path/data/data_balancesheetsUSEA.dta, replace
xtset firmID tq


** Firm data cleaning
* exclude from estimation sample: Ireland
tab country
drop if country == 8
* exclude from estimation sample: real estate, utilities, some fintech sectors
drop if sector_broad == 7 | sector_broad == 9 | sector == 22 | sector == 52 
* generate leverage (= liabilities / assets)
gen leverage = bs_LiabTotal/bs_AssetsTotal*100
* interaction variable: leverage in deviation from firm-level mean (only use
* within-firm variation)
bys firmID (tq): egen leverage_mean = mean(leverage)
gen leverage_devfrommean = leverage - leverage_mean
sum leverage_devfrommean, det
ssc install winsor2 
winsor2 leverage_devfrommean, cuts(1 99) replace  // remove strong outliers
sum leverage_devfrommean, det
* variables in logs
gen lAssetsFixed = 100*log(bs_AssetsFixed)  // log of fixed assets (=capital, outcome variable)
gen lAssetsTotal = 100*log(bs_AssetsTotal)	// log of total assets (control variable)
gen lRevenue = 100*log(rs_Revenue)   // log of revenue


** some descriptive statistics (know your data!)
* number of firms and number of observations by firm 
ssc install distinct 
distinct firmID 
bys firmID (tq): gen obsid = _n
bys firmID (tq): gen nobs = _N
hist nobs if obsid == 1
* distribution of leverage 
hist leverage if inrange(leverage, 0, 150)
tabstat leverage, by(country) s(mean p25 p50 p75 N)
gen year = yofd(dofq(tq))
reg leverage i.year
drop year
hist leverage_devfrommean
* other variables
hist ra_ebit_to_revenue if inrange(ra_ebit_to_revenue, -100, 100)
hist gr_revenue if inrange(gr_revenue, -100, 100)
ssc install binscatter 
binscatter lRevenue lAssetsFixed
ssc install reghdfe 
reghdfe lRevenue lAssetsFixed, absorb(firmID tq)
binscatter leverage lAssetsFixed
reg leverage lAssetsFixed
reghdfe leverage lAssetsFixed, absorb(firmID)


** Merge time series data (shocks and GDP growth) onto panel
merge m:1 tq region using $path/data/data_JK20shocksUSEA_long.dta, keep(match master)
tab tq if _merge == 1  // the master-only data are post-2019, drop
drop if _merge == 1
drop _merge 
merge m:1 tq region using $path/data/data_GDPUSEA_long.dta, keep(match) nogen
xtset firmID tq


** Local projections 
global maxhorizon 16
gen h = _n-1 if _n <= $maxhorizon + 1
gen beta_int = .
gen se_int = .
gen pval_int = .

foreach h of numlist 0/$maxhorizon {
		
	* prepare left-hand side variable, which changes in each iteration	
	cap drop lhs   // drop left-hand side variable only if it already exists
	gen lhs = F`h'.lAssetsFixed - L1.lAssetsFixed
	
	* regression
	reghdfe lhs c.mpspr_JK#cL.leverage_devfrommean L1.leverage L1.lAssetsTotal L1.gr_revenue L1.ra_current_to_totalassets cL.leverage_devfrommean#cL.GDPgrowth if tq < yq(2020,1), absorb(firmID i.sector#i.tq i.country#i.tq) vce(cl firmID tq)
	
	* save coefficients and standard errors in dataset
	qui replace beta_int = _b[c.mpspr_JK#cL.leverage_devfrommean] if _n == `h'+1
	qui replace se_int = _se[c.mpspr_JK#cL.leverage_devfrommean] if _n == `h'+1
	matrix b = r(table)
	qui replace pval_int = b[4,1] if _n == `h'+1
	
	* save estimation output and number of firms for each regression within it
	estadd scalar nfirms = e(N_clust1)
	quietly est sto reg_h`h'
	
}


** Presentation of estimation results
* print in log file 
cap mkdir $path/R08_OttonelloWinberry2020
log using $path/R08_OttonelloWinberry2020/estimates.log, replace
list h beta_int se_int pval_int if h <= $maxhorizon
log close 

* figure with confidence intervals
drop if beta_int == .
gen up_int = beta_int + 1.96*se_int 
gen lo_int = beta_int - 1.96*se_int
	
graph twoway (line beta_int h, lc(black)) ///
	(line up_int h, lc(black) lp(dash)) /// 
	(line lo_int h , lc(black) lp(dash)), ///
	yline(0, lc(gs10) lp(solid)) ///
	xsc(r(0 $maxhorizon)) xlab(0(4)$maxhorizon) ///
	title("MP shock x leverage: Interaction term") ///
	legend(r(1) pos(6) order(1 "Mean estimate" 2 "95% CI"))
graph export $path/R08_OttonelloWinberry2020/irf.png, replace
	
* regression outputs
esttab reg_h0 reg_h2 reg_h4 reg_h6 reg_h12 reg_h16, se ar2 scalar(nfirms) ///
	rename(c.mpspr_JK#cL.leverage_devfrommean  interact) keep(interact) coeflabel(interact "MP shock x leverage") ///
	var(20) star(* .10 ** .05 *** .01) mtitle("h=0" "h=2" "h=4" "h=6" "h=12" "h=16")


