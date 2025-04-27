
/* 
This script estimates the effects of monetary policy shocks on firm-level investment. 
The dataset used contains balance sheet information on US and euro area firms. 
All firm-level variables have been anonomyzed and contain random noise to mask
the identity of the firm. In the last part, it replicates Ottonello & Winberry (2020)'s 
result that high-leverage firms are less responsive to monetary policy.
*/


// Setup
capture log close
clear all
macro drop _all 
set more off 
set segmentsize 2g
program drop _all  
global path "/Users/gabrielzullig/Dropbox/Teaching/Basel2025/macrometrics/"
cd $path


********************************************************************************
** Data preparation and exploration
********************************************************************************

// Load and format data
use $path/data/data_balancesheetsUSEA_in.dta, replace  // Load raw (anonymized) data
xtset firmID tq  // xtset formats the data as a panel (identifier: firmID, time: tq)

// Some data cleaning is necessary
tab country
drop if country == 8   // Exclude from estimation sample: Ireland
// Exclude from estimation sample: real estate, utilities, some fintech sectors:
drop if sector_broad == 7 | sector_broad == 9 | sector == 22 | sector == 52   

// Some descriptive statistics (part of "KYD!")
ssc install distinct 
tab tq    // Data from 1993 to 2023 (121 quarters)
distinct firmID    // Total of 3,811 firms, but not all firms in all periods
bys firmID (tq): gen obsid = _n
bys firmID (tq): gen nobs = _N
hist nobs if obsid == 1   // How many quarters do we data on for each firm?
// Descriptive statistics on size and growth
sum gr_revenue, det  // between -100 and 100%, with some crazy outliers
hist gr_revenue if inrange(gr_revenue, -100, 100)
ssc install binscatter   // binscatter is a good tool to quickly explore correlations between variables
gen lRevenue = 100*log(rs_Revenue)   // log of revenue
gen lAssetsFixed = 100*log(bs_AssetsFixed)  // log of fixed assets (=capital, outcome variable)
gen lAssetsTotal = 100*log(bs_AssetsTotal)	// log of total assets (control variable)
binscatter lRevenue lAssetsFixed    // positive relationship
// Compute leverage (= liabilities / assets) and show descriptive statistics
gen leverage = bs_LiabTotal/bs_AssetsTotal*100
sum leverage, det   // between 0 and around 140%, with some crazy outliers
hist leverage if inrange(leverage, 0, 150)   // mode is around 55%
tabstat leverage, by(country) s(mean p25 p50 p75 N)    // European firms have higher leverage than US
gen year = yofd(dofq(tq))
reg leverage i.year    // Leverage had an increasing trend from 1993-2013, since then roughly flat (with peak in 2020)
drop year


********************************************************************************
** Effects of monetary policy on investment 
********************************************************************************

// prepare variables 
gen lI = lAssetsFixed - L1.lAssetsFixed   // "L1." refers to the first lag of a variable; works only if the data is xtset.
gen epsM = mpspr_JK_US
global controls = "L1.leverage L1.lAssetsTotal L1.gr_revenue L1.ra_current_to_totalassets"  // variables we will use as control variables. They can be stored in a "global" X once, and then every time I refer to "$X", it calls all variables specified. Choice is motivated by OW20 paper, but could be improved.

// raw regression with no fixed effects 
reg lI epsM $controls  // b1 = -0.15 (p = 0.84)

// add firm fixed effects 
ssc install reghdfe   // reghdfe allows for efficient regressions with many fixed effects 
reghdfe lI epsM $controls, absorb(firmID)    // b1 = -0.89 (p = 0.22)

// add time fixed effects 
reghdfe lI epsM $controls, absorb(firmID tq)    // b1 cannot be estimated because of collinearity with time fixed effects 

// solve via interaction with US dummy, i.e., use European firms as a control group
gen USdummy = 0 
replace USdummy = 1 if country == 1
tab country USdummy
reghdfe lI epsM c.epsM#ib0.USdummy $controls, absorb(firmID i.sector#i.tq)     
// b2 (on epsM x USdummy) is now the coefficient of interest. It is unfortunately still insignificant.

// add dynamics and try with local projections 
global maxhorizon 16 
gen h = _n-1 if _n <= $maxhorizon + 1
gen beta2 = .   // prepare container for estimates
gen se2 = .
foreach h of numlist 0/$maxhorizon {   // loop over horizons h from 0 to 16; in each iteration, the index is stored in a local `h'
	cap drop lhs // drop the variable lhs if it exists 
	gen lhs = F`h'.lAssetsFixed - L1.lAssetsFixed  // compute left-hand side variable for this iteration 
	reghdfe lhs epsM c.epsM#ib0.USdummy $controls if tq < yq(2020,1), absorb(firmID i.sector#i.tq)   // actual regressions  
	qui replace beta2 = _b[1.USdummy#c.epsM] if _n == `h'+1   // save coefficients as a variable
	qui replace se2 = _se[1.USdummy#c.epsM] if _n == `h'+1
}
list h beta2 se2 if h <= $maxhorizon  // read out and plot estimates
gen up2 = beta2 + 1.96*se2 
gen lo2 = beta2 - 1.96*se2
graph twoway (line beta2 h, lc(black)) ///
	(line up2 h, lc(black) lp(dash)) /// 
	(line lo2 h , lc(black) lp(dash)), ///
	yline(0, lc(gs10) lp(solid)) ///
	xsc(r(0 $maxhorizon)) xlab(0(4)$maxhorizon) ///
	title("MP shock") ///
	legend(r(1) pos(6) order(1 "Mean estimate" 2 "95% CI"))

	
/* Using the tools so far, we cannot document a negative effect of monetary policy
	on investment in our sample of firms, at least not at the horizon where we 
	would normally expect it (ca. 1-2 years). Possible explanations:
	 - It does not exist (unlikely, given time series evidence)
	 - There is too much noise in the data (possible, but here unlikely, given results in next section)
	 - All these firms are multinational firms, so even the EU firms are exposed 
	   to US monetary policy shocks. They are a bad control group; thus a different 
	   specification is needed, either without time fixed effects or using a different type of variation
	 - .... 
*/


********************************************************************************
** OW20 replication: How leverage affects to monetary policy effects on investment
********************************************************************************

// prepare shocks 
drop epsM 
gen epsM = mpspr_JK_US if country == 1    // use US shocks for US firms...
replace epsM = mpspr_JK_EA if country != 1    // ... and EA shocks for EA firms

// prepare interaction variable: leverage in deviation from firm-level mean (only use
* within-firm variation)
bys firmID (tq): egen leverage_mean = mean(leverage)
gen leverage_devfrommean = leverage - leverage_mean
sum leverage_devfrommean, det    // crazy outliers: winsorize the variable at the 1st and 99th percentiles 
ssc install winsor2 
winsor2 leverage_devfrommean, cuts(1 99) replace
sum leverage_devfrommean, det
hist leverage_devfrommean

// add interaction of leverage with country's GDP growth to vector of control variables, as in OW20 paper
global controls = "$controls cL.leverage_devfrommean#cL.GDPgrowth"

// run local projections 
drop h beta2 se2
gen h = _n-1 if _n <= $maxhorizon + 1
gen beta2 = .   // prepare container for estimates
gen se2 = .
foreach h of numlist 0/$maxhorizon {
	cap drop lhs   // drop left-hand side variable only if it already exists
	gen lhs = F`h'.lAssetsFixed - L1.lAssetsFixed
	reghdfe lhs c.epsM#cL.leverage_devfrommean $controls if tq < yq(2020,1), absorb(firmID i.sector#i.tq i.country#i.tq) vce(cl firmID tq)
	qui replace beta2 = _b[c.epsM#cL.leverage_devfrommean] if _n == `h'+1
	qui replace se2 = _se[c.epsM#cL.leverage_devfrommean] if _n == `h'+1	
}
list h beta2 se2 if h <= $maxhorizon
drop up2 lo2
gen up2 = beta2 + 1.96*se2 
gen lo2 = beta2 - 1.96*se2
graph twoway (line beta2 h, lc(black)) ///
	(line up2 h, lc(black) lp(dash)) /// 
	(line lo2 h , lc(black) lp(dash)), ///
	yline(0, lc(gs10) lp(solid)) ///
	xsc(r(0 $maxhorizon)) xlab(0(4)$maxhorizon) ///
	title("MP shock x leverage: Interaction term") ///
	legend(r(1) pos(6) order(1 "Mean estimate" 2 "95% CI"))


