
/* 
This script replicates Autor/Dorn/Hansen (2013) on the "China shock" based on
localities in the US and industry exposure to the rise in Chinese imports.
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

use $path/data/data_commzones.dta, replace
xtset czone period


** Some descriptive statistics
distinct czone  // 722 commuting zones, 2 observations each
tabstat l_sh_empl_mfg, by(period) s(p25 mean p75)
tabstat d_sh_empl_mfg, by(period) s(p25 mean p75)
hist d_sh_empl_mfg


** Different IV regressions to replicate column (1) of Table 3 in ADH2013
// Y = d_sh_empl_mfg
// X = d_IPW
// Z = d_IPWO
cap ssc install ivreg2
cap ssc install ranktest


// OLS 
reg d_sh_empl_mfg d_IPW period [aw=timepwt48]

// 2SLS with Stata command
ivreg2 d_sh_empl_mfg (d_IPW=d_IPWO) period [aw=timepwt48]

// 2SLS with manual steps
// 1st stage
reg d_IPW d_IPWO period [aw=timepwt48]
global gamma = _b[d_IPWO]
predict d_IPW_hat, xb  // calculate fitted values 
// 2nd stage 
reg d_sh_empl_mfg d_IPW_hat period [aw=timepwt48]

// Wald estimator 
// reduced-form 
reg d_sh_empl_mfg d_IPWO period [aw=timepwt48]
global delta = _b[d_IPWO]
global beta = $delta / $gamma 
disp "$beta"


** Other outcome variables

// total employment: less responsive because non-manufacturing slightly increases
ivreg2 d_sh_empl (d_IPW=d_IPWO) period [aw=timepwt48]
ivreg2 d_sh_empl_nmfg (d_IPW=d_IPWO) period [aw=timepwt48]
// unemployment: increases slightly
ivreg2 d_sh_unempl (d_IPW=d_IPWO) period [aw=timepwt48]
// wages in manufacturing: decrease
ivreg2 d_avg_lnwkwage_mfg (d_IPW=d_IPWO) period [aw=timepwt48]

