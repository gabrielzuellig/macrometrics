
/* 
This script replicates Hazell/Herreno/Nakamura/Steinsson (2022) on the slope 
of the Phillips curve, but for the euro area (at the country level) instead of the US
*/

*** Program setup
capture log close
clear all
macro drop _all 
set more off 
set segmentsize 2g
program drop _all  
global path "/Users/gabrielzullig/Dropbox/Teaching/Basel2025/macrometrics/"
cd $path


** Load raw data

use $path/data/data_EAcountries.dta, replace


** Generate necessary variables and subset
gen service_infl = (shicp_sa / L1.shicp_sa - 1 )*100
gen rel_service_price = shicp_sa / chicp_sa * 100
keep if time_period >= yq(2000,1)

** OLS 
reghdfe service_infl unemp rel_service_price, absorb(country time_period) cluster(country time_period)   // -0.050**

** IV
ivreg2 service_infl (unemp = L4.unemp) rel_service_price ib1.country ib164.time_period, cluster(country time_period)   // -0.054***

