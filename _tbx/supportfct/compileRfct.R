
# Load packages
if (!require('gdata')) install.packages('gdata')
library(gdata)
if (!require('mFilter')) install.packages('mFilter')
library(mFilter)
if (!require('sandwich')) install.packages('sandwich')
library(sandwich)
if (!require('R.utils')) install.packages('R.utils')
library(R.utils)

# Load own functions
source('_tbx/lp_tbx/estimateLP.R')
source('_tbx/lp_tbx/estimateLPz.R')
source('_tbx/lp_tbx/estimateLPnonlin.R')
source('_tbx/lp_tbx/estimateLPznonlin.R')
source('_tbx/var_tbx/estimateVAR.R')
source('_tbx/var_tbx/dyn_multipliers.R')
source('_tbx/var_tbx/bootstrapVAR.R')
source('_tbx/var_tbx/testIRFsign.R')
#source('_tbx/ivar_tbx/simGIRF.R')
#source('_tbx/ivar_tbx/bootstrapIVAR.R')
source('_tbx/supportfct/plotirf1.R')
source('_tbx/supportfct/plotirf2.R')
source('_tbx/supportfct/plotirf3.R')

# Other settings
par(mar=c(4.6, 4.1, 2.1, 2.1))
