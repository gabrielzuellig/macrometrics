
%% Preamble
% This script replicates our small VAR where the IRFs are estimated with
% local projections à la Jordà (2005) instead.


clear all; clc;

addpath('./_tbx')   % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 
addpath('./_tbx/var_tbx') 
addpath('./_tbx/lp_tbx') 


%% Settings
folder = 'R05a_MPLP_Cholesky/';   
if ~exist('folder', 'dir')
    mkdir(folder)
end
p = 12;
h = 48;
c_case = 1;
vars = {'GDPgap','Unemp','CoreCPIGr12','FFR'};
exvars = {};
n = length(vars);
ident = 'chol';
shockpos = 4; % position of the shock in the Cholesky ordering
shocksize = 0; 
state.nonlinear = 'no';
alpha = 90; %confidence level


%% Load data, generate lags, subset relevant subsample, etc.
load('data/data_m_ready.mat')
subset_data;


%% Estimate dynamic multipliers
LP = estimateLP(data, h, p, c_case, [], alpha);


%% Identification: Estimate VAR to get S
VAR = estimateVAR(data, p, c_case, []);
LP.ident = ident;
LP.shock = zeros(n,1);
LP.shock(shockpos) = 1;
if shocksize ~= 0  % absolute values, e.g. 25bp = 0.25
    LP.shock = LP.shock ./ VAR.S(shockpos, shockpos) .* shocksize;
end
% structural impulse responses
LP.IRF = zeros(h,n);
LP.IRFbands = zeros(h,n,1);
for hh=1:h
    LP.IRF(hh,:) = (LP.Gamma(:,:,hh)*VAR.S*LP.shock)';
    LP.IRFbands(hh,:,1) = (LP.Gammalo(:,:,hh)*VAR.S*LP.shock)';
    LP.IRFbands(hh,:,2) = (LP.Gammaup(:,:,hh)*VAR.S*LP.shock)';
end


%% Plot impulse responses
plotirf1(LP.IRF, LP.IRFbands, printvars, strcat(folder,'lin'))
save(strcat(folder, 'out.mat'), 'LP')
