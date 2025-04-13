
%% Preamble
% This script estimates local projections with Romer&Romer shocks in
% expansions vs. recessions.


clear all; clc;

addpath('./_tbx')   % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 
addpath('./_tbx/var_tbx') 
addpath('./_tbx/lp_tbx') 


%% Settings
folder = 'R05d_MPLP_nonlin_extinst/';   
if ~exist('folder', 'dir')
    mkdir(folder)
end
p = 6;
h = 48;
c_case = 1;
vars = {'GDPgap','Unemp','CoreCPIGr12','FFR'};
exvars = {};
n = length(vars);
ident = 'proxy';
proxyvar = {'mpsprRR'};
shockpos = 4; % position of the shock in the Cholesky ordering
shocksize = 0; 
state.nonlinear = 'yes'; 
state.logistic = 'yes';
state.interacted = 'no';
state.statevar = {'Unemp'}; 
state.gamma = -5;  % because gamma enters negatively, a double negative indicates positive values = regime 2
state.cq = 60;
alpha = 90; %confidence level
ignoreyears = [2020, 2021];  % years to ignore in regression (e.g. covid)


%% Load data, generate lags, subset relevant subsample, etc.
load('data/data_m_ready.mat')
subset_data; 


%% Estimate dynamic multipliers
LP = estimateLPznonlin(data, z, state.Fs, h, p, c_case, [], alpha, timeinreg);


%% Plot impulse responses
plotirf2(LP.gamma1, LP.gamma1bands, LP.gamma2, LP.gamma2bands, printvars, strcat(folder,'nonlin'), {'Expansion','Recession'},'NorthEast')
save(strcat(folder, 'out.mat'), 'LP')
