
%% Preamble
% This script estimates local projections with Miranda-Agrippino shocks.
% As can be seen, the resulting impulse response functions are rather
% unstable, which is a frequent observation when estimated as a
% reduced-form model.


clear all; clc;

addpath('./_tbx')   % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 
addpath('./_tbx/var_tbx') 
addpath('./_tbx/lp_tbx') 


%% Settings
folder = 'R05c_MPLP_proxy/';   
if ~exist('folder', 'dir')
    mkdir(folder)
end
p = 6; % easily overfitted
h = 48;
c_case = 1;
vars = {'GDPgap','Unemp','CoreCPIGr12','FFR'};
exvars = {};
n = length(vars);
ident = 'proxy';
proxyvar = {'mpsprMA'};
shockpos = 4; % position of the shock in the Cholesky ordering
shocksize = 0; 
state.nonlinear = 'no';
alpha = 90; %confidence level
ignoreyears = [2020, 2021];  % years to ignore in regression (e.g. covid)


%% Load data, generate lags, subset relevant subsample, etc.
load('data/data_m_ready.mat')
subset_data;


%% Estimate dynamic multipliers
LP = estimateLPz(data, z, h, p, c_case, [], alpha, timeinreg);


%% Plot impulse responses
plotirf1(LP.gamma, LP.gammabands, printvars, strcat(folder,'lin'))
save(strcat(folder, 'out.mat'), 'LP')
