
%% Preamble
% This script replicates a small VAR in which a monetary policy shock is 
% identified using the Cholesky decomposition, similar to Christiano, 
% Eichenbaum & Evans (1999), but with notable differences (e.g. monthly frequency)


clear all; clc;

addpath('./_tbx')   % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 
addpath('./_tbx/var_tbx') 


%% Settings
folder = 'R02_MPVAR_Cholesky/';   
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
shocksize = 0; % 0 = one standard deviation, all else: absolute values
state.nonlinear = 'no';
nboot = 1000;
alpha = 90; %confidence level
ignoreyears = [2020, 2021];  % years to ignore in regression (e.g. covid)


%% Load data, generate lags, subset relevant subsample, etc.
load('data/data_m_ready.mat')
subset_data;

    
%% Estimate matrices A, Omega, S, dynamic multipliers
VAR = estimateVAR(data, p, c_case, exdata, timeinreg);
VAR.C = dyn_multipliers(VAR, h); % not identified


%% Identification: Cholesky
% organization
VAR.ident = ident;
VAR.shock = zeros(n,1);
VAR.shock(shockpos) = 1;
if shocksize ~= 0  % absolute values, e.g. 25bp = 0.25
    VAR.shock = VAR.shock ./ VAR.S(shockpos, shockpos) .* shocksize;
end
% Cholesky decomposition (already done in estimation function)
disp(VAR.S)
eps = VAR.u*inv(VAR.S);
VAR.eps = eps(:, shockpos);
clear eps
% impulse response functions identified
VAR.IRF = zeros(h,n);
for hh=1:h
    VAR.IRF(hh,:) = (VAR.C(:,:,hh)*VAR.S*VAR.shock)';
end


%% Bootstrapping
[VAR.IRFbands] = bootstrapVAR(VAR, nboot, alpha, 'residual');


%% Plot impulse responses and save VAR structure
plotirf1(VAR.IRF, VAR.IRFbands, printvars, strcat(folder,'lin'))
save(strcat(folder, 'out.mat'), 'VAR')

