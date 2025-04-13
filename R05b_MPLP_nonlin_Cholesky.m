
%% Preamble
% This script estimates a nonlinear model of monetary policy whose effects
% depend on whether or not the interest rate in the previous 12 months has
% increased or decreased, inspired by Berger et al (2021)


clear all; clc;

addpath('./_tbx')   % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 
addpath('./_tbx/var_tbx') 
addpath('./_tbx/lp_tbx') 


%% Settings
folder = 'R05b_MPLP_nonlin_Cholesky/';   
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
shocksize = 1; 
state.nonlinear = 'yes'; 
state.logistic = 'yes';
state.interacted = 'no';
state.statevar = {'dFFR'}; 
state.gamma = -200;  % because gamma enters negatively, a double negative indicates positive values = regime 2
state.cq = 50;
alpha = 90; %confidence level
ignoreyears = [2020, 2021];  % years to ignore in regression (e.g. covid)


%% Load data, generate lags, subset relevant subsample, etc.
load('data/data_m_ready.mat')
subset_data; 


%% Overwrite indicator function with interest level change of past 12 months
scum = state.s;
for pp = 1:11
    scum = scum + [repmat(NaN, pp, 1); state.s(1:end-pp)];
end
state.s = scum(12:end); % shorten by first 11 periods (NA's)
state.Fs = state.s > 0;
data = data(12:end,:); % this shortens the data by 11 periods
timeinreg = timeinreg(12:end,:);


%% Estimate dynamic multipliers
LP = estimateLPnonlin(data, state.Fs, h, p, c_case, [], alpha, timeinreg);


%% Identification: Estimate VAR to get S
VAR = estimateVAR(data, p, c_case, exdata, timeinreg);  % could/should do an SVAR here
LP.ident = ident;
LP.shock = zeros(n,1);
LP.shock(shockpos) = 1;
if shocksize ~= 0  % absolute values, e.g. 25bp = 0.25
    LP.shock1 = LP.shock ./ VAR.S(shockpos, shockpos) .* shocksize;
    LP.shock2 = LP.shock ./ VAR.S(shockpos, shockpos) .* shocksize;
end
% impulse responses
LP.IRF1 = zeros(h,n);
LP.IRF1bands = zeros(h,n,1);
LP.IRF2 = zeros(h,n);
LP.IRF2bands = zeros(h,n,1);
for hh=1:h
    LP.IRF1(hh,:) = (LP.Gamma1(:,:,hh)*VAR.S*LP.shock1)';
    LP.IRF1bands(hh,:,1) = (LP.Gamma1lo(:,:,hh)*VAR.S*LP.shock1)';
    LP.IRF1bands(hh,:,2) = (LP.Gamma1up(:,:,hh)*VAR.S*LP.shock1)';
    LP.IRF2(hh,:) = (LP.Gamma2(:,:,hh)*VAR.S*LP.shock2)';
    LP.IRF2bands(hh,:,1) = (LP.Gamma2lo(:,:,hh)*VAR.S*LP.shock2)';
    LP.IRF2bands(hh,:,2) = (LP.Gamma2up(:,:,hh)*VAR.S*LP.shock2)';
end


%% Plot impulse responses
plotirf2(LP.IRF1, LP.IRF1bands, LP.IRF2, LP.IRF2bands, printvars, strcat(folder,'nonlin'), {'Easing cycle','Tightening cycle'},'NorthEast')
save(strcat(folder, 'out.mat'), 'LP')
