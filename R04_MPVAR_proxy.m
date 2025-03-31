
%% Preamble
% This script replicates a small VAR in which a monetary policy shock is 
% identified using Gertler & Karadi high-frequency shocks as an external
% instrument


clear all; clc;

addpath('./_tbx')   % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 
addpath('./_tbx/var_tbx') 


%% Settings
folder = 'R04_MPVAR_proxy/';   
if ~exist('folder', 'dir')
    mkdir(folder)
end
p = 12;
h = 48;
c_case = 1;
vars = {'GDPgap','Unemp','CoreCPIGr12','EBP','FFR'};
exvars = {};
n = length(vars);
ident = 'proxy';
proxyvar = {'mpsprGK4'};
shockpos = 5; % position of the variable to instrument
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


%% Identification: Instrumental variable (z)
% prepare
VAR.ident = ident;
VAR.shock = zeros(n,1);
VAR.shock(shockpos) = shocksize;
VAR.z = z(p+1:end,:);
VAR.shockpos = shockpos;
s = NaN*ones(1,n);
b = logical(1:n == shockpos);
% save residuals
u_p = VAR.u(:,b);  % reduced-form residuals to be instrumented
u_q = VAR.u(:,~b); % all other residuals
% exclude covid period (if applicable)
u_p = u_p(timeinreg(p+1:end) == 1);
u_q = u_q(timeinreg(p+1:end) == 1,:);
z = VAR.z(timeinreg(p+1:end) == 1);
% only use periods for which z is not NA
nona = logical(isnan(z) == 0);
u_p = u_p(nona,1);
u_q = u_q(nona,:);
z = z(nona,1);
% 1st-stage of IV: u_p on z
beta1 = [ones(size(z,1),1),z]\u_p; % elast.: +1.2
u_p_hat = [ones(size(z,1),1), z]*beta1;  % fitted values of residuals
xi = u_p - u_p_hat;
disp('First-stage F-statistic is: ')
Fstat = ((u_p'*u_p-xi'*xi)/(length(beta1)))/((xi'*xi)/(length(z)-length(beta1)))
figure()
plot([u_p_hat, z, u_p])
% 2nd-stage
s(~b) = u_p_hat\u_q;
s(b) = 1; % normalize s at shockpos to 1
% impulse responses
VAR.IRF = zeros(h,n);
for hh=1:h
    VAR.IRF(hh,:) = (VAR.C(:,:,hh)*s')';
end


%% Bootstrap
[VAR.IRFbands] = bootstrapVAR(VAR, nboot, alpha, 'wild');


%% Scaling and cleaning up
%so far, s gives us the vector of relative responses. We want absolute values for a 1sd
%shock. To get that, follow Piffer (2020):
if shocksize == 0   
    sigma11 = VAR.Omega(b,b);
    sigma12 = VAR.Omega(b,~b);
    sigma21 = VAR.Omega(~b,b);
    sigma22 = VAR.Omega(~b,~b);
    mu = s(~b)';
    Gamma = sigma22 + mu*sigma11*mu' - sigma21*mu' - mu*sigma21';
    b11b11prime = sigma11 - (sigma21-mu*sigma11)'*Gamma^(-1)*(sigma21-mu*sigma11);
    b11 = chol(b11b11prime)';
    s = s*b11;
    VAR.shock(shockpos) = b11;
    clear sigma11 sigma12 sigma21 sigma22 mu Gamma b11b11prime
else
    b11 = shocksize;
    s = b11 * s;
end
VAR.IRF = VAR.IRF*b11;
VAR.IRFbands = VAR.IRFbands*b11;
VAR.s = s;
VAR.Fstat = Fstat;
clear u_p u_q nona stage1 beta1 u_p_hat regsummary Fstat s b b11


%% Plot impulse responses
plotirf1(VAR.IRF, VAR.IRFbands, printvars, strcat(folder,'lin'))
save(strcat(folder, 'out.mat'), 'VAR')
