
%% Preamble
% This script replicates Born, Müller, Schularick & Sedlacek (2019) and
% provides an update to the cost of Brexit with updated data.
% Unfortunately, the code (for the estimation of the coefficients) only works 
% when the Global Optimization Toolbox is activated. 


clear all; clc;

addpath('./_tbx/optimfct') % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 


%% Settings
folder = 'R08_Brexit_SCM/';   
if ~exist('folder', 'dir')
    mkdir(folder)
end
start_date = 1995;         % let estimation sample start in 1995q1
treatment_date = 2016.25;  % include up to 2016q2 in estimation sample
end_date = 2018.75;        % original end date of data by BMSS19
end_date_extd = 2023.5;    % extended end date


%% Load and examine pre-prepared data

load('./data/data_GDPOECDwide_in.mat');

% GDP_dep ("dependent" variable): real GDP of UK in levels, indexed to 1995 = 1
% This is the main variable we are trying to match for the pre-Brexit period
plot(time,GDP_dep)
% GDP_j: real GDP of 23 OECD economies in levels, indexed to 1995 = 1
plot(time,GDP_j)
ncandidates = size(GDP_j,2);

% calculate log levels 
lGDP_dep = log(GDP_dep);
lGDP_j = log(GDP_j);

% other characteristics for the UK and all other economies
T = table(othercharmeans_dep, 'VariableNames', {countryname_dep}, ...
    'RowNames',otherchar_names);
disp(T)   % these are the additional characteristics we are trying to match.


%% Build the data for the minimization routine
% combine GDP and other characteristics
Y = [lGDP_dep(time >= start_date & time <= treatment_date, 1); othercharmeans_dep];
X = [lGDP_j(time >= start_date & time <= treatment_date, :); othercharmeans_j];


%% Try to estimate weights with unrestricted OLS
w_OLS = X\Y;  
T = table(round(w_OLS,3), 'VariableNames', {'OLS weights'},'RowNames',countrynames_j);
disp(T)  % some (quite many) countries have negative weights
sum(w_OLS)   % and they don't sum to 1 (that the sum is close to 1 is a coincidence)


%% Regression with constraints that w_j >= 0 and sum_j w_j = 1
% Preparations: starting values       
s = std([Y X],[],2);
s1 = s(1);
s2 = s(2:end);
v20 =log((s1./s2).^2);
clear s s1 s2
H0 = 1e-2*eye(length(v20)); %Initial Hessian for csminwel function
crit = 1e-8; %Tolerance
nit = 10000;  %Number of iterations
        
% optimizer in two steps
[fminv1,v2] = csminwel(@loss_function_transformed,v20,H0,[],crit,nit,Y,X,Y,X);
clear v20 H0 crit nit fminv1
options = optimset('MaxFunEvals', 100000, 'Display', 'iter', 'MaxIter', 100000, 'TolFun', 1e-11, 'TolX', 1e-3);
[v22,fminv,exitflag] = fminsearch('loss_function_transformed',v20,options,Y,X,Y,X);
v = [1;exp(v22)]; % retransform
clear v2 v22 fminv exitflag
% Recover weights
D = diag(v);
H = X'*D*X;
H = (H+H')/2;
f = - Y'*D*X;
options = optimoptions('quadprog', 'Display', 'off');
[w_SCM,fval,e]=quadprog(H,f,[],[],ones(1,ncandidates),1,zeros(ncandidates,1),ones(ncandidates,1),[],options);
clear X Y D v d H f fval options e
T = table(round(w_SCM,3), 'VariableNames', {'weight'},'RowNames',countrynames_j);
disp(T)  % all weights are positive, some (quite many) have zero weights 
sum(w_SCM)   % they sum to exactly 1.


%% Verify estimation results

% other characteristics and how the doppelganger fits the UK
othercharmeans_dep_hat = othercharmeans_j * w_SCM;  
T = table(othercharmeans_dep, othercharmeans_dep_hat, ...
    'VariableNames',{countryname_dep,'Doppelganger'}, ...
    'RowNames',otherchar_names);
disp(T)

% Fit of doppelganger over estimation sample (pre-treatment)
time_estsample = time(time >= start_date & time <= treatment_date, 1);
lGDP_dep_estsample = 100*lGDP_dep(time >= start_date & time <= treatment_date, :);
lGDP_dep_hat_estsample = 100*lGDP_j(time >= start_date & time <= treatment_date, :) * w_SCM;
u_dep_estsample = lGDP_dep_estsample - lGDP_dep_hat_estsample;
sigma = std(u_dep_estsample);  % to construct confidence intervals
lo = lGDP_dep_hat_estsample - sigma;
up = lGDP_dep_hat_estsample + sigma;

figure()
shadedplot(time_estsample,lo',up',[.9 .9 .9],[.9 .9 .9] );
hold on
p1=plot(time_estsample, lGDP_dep_estsample,'b-','LineWidth', 1.5);
p2=plot(time_estsample, lGDP_dep_hat_estsample,'r-','LineWidth', 1.5);
xline(treatment_date,'k:','Brexit Vote')
ylabel('deviation from 1995 (percent)')
ylim([-5 60])
grid on
title('Pre-treatment fit')
legend([p1 p2], {countryname_dep,'Doppelganger'}, 'Location', 'NorthWest')
set(gcf, 'position', [0 0 600 300]);
print(gcf,'-dpng','-loose',strcat(folder,'doppelganger_pre.png'));


%% Post-treatent comparison as in paper (until end-2018)
plot_start_date = 2015;  % zoomed in, i.e. only start 2015
time_post = time(time >= plot_start_date & time <= end_date);
lGDP_dep_post = 100*lGDP_dep(time >= plot_start_date & time <= end_date, :);
lGDP_dep_hat_post = 100*lGDP_j(time >= plot_start_date & time <= end_date, :) * w_SCM;
% compute everything relative to 2016q2
lGDP_dep_hat_post = lGDP_dep_hat_post - lGDP_dep_post(find(time_post == treatment_date));
lGDP_dep_post = lGDP_dep_post - lGDP_dep_post(find(time_post == treatment_date));
lo = lGDP_dep_hat_post - sigma;
up = lGDP_dep_hat_post + sigma;

figure()
shadedplot(time_post,lo',up',[.9 .9 .9],[.9 .9 .9] );
hold on
p1=plot(time_post, lGDP_dep_post,'b-','LineWidth', 1.5);
p2=plot(time_post, lGDP_dep_hat_post,'r-','LineWidth', 1.5);
xline(treatment_date,'k:','Brexit Vote')
ylabel('deviation from 2016q2 (percent)')
ylim([-4 8])
grid on
title('Treatment effect')
legend([p1 p2], {countryname_dep,'Doppelganger'}, 'Location', 'NorthWest')
set(gcf, 'position', [0 0 600 300]);
print(gcf,'-dpng','-loose',strcat(folder,'doppelganger_post.png'));


%% Post-treatment comparison extended until 2023

% show extended GDP data: there have been some revisions in the dynamic in the meantime
figure()
plot(timeextd, GDPextd_dep)
hold on
plot(time, GDP_dep)

% show extended GDP series of candidate pool
figure()
plot(timeextd, GDPextd_j)

% compute post-treatment doppelganger based on existing weights, but extended data
time_post_extd = timeextd(timeextd >= plot_start_date & timeextd <= end_date_extd);
lGDP_dep_postextd = 100*log(GDPextd_dep(timeextd >= plot_start_date & timeextd <= end_date_extd, :));
lGDP_dep_hat_postextd = 100*log(GDPextd_j(timeextd >= plot_start_date & timeextd <= end_date_extd, :)) * w_SCM;
% compute everything relative to 2016q2
lGDP_dep_hat_postextd = lGDP_dep_hat_postextd - lGDP_dep_postextd(find(time_post_extd == treatment_date));
lGDP_dep_postextd = lGDP_dep_postextd - lGDP_dep_postextd(find(time_post_extd == treatment_date));
lo = lGDP_dep_hat_postextd - sigma;
up = lGDP_dep_hat_postextd + sigma;

figure()
shadedplot(time_post_extd,lo',up',[.9 .9 .9],[.9 .9 .9] );
hold on
p1=plot(time_post_extd, lGDP_dep_postextd,'b-','LineWidth', 1.5);
p2=plot(time_post_extd, lGDP_dep_hat_postextd,'r-','LineWidth', 1.5);
xline(treatment_date,'k:','Brexit Vote')
ylabel('deviation from 2016q2 (percent)')
ylim([-10 20])
grid on
title('Treatment effect with extended data')
legend([p1 p2], {countryname_dep,'Doppelganger'}, 'Location', 'SouthWest')
set(gcf, 'position', [0 0 600 300]);
print(gcf,'-dpng','-loose',strcat(folder,'doppelganger_postextd.png'));


