
%% Preamble
% This script reads in monthly time series data from an .xlsx, computes some
% transformations and saves all series in a database called data_m_ready.mat, 
% later to be used for the estimation of time series models.


clear all; clc;

addpath('./_tbx') 
addpath('./_tbx/supportfct') 
if ~exist('data/tsm', 'dir')
    mkdir('data/tsm')
end


%% Import
xlsdata = readtable('data/data_m_in.xlsx', 'ReadRowNames', 0);
data = xlsdata{:,:};
labels = xlsdata.Properties.VariableNames;


%% Timing and labelling
time = data(:,1)+(1/12)*(data(:,2)-1);
data = data(:,3:end);
labels = labels(1,3:end);

labelmat = readtable('data/data_m_in.xlsx','Sheet','labels','ReadVariableNames',0);
labels_print = labels;
for vv = 1:length(labels_print)
    col = find(strcmp(labelmat{:,1}, labels(vv)));
    if ~isempty(col)
        labels_print(vv) = labelmat{col,2};
    end
end
clear xlsdata labelmat


%% Data treatment

% 100*log of variables in levels (replace series)
var = {'IP','CPI','CoreCPI','GDPIHS','PotGDP','PPI'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, col) = 100*log(data(:, col));
    end
end

% % y/y growth rates (generate new series)
var = {'IP','CPI','CoreCPI','PPI'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, end+1) = NaN;
        data(13:end, end) = data(13:end, col) - data(1:end-12, col);
    labels{length(labels) + 1} = strcat(var{vv},'Gr12');
    labels_print{length(labels_print) + 1} = strcat(labels_print{find(strcmp(labels, var{vv}))},' growth');
    end
end
labels_print(find(strcmp(labels, 'CPIGr12'))) = {'Inflation'};
labels_print(find(strcmp(labels, 'CoreCPIGr12'))) = {'Inflation'};
labels_print(find(strcmp(labels, 'PPIGr12'))) = {'Commodity price inflation'};

% first differences (generate new series)
var = {'FFR','FFRshadow','CPI','CoreCPI','GDPgap'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, end+1) = NaN;
        data(2:end, end) = data(2:end, col) - data(1:end-1, col);
    labels{length(labels) + 1} = strcat('d',var{vv});
    labels_print{length(labels_print) + 1} = strcat(labels_print{find(strcmp(labels, var{vv}))},' growth');
    end
end
labels_print(find(strcmp(labels, 'dFFR'))) = {'Interest rate change'};
labels_print(find(strcmp(labels, 'dCPI'))) = {'Inflation'};
labels_print(find(strcmp(labels, 'dCoreCPI'))) = {'Inflation'};

% moving average
var = {'dGDPgap'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, end+1) = movmean(data(:,col), [19 0]); % 20-period backward moving average
    labels{length(labels) + 1} = strcat(var{vv},'MA20');
    labels_print{length(labels_print) + 1} = strcat(labels_print{find(strcmp(labels, var{vv}))},' Trend');
    end
end

% change the volatility of the monthly Chicago index implies the same std of 
% its quarterly aggregate and the real GDP growth rate
col = find(strcmp(labels, 'oriChicago'));
chicagom = data(:,col);
year = floor(time);
month = (time - year)*12+1;
quarter = NaN*ones(size(month));
quarter(month <= 3.1) = 1;
quarter(month >= 3.9 & month <= 6.1) = 2;
quarter(month >= 6.9 & month <= 9.1) = 3;
quarter(month >= 9.9) = 4;
quarter = year + (quarter/4-0.25);
quarter = quarter(year < 1992);  % only need period prior to 1992
qvalues = unique(quarter);
chicagoq = NaN*ones(length(qvalues), 1);
for qq = 1:length(qvalues)
    tmp = chicagom(quarter == qvalues(qq));
    chicagoq(qq,1) = sum(tmp(~isnan(tmp)));
end

warning('off')
xlsdata = readtable('data/data_m_in.xlsx', 'sheet', 'help_quarterly_to_monthly','NumHeaderLines',0);
warning('on')
dataq = xlsdata{:,:}; clear xlsdata; 
GDPgrowth = dataq(1:size(chicagoq,1),4); 
plot(chicagoq)
hold on
plot(GDPgrowth)
chicagoq = chicagoq(2:end, 1);
%chicagoq = [ones(size(chicagoq)), chicagoq];
GDPgrowth = GDPgrowth(2:end, 1);
chicagoq'*chicagoq \ chicagoq'*GDPgrowth
clear chicagom chicagoq factor GDPgrowth month year quarter qq qvalues

% monthly output gap: compare own measure to quarterly CBO measure and
% hp-filtered version
dataq = dataq(1:215,[3 6]); % real GDP level; CBO output gap
dataq(:,1) = 100*log(dataq(:,1));  % 100*log(GDP)
dataq(:,3) = hpfilter(dataq(:,1), 1600); % hp-filtered potential
dataq(:,4) = dataq(:,1) - dataq(:,3); % hp-filtered gap
dataq = [ [1:1:size(dataq,1)]', dataq]; % copy quarterly series three times and order to get monthly
dataq3 = [dataq; dataq; dataq];
[~, idx] = sort(dataq3(:,1));
dataq3 = dataq3(idx, 2:end);
if size(dataq3,1) < size(data,1)
    dataq3 = [dataq3; NaN*ones(size(data,1)-size(dataq3,1),size(dataq3,2)) ];
end
col = find(strcmp(labels, 'GDPgap'));


%% Plot descriptive time series
xt00 = find(time==2000);
xt = [flip(xt00:-60:1) (xt00+60):60:size(data,1)];

% output gap measures
figure()
plot(1:size(data,1), data(:,col),'k','LineWidth',1.5)
hold on
grid on
plot(1:size(dataq3,1), dataq3(:,2),'-b','LineWidth',2)
plot(1:size(dataq3,1), dataq3(:,4),':r','LineWidth',2)
plot(1:size(data,1), zeros(size(data,1),1), 'k')
axis('tight')
set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
title('Output gap')
set(gca,'FontSize',14)
legend('Monthly own measure','Quarterly CBO estimated','Quarterly hp-filtered series (lambda = 1600)','Location','SouthWest')
legend boxoff
set(gcf,'paperpositionmode','auto')
set(gcf, 'position', [0 0 800 400]);
print(gcf,'-dpng','-loose',strcat('data/tsm/outputgaps.png'));
clear dataq dataq3 idx

% 'state' graph
s = data(:,find(strcmp(labels, 'GDPgap')));
s2 = (s - mean(s(~isnan(s))))/std(s(~isnan(s)));
s3 = (s - prctile(s,5))/std(s(~isnan(s)));

figure()
Fs = exp(-1.5*s2)./(1+exp(-1.5*s2));
plot(1:size(data,1), Fs, 'k', 'LineWidth', 2);
hold on
grid on
Fs = exp(-5*s2)./(1+exp(-5*s2));
plot(1:size(data,1), Fs, '--b', 'LineWidth', 2);
Fs = exp(-1.5*s3)./(1+exp(-1.5*s3));
plot(1:size(data,1), Fs, ':r', 'LineWidth', 2);
legend('\gamma = 1.5','\gamma = 5','\gamma = 1.5, \mu = q_{5}',...
    'Location','SouthOutside','Orientation','horizontal','FontSize',14)
legend boxoff
axis('tight')
set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
set(gca,'FontSize',12)
set(gcf,'paperpositionmode','auto')
set(gcf, 'position', [0 0 800 300]);
print(gcf,'-dpng','-loose',strcat('data/tsm/state.png'));
clear s s2 s3 Fs

for vv = 1:size(data,2)
    figure()
    plot(1:size(data,1), data(:,vv),'k','LineWidth',2)
    axis('tight')
    grid on
    set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
    title(labels_print{vv})
    set(gca,'FontSize',12)
    set(gcf,'paperpositionmode','auto')
    print(gcf,'-dpng','-loose',strcat('data/tsm/',labels{vv},'.png'));
end

close all


%% Export, housekeeping
data_lib = data;
labels_lib = labels;
printlabels_lib = labels_print;
clear col vv var data labels labels_print labelmat xt xt00
save(strcat('data/data_m_ready.mat'), 'data_lib', 'labels_lib', 'printlabels_lib','time')

