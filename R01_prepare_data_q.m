
%% Preamble
% This script reads in quarterly time series data from an .xlsx, computes some
% transformations and saves all series in a database called data_q_ready.mat, 
% later to be used for the estimation of time series models.


clear all; clc;

addpath('./_tbx') 
addpath('./_tbx/supportfct') 
if ~exist('data/tsq', 'dir')
    mkdir('data/tsq')
end


%% Import
xlsdata = readtable('data/data_q_in.xlsx', 'ReadRowNames', 0);
data = xlsdata{:,:};
labels = xlsdata.Properties.VariableNames;


%% Timing and labelling
time = data(:,1)+(1/4)*(data(:,2)-1);
data = data(:,3:end);
labels = labels(1,3:end);

labelmat = readtable('data/data_q_in.xlsx','Sheet','labels','ReadVariableNames',0);
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
var = {'GDP','PotGP','CPI','CoreCPI','NCONS','GDPDEF','LF','CONS',...
    'NINV','INV','GDPPC','HOUR','EMP','THOURS','WAGE','Mortgages','FedAss','Credit'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, col) = 100*log(data(:, col));
    end
end

% % y/y growth rates (generate new series)
var = {'GDP','CPI','CoreCPI','MortgGDP'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, end+1) = NaN;
        data(5:end, end) = data(5:end, col) - data(1:end-4, col);
    labels{length(labels) + 1} = strcat(var{vv},'Gr4');
    labels_print{length(labels_print) + 1} = strcat(labels_print{find(strcmp(labels, var{vv}))},' growth');
    end
end
labels_print(find(strcmp(labels, 'CPIGr4'))) = {'Inflation'};
labels_print(find(strcmp(labels, 'CoreCPIGr4'))) = {'Inflation'};

% first difference
var = {'GDP','CPI','CoreCPI','FFR','FFRshadow','GDPDEF','GDPPC','CONS','INV','WAGE'};
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
labels_print(find(strcmp(labels, 'dFFRshadow'))) = {'Interest rate change'};
labels_print(find(strcmp(labels, 'dCPI'))) = {'Inflation'};
labels_print(find(strcmp(labels, 'dCoreCPI'))) = {'Inflation'};
labels_print(find(strcmp(labels, 'dGDPDEF'))) = {'Inflation'};

% hp filtering
var = {'GDP'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:,end+1) = NaN;
        data(~isnan(data(:,col)),end) = data(~isnan(data(:,col)),col) - hpfilter(data(~isnan(data(:,col)),col), 1600);
        labels{length(labels) + 1} = strcat(var{vv},'hpgap');
        labels_print{length(labels_print) + 1} = strcat(labels_print{find(strcmp(labels, var{vv}))});
    end
end

% interest rate in deviation from natural rate
var = {'FFR','FFRshadow'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    col2 = find(strcmp(labels, 'Rstar'));
    if ~isempty(col)
        data(:,end+1) = NaN;
        data(~isnan(data(:,col)),end) = data(~isnan(data(:,col)),col) - data(~isnan(data(:,col)),col2);
        labels{length(labels) + 1} = strcat(var{vv},'dev');
        labels_print{length(labels_print) + 1} = strcat(labels_print{find(strcmp(labels, var{vv}))});
    end
end
    

%% Plot descriptive time series
xt00 = find(time==2000);
xt = [flip(xt00:-20:1) (xt00+20):20:size(data,1)];

for vv = 1:size(data,2)
    figure()
    plot(1:size(data,1), data(:,vv),'k','LineWidth',2)
    axis('tight')
    grid on
    set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/4 + 2000)
    title(labels_print{vv})
    set(gca,'FontSize',12)
    set(gcf,'paperpositionmode','auto')
    print(gcf,'-dpng','-loose',strcat('data/tsq/',labels{vv},'.png'));
end

close all

%% Export, housekeeping
data_lib = data;
labels_lib = labels;
printlabels_lib = labels_print;
clear col col2 vv var data labels labels_print labelmat xt xt00
save(strcat('data/data_q_ready.mat'), 'data_lib', 'labels_lib', 'printlabels_lib','time')

