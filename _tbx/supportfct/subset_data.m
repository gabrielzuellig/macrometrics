
% select T x N matrix of endogenous variables
data = NaN*ones(size(data_lib,1), n);
printvars = {};
for vv = 1:length(vars)
    col = find(strcmp(labels_lib, vars(vv)));
    data(:,vv) = data_lib(:, col);
    printvars = [printvars, printlabels_lib(col)];
end

% select T x N matrix of exogenous variables
exdata = NaN*ones(size(data_lib,1), length(exvars));
if ~isempty(exvars)
    for vv = 1:length(exvars)
        col = find(strcmp(labels_lib, exvars(vv)));
        exdata(:,vv) = data_lib(:, col);
    end
end

% select instrument (if necessary)
z = [];
if strcmp(ident, 'proxy')
    for vv = 1:length(proxyvar)
        col = find(strcmp(labels_lib, proxyvar(vv)));
        z = [z, data_lib(:, col)];
    end
end

% select state (if necessary)
s = [];
if strcmp(state.nonlinear, 'yes') && strcmp(state.logistic, 'yes')
    col = find(strcmp(labels_lib, state.statevar));
    s = data_lib(:, col);
end


% select interaction (if necessary)
if strcmp(state.nonlinear, 'yes') && strcmp(state.interacted, 'yes')
    exdata = NaN*ones(size(data_lib,1), 1);
    if max(strcmp(vars, state.shockvar)) == 0 || max(strcmp(vars, state.statevar)) == 0
        error('Both shockvar and statevar need to be part of the vector of endogenous variables.')
    end
    col1 = find(strcmp(labels_lib, state.shockvar));
    col2 = find(strcmp(labels_lib, state.statevar));
    exdata = data_lib(:, col1) .* data_lib(:, col2);
    state.shockpos = find(strcmp(state.shockvar, vars));
    state.statepos = find(strcmp(state.statevar, vars));
end


% ignore non-NA columns
nona = logical(sum(isnan([data, exdata, s]),2) == 0);
data = data(nona, :);
exdata = exdata(nona, :);
time = time(nona, :);
if strcmp(ident, 'proxy')
    z = z(nona, :);
end
if strcmp(state.nonlinear, 'yes') && strcmp(state.logistic, 'yes')
    s = s(nona, :);
end


% ignore data before 1965
data = data(time >= 1965, :);
exdata = exdata(time >= 1965, :);
if strcmp(ident, 'proxy')
    z = z(time >= 1965, :);
end
if strcmp(state.nonlinear, 'yes') && strcmp(state.logistic, 'yes')
    s = s(time >= 1965, :);
end
time = time(time >= 1965, :);


% end period: ignore data from 2020 onwards (hard-coded for now, more
% flexible in next version of toolbox)
data = data(time < 2020, :);
time = time(time < 2020);
if strcmp(ident, 'proxy')
    z = z(time < 2020, :);
end
if strcmp(state.nonlinear, 'yes') && strcmp(state.logistic, 'yes')
    s = s(time < 2020, :);
end
exdata = exdata(time < 2020, :);


% Logistic transformation of s-variable
if strcmp(state.nonlinear, 'yes') && strcmp(state.logistic, 'yes')
    state.s_orig = s;
    state.s_orig_cons = prctile(s, state.cq);
    state.s_orig_std = std(s(~isnan(s)));
    state.s = (s - state.s_orig_cons) / state.s_orig_std;
    state.Fs = exp(-state.gamma*state.s)./(1+exp(-state.gamma*state.s));   
end

% dummy of state
if strcmp(state.nonlinear, 'yes') && strcmp(state.interacted, 'yes')
    state.absval = prctile(data(:,state.statepos), state.cq);
    state.s = logical(data(:,state.statepos) <= state.absval);
end


% generate plot of actual data that goes into model
frequency = sum(floor(time) == 2000); % 4 if quarterly, 12 if monthly
xt00 = find(time==2000);
if frequency == 4
    xt = [flip(xt00:-40:1) (xt00+40):40:size(data,1)];
elseif frequency == 12
    xt = [flip(xt00:-120:1) (xt00+120):120:size(data,1)];
end
nrow = ceil(size(data,2)/3);

figure()
for vv = 1:size(data,2)
    subplot(nrow, 3, vv)
    plot(1:size(data,1), data(:,vv),'k','LineWidth',2)
    axis('tight')
    grid on
    hold on
    plot(1:size(data,1), zeros(size(data,1),1), 'k')
    if frequency == 4
        set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/4 + 2000)
    elseif frequency == 12
        set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
    end
    set(gca,'FontSize',10)
    title(printvars(vv),'FontSize',14)
end
set(gcf,'paperpositionmode','auto')
set(gcf, 'position', [0 0 800 nrow*200]);
print(gcf,'-dpng','-loose',strcat(folder,'series.png'));

if strcmp(state.nonlinear, 'yes')
    figure()
    if strcmp(state.logistic, 'yes')
        subplot(3, 1, 1)
        plot(1:size(data,1), state.s_orig, 'k', 'LineWidth', 2)
        if min(state.s(~isnan(state.s))) < 0 && max(state.s(~isnan(state.s))) > 0
            hold on
            plot(1:size(data,1), zeros(size(data,1),1), 'k', 'LineWidth', .5)
        end
        axis('tight')
        grid on
        if frequency == 4
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/4 + 2000)
        elseif frequency == 12
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
        end
        ylabel('s, not standardized')
        set(gca,'FontSize',10)
        title('State variable','FontSize',14)
        subplot(3, 1, 2)
        temps = prctile(state.s_orig, 5):((prctile(state.s_orig, 95)-prctile(state.s_orig, 5))/99):prctile(state.s_orig, 95);
        tempshat = (temps - state.s_orig_cons) / state.s_orig_std;
        tempFs = exp(-state.gamma*tempshat)./(1+exp(-state.gamma*tempshat));
        scatter(state.s_orig(state.s_orig > min(temps) & state.s_orig < max(temps)), state.Fs(state.s_orig > min(temps) & state.s_orig < max(temps)), 'x')
        hold on
        plot(temps, tempFs, 'k', 'LineWidth', 2)
        axis('tight')
        grid on
        ylim([0 1 ])
        xlabel('s, not standardized')
        ylabel('F(s)')
        set(gca,'FontSize',10)
        title('Transition function','FontSize',14)
        clear temps tempshat tempFs
        subplot(3, 1, 3)
        plot(1:size(data,1), state.Fs, 'k', 'LineWidth', 2)
        axis('tight')
        grid on
        ylim([0 1 ])
        if frequency == 4
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/4 + 2000)
        elseif frequency == 12
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
        end
        ylabel('F(s)')
        set(gca,'FontSize',10)
        title('Regime indicator / logistic transformation','FontSize',14)
        set(gcf,'paperpositionmode','auto')    
        set(gcf, 'position', [0 0 800 600]);
    elseif strcmp(state.interacted, 'yes')
        subplot(2, 1, 1)
        plot(1:size(data,1), data(:,state.shockpos), 'k', 'LineWidth', 2)
        if min(data(~isnan(data(:,state.shockpos)),state.shockpos)) < 0 && max(data(~isnan(data(:,state.shockpos)),state.shockpos)) > 0
            hold on
            plot(1:size(data,1), zeros(size(data,1),1), 'k', 'LineWidth', .5)
        end
        axis('tight')
        grid on
        if frequency == 4
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/4 + 2000)
        elseif frequency == 12
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
        end
        set(gca,'FontSize',10)
        title('Interacted variable 1: Shock','FontSize',14)
        subplot(2, 1, 2)
        area([state.s*min(data(:,state.statepos)), ...
            state.s*(max(data(:,state.statepos))-min(data(:,state.statepos)))], ...
            'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [0.9 0.9 0.9])
        hold on
        plot(1:size(data,1), data(:,state.statepos), 'k', 'LineWidth', 2)
        if min(data(~isnan(data(:,state.shockpos)),state.shockpos)) < 0 && max(data(~isnan(data(:,state.shockpos)),state.shockpos)) > 0
            plot(1:size(data,1), zeros(size(data,1),1), 'k', 'LineWidth', .5)
        end
        axis('tight')
        legend('Regime 2', 'Location', 'NorthWest')
        legend boxoff
        grid on
        if frequency == 4
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/4 + 2000)
        elseif frequency == 12
            set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
        end
        set(gca,'FontSize',10)
        title('Interacted variable 2: State','FontSize',14)
        set(gcf,'paperpositionmode','auto')    
        set(gcf, 'position', [0 0 800 400]);
    end
    print(gcf,'-dpng','-loose',strcat(folder,'state.png'));
end

clear vv col nona xt00 xt nrow s col1 col2
