function [LP] = estimateLP(data, h, p, c_case, exdata, alpha, timeinreg)

% Prep parameters/containers
t = size(data,1);
n = size(data,2);
Gamma = NaN*ones(n,n,h);
Gammalo = NaN*ones(n,n,h);
Gammaup = NaN*ones(n,n,h);
if isempty(alpha)
    alpha = 90;
end
zscore = abs(mynorminv((1-alpha/100)/2));
timeinreg_all = timeinreg;


%% Y matrix
Y = data;

%% X matrix: Generate lags
X = [];
for pp = 0:p
    X = [X, [NaN*ones(pp, n); data(1:end-pp,:)] ];
end
% add constant and/or trend to X matrix
if c_case == 0
    c = [];
elseif c_case == 1
    c = ones(size(X,1),1); 
elseif c_case == 2
    c = [ones(size(X,1),1), [1:1:size(X)]' ];
else
    print('c_case variable needs to be set to 0, 1, or 2.')
end
% add exogenous variables
if ~isempty(exdata)
    Xex = [];
    for pp = 1:p
        Xex = [Xex, [NaN*ones(pp, size(exdata,2)); exdata(1:end-pp,:)]  ];
    end
else
    Xex = [];
end
X = [c, X, Xex];

%% Estimation
disp('Progress \n')
counter = 1;
nloops = n*h;
for nn = 1:1:n
    for hh = 0:1:(h-1)

        % prepare LHS variable (changes in each loop while RHS stays
        % constant)
        lhs = [Y(hh+1:end,nn); NaN*ones(hh,1)];

        % drop observations with any NA's
        nona = logical(sum(isnan([lhs,X]), 2) == 0);
        lhs = lhs(nona,:);
        rhs = [X(nona,:)];

        % drop covid period
        timeinreg = timeinreg_all(nona,:);
        lhs = lhs(logical(timeinreg),:);
        rhs = rhs(logical(timeinreg),:);
        
        % reduced-form estimation (OLS with Newey-West standard errors)
        reg = nwest(lhs,rhs,h); % regression step
        vcov = rhs'*rhs;
        Gamma(nn,:,hh+1) = reg.beta(c_case+1:c_case+n);
        Gammalo(nn,:,hh+1) = reg.beta(c_case+1:c_case+n) - zscore*reg.se(c_case+1:c_case+n);
        Gammaup(nn,:,hh+1) = reg.beta(c_case+1:c_case+n) + zscore*reg.se(c_case+1:c_case+n);
        
        % show progress report
        fprintf(1,'\b\b\b\b%3.0f%%',100*(counter/nloops)); 
        counter = counter + 1;

    end
end
fprintf('\n');

%% Housekeeping
LP.data = data;
LP.X = X;
LP.Xex = Xex;
LP.c_case = c_case;
LP.p = p;
LP.t = t;
LP.n = n;
LP.Gamma = Gamma;
LP.Gammalo = Gammalo;
LP.Gammaup = Gammaup;
LP.h = h;
