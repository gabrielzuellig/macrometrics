function [LP] = estimateLPznonlin(data, z, Fs, h, p, c_case, exdata, alpha)

% Prep parameters/containers
t = size(data,1);
n = size(data,2);
gamma1 = NaN*ones(h,n);
gamma1bands = NaN*ones(h,n,2);
gamma2 = NaN*ones(h,n);
gamma2bands = NaN*ones(h,n,2);
if isempty(alpha)
    alpha = 90;
end
zscore = abs(mynorminv((1-alpha/100)/2));
Fs = [NaN; Fs(1:end-1)]; % lagged by one period

%% Y matrix
Y = data;

%% X matrix: Generate lags
X = [];
for pp = 1:p
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
% right-hand side interaction
nX1 = size(X,2); % number of columns before interaction
X = [(1-Fs).*X, Fs.*X];


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
        nona = logical(sum(isnan([lhs,X,z]), 2) == 0);
        lhs = lhs(nona,:);
        rhs = [(1-Fs(nona)).*z(nona), Fs(nona).*z(nona), X(nona,:)];
        
        % reduced-form estimation (OLS with Newey-West standard errors)
        reg = nwest(lhs,rhs,h); % regression step
        vcov = rhs'*rhs;
        gamma1(hh+1,nn) = reg.beta(1);
        gamma2(hh+1,nn) = reg.beta(2);
        gamma1bands(hh+1,nn,:) = reg.beta(1) + [-1 1]*zscore*reg.se(1);
        gamma2bands(hh+1,nn,:) = reg.beta(2) + [-1 1]*zscore*reg.se(2);

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
LP.gamma1 = gamma1;
LP.gamma1bands = gamma1bands;
LP.gamma2 = gamma2;
LP.gamma2bands = gamma2bands;
LP.h = h;
