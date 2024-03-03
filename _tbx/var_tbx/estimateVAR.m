function [VAR] = estimateVAR(data, p, c_case, exdata)

% Prep
t = size(data, 1);
n = size(data, 2);

%% Y matrix
Y = data;

%% X matrix: Generate lags
X = [];
for pp = 1:p
    X = [X, [NaN*ones(pp, n); data(1:end-pp,:)] ];
end
% lags introduce some NA's => truncate
X = X(p+1:end,:); 
Y = Y(p+1:end,:);
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
    Xex = Xex(p+1:end,:);
else
    Xex = [];
end
X = [c, X, Xex];

%% Estimation
A = X\Y;  % a_ij is reaction of column j to row i
u = Y - X*A;  % reduced-form residuals
Omega = cov(u);  % variance-covariance matrix
S = chol(Omega)';  % lower-triangular matrix (Cholesky decomposition)


%% Housekeeping
VAR.data = data;
VAR.exdata = exdata;
VAR.X = X(:,(c_case+1):(c_case+n*p));
VAR.Xex = Xex;
VAR.c_case = c_case;
VAR.p = p;
VAR.t = t;
VAR.n = n;
VAR.A = A;
VAR.u = u;
VAR.Omega = Omega;
VAR.S = S;

