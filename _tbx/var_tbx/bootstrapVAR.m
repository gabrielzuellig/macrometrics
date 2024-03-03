function [IRFbands, Boot] = bootstrapVAR(VAR, nboot, alpha, method)

% inputs: 1. VAR structure incl. data and estimated A, 2. number of
% bootstraps, 3. confidence level, 4. bootstrap method ('residual' or
% 'wild')

% unpacking
t = size(VAR.u,1);
n = VAR.n;
u = VAR.u;
A = VAR.A;
c_case = VAR.c_case;
X = VAR.X;
Xex = VAR.Xex;
p = VAR.p;
IRFb = NaN*ones(size(VAR.IRF,1), size(VAR.IRF,2), nboot);
h = size(VAR.IRF,1);
Boot.Ab = NaN*ones(size(A, 1), size(A, 2), nboot);
Boot.Cb = NaN*ones(size(VAR.C, 1), size(VAR.C, 2), size(VAR.C, 3), nboot);
Boot.ub = NaN*ones(size(u, 1), size(u, 2), nboot);
Boot.S = NaN*ones(n, n, nboot);

% define upper and lower thresholds
if isempty(alpha)
    alpha = 90;
end
lo = (100-alpha)/2;
up = 100 - lo;

% combine constant, endogenous and exogenous variables
if c_case == 0
    c = [];
elseif c_case == 1
    c = ones(size(X,1),1); 
elseif c_case == 2
    c = [ones(size(X,1),1), [1:1:size(X)]' ];
end
endogvars = c_case + [1:1:size(X,2)];
X = [c, X, Xex];
if ~isempty(Xex) && strcmp(method, 'residual')
    method = 'wild';
    disp('Changing method to Wild bootstrapping (due to exogenous variables).')
end
if strcmp(VAR.ident, 'proxy') && strcmp(method, 'residual')
    method = 'wild';
    disp('Changing method to Wild bootstrapping (due to external instrument).')
end


disp('Progress \n')
for bb = 1:1:nboot  % 'nboot' bootstraps

    % generate pseudo-disturbance 'ub'
    if strcmp(method, 'residual')
        % residual bootstrapping assumes E(y|Xb) = X*b and u are iid (i.e.
        % homoskedasticity => Draw from the empirical distribution of u
        
        segment = (1:t)/t;
        ub = zeros(ceil(1.25*t)+p,n);
        for i=1:size(ub,1)
            draw = rand(1,1);  % uniform distribution
            ub(i,:) = u(min(find(segment>=draw)),:);
        end
        
    elseif strcmp(method, 'wild')
        % allows for heteroskedasticity

        fu = 1-2*(rand(t,1)>0.5);
        ub = zeros(size(VAR.u));
        for i=1:n
            ub(:,i) = u(:,i).*fu;  % flip sign of randomly selected 50%
        end

    end

    % generate pseudo-sample based on drawn u's
    Yb = zeros(size(ub,1),n);
    r = X(1,:);
    for i=1:size(ub,1)
        Yb(i,:) = r*A + ub(i,:);
        r = [Yb(i,:), r(endogvars(1:end-n))];
        if c_case == 1
            r = [1, r];
        elseif c_case == 2
            r = [1, i, r];
        end
        if ~isempty(Xex)
            r = [r, Xex(i,:)];
        end
    end
    data = Yb(end-size(VAR.u)+1:end,:);
    if strcmp(VAR.ident, 'proxy')
        zb = VAR.z.*fu;
        zb = zb(p+1:end,:);
    end
        
    % estimate VAR and IRF
    VARb = estimateVAR(data, VAR.p, c_case, VAR.exdata(p+1:end,:));
    VARb.C = dyn_multipliers(VARb, h);
    Boot.Ab(:,:,bb) = VARb.A;
    Boot.Cb(:,:,:,bb) = VAR.C;
    Boot.ub(:,:,bb) = [NaN*ones(p, n); VARb.u ];
    Boot.S(:,:,bb) = VARb.S;

    % identify shock
    if strcmp(VAR.ident, 'chol')
        for hh=1:h
            IRFb(hh,:,bb) = (VARb.C(:,:,hh)*VARb.S*VAR.shock)';
        end
    elseif strcmp(VAR.ident, 'proxy')
        b = logical(1:n == VAR.shockpos);
        u_p = VARb.u(:,b);  % reduced-form residuals to be instrumented
        u_q = VARb.u(:,~b); % all other residuals
        nona = logical(isnan(zb) == 0);
        u_p = u_p(nona,1);
        u_q = u_q(nona,:);
        zb = zb(nona,1);
        % 1st-stage of IV: u_p on z
        zb = [ones(size(zb,1),1), zb];
        beta1 = zb\u_p;
        u_p_hat = zb*beta1;
        % 2nd-stage
        sb(~b) = u_p_hat\u_q;
        sb(b) = 1; % normalize s at shockpos to 1
        % impulse responses
        for hh=1:h
            IRFb(hh,:,bb) = (VARb.C(:,:,hh)*sb')';
        end    
    end

    % progress bar
    fprintf(1,'\b\b\b\b%3.0f%%',100*(bb/nboot)); 

end
fprintf('\n');

% retrieve intervals
IRFbands=prctile(IRFb, [lo up], 3);
