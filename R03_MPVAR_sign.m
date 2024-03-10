
%% Preamble
% This script replicates a small VAR in which a monetary policy shock is 
% identified using sign restrictions, similar to Uhlig (2005), but with 
% notable differences (e.g. restriction on GDP)

clear all; clc;

addpath('./_tbx')   % if script is run outside of the toolbox (recommended), make adjustments
addpath('./_tbx/supportfct') 
addpath('./_tbx/var_tbx') 


%% Settings
folder = 'R03_MPVAR_sign/';   
if ~exist('folder', 'dir')
    mkdir(folder)
end
p = 12;
h = 48;
c_case = 1;
vars = {'GDPgap','Unemp','CoreCPIGr12','FFR'};
exvars = {};
n = length(vars);
ident = 'sign';
shockpos = 1; % position of the shock in the Rest matrix columns
shocksize = 0; % 0 = one standard deviation, all else: absolute values
state.nonlinear = 'no';
nboot = 1000;
alpha = 90; %confidence level


%% Load data, generate lags, subset relevant subsample, etc.
load('data/data_m_ready.mat')
subset_data;


%% Estimate matrices A, Omega, S, dynamic multipliers
VAR = estimateVAR(data, p, c_case, []);
VAR.C = dyn_multipliers(VAR, h); % not identified


%% Identification: Sign restrictions
VAR.ident = ident;
% define restrictions // each column is one of q identified shocks 
% and each of n rows contains -1 and 1 for negative/positive restrictions
% of the row-variable; 
Rest = [-1,  +1, -1, NaN; 
        NaN, NaN, NaN, NaN; 
        -1,  +1, +1, NaN;
        1,   +1, +1, NaN];
if size(Rest,1) ~= n
    error('Number of rows in Rest matrix has to be equal to n.')
end
Rest = Rest(:, sum(~isnan(Rest), 1) > 0); % select only columns that actually contain identifying restrictions
Start = [5,   5,  5, NaN; 
        NaN, NaN, NaN, NaN; 
         5,   5,  5, NaN;
         0,   0,  0, NaN]; % first period for binding restriction
Start = Start(1:size(Rest,1),1:size(Rest,2));
End = 5*ones(size(Rest)); % final period for binding restriction
End = End(1:size(Rest,1),1:size(Rest,2));
nshocks = size(Rest, 2);
Rest_check=zeros(nshocks,1);
perc_accepted =NaN*ones(1,nboot);
Qb = NaN*ones(n,n,nboot);

% while-loop to look for nboot different Q's that satisfy the restrictions 
accepteddraws = 0;
attempts = 1;
disp('Progress \n')
while accepteddraws < nboot

    % draw a rotation of Q;
    draw=randn(n,n);
    [Q,R]=qr(draw);
    for nn=1:n;
        if R(nn,nn)<0
            Q(:,nn)=-Q(:,nn);
        end
    end

    ii = 0;
    while ii<nshocks  % loop over each of the shocks (columns in 'Rest')
        ii=ii+1;
        shock=zeros(n,1);
        shock(ii,1)=1;
        IRFcandidate = NaN*ones(max(max(End))+1,n);
        for hh=1:(max(max(End))+1)
            IRFcandidate(hh,:) = (VAR.C(:,:,hh)*VAR.S*Q*shock)';
        end

        junk=[1:n]';
        Which_Variab=junk(isfinite(Rest(:,ii)));
        Which_Sign=Rest(isfinite(Rest(:,ii)),ii);
        Which_Start=Start(isfinite(Rest(:,ii)),ii);
        Which_End=End(isfinite(Rest(:,ii)),ii);

        jj=1;
        while jj<=size(Which_Variab,1)
            a=1;
            rev_cond= 1; % check all restrictions the way the are written 
            a_temp = testIRFsign(IRFcandidate(Which_Start(jj,1)+1:Which_End(jj,1)+1,Which_Variab(jj,1)),Which_Sign(jj,1),rev_cond);
            a=a*a_temp; 
            if a_temp==1; 
                jj=jj+1;
            elseif a_temp==0; 
                jj=size(Which_Variab,1)+100; 
            end    
        end

        if a==0 % try to invert all the signs 
            jj=1; 
            while jj<=size(Which_Variab,1)
                a=1;
                rev_cond= 0; % check all restrictions the way the are written 
                a_temp = testIRFsign(IRFcandidate(Which_Start(jj,1)+1:Which_End(jj,1)+1,Which_Variab(jj,1)),Which_Sign(jj,1),rev_cond);
                a=a*a_temp; 
                if a==1; 
                    jj=jj+1; 
                elseif a==0; 
                    jj=size(Which_Variab,1)+100; 
                end    
            end
            if a==1; 
                Q(:,ii)=-Q(:,ii); 
            end
        end
           
        if a==0; 
            Rest_check(ii,1)=0; 
            ii=nshocks+100;   %The restrictions are not satisfied.  Go the beginning to redraw.
        else
            Rest_check(ii,1)=1;
        end

    end

    if mean(Rest_check) == 1
        % save accepted Q's and show progress bar
        accepteddraws = accepteddraws + 1;
        perc_accepted(accepteddraws) = (accepteddraws / attempts)*100;
        Qb(:,:,accepteddraws) = Q;
        fprintf(1,'\b\b\b\b%3.0f%%',100*(accepteddraws/nboot)); 
    end
    attempts = attempts + 1;
end
fprintf('\n');

% figure: share of accepted draws (in %)
figure()
plot(1:nboot, perc_accepted, 'k', 'LineWidth', 2)
grid on
title('Percent of draws accepted','FontSize',14)
set(gcf,'paperpositionmode','auto')
set(gcf, 'position', [0 0 600 200]);
xlabel('Draw')
print(gcf,'-dpng','-loose',strcat(folder,'/drawsaccepted.png'));

clear accepteddraws attempts Rest_check junk Which_End Which_Sign ...
    Which_Start Which_Variab ii jj a a_temp draw H R rev_cond perc_accepted ...
    nn IRFcandidate

% Produce IRFs with accepted Q's for each shock/column in 'Rest'
for s = 1:nshocks
    
    IRFb = NaN*ones(h,n,nboot);
    shock = zeros(n,1);
    shock(s) = 1;
    for bb = 1:nboot
        P = VAR.S*Qb(:,:,bb);
        for hh=1:h
            IRFb(hh,:,bb) = (VAR.C(:,:,hh)*P*shock)';
        end
    end
    IRF = mean(IRFb, 3);
    lo = (100-alpha)/2;
    up = 100 - lo;
    IRFbands = prctile(IRFb, [lo up],3);

    if s == shockpos
        VAR.shock = shock;
        VAR.Rest = Rest;
        VAR.Start = Start;
        VAR.End = End;
        VAR.Q = mean(Qb, 3);
        VAR.IRF = IRF;
        VAR.IRFbands = IRFbands;
    end
    
    % Plot
    plotirf1(IRF, IRFbands, printvars, strcat(folder, 'shock', num2str(s),'_lin'))
    
end


%% Historical decompositions
Acomp = [VAR.A(1+c_case:n*p+c_case,:)'; eye(n*(p-1)) zeros(n*(p-1),n)];
Icomp = [eye(n) zeros(n,(p-1)*n)];
% companion matrix: 3 x (3x4), lower rows filled with identity
t = VAR.t-p;
HD = zeros(t+p,n,n,nboot);
for bb = 1:nboot
    P = VAR.S*Qb(:,:,bb);
    eps = P\VAR.u'; % structural errors (nxT)
    P_big = zeros(n*p,n);
    P_big(1:n,:) = P;
    HDshock_big = zeros(p*n,t+1,n);
    HDshock = zeros(n,t+1,n);
    for j=1:n % for each variable
        eps_big = zeros(n,t+1); % matrix of shocks conformable with companion
        eps_big(j,2:end) = eps(j,:);
        for i = 2:t+1
            HDshock_big(:,i,j) = P_big*eps_big(:,i) + Acomp*HDshock_big(:,i-1,j);
            HDshock(:,i,j) =  Icomp*HDshock_big(:,i,j);
        end
    end
    hd = zeros(t+p,n,n);  % [nobs x shock x var]
    for i=1:n
        for j=1:n
            hd(:,j,i) = [nan(p,1); HDshock(i,2:end,j)'];
        end
    end
    HD(:,:,:,bb) = hd;
end
VAR.HD = mean(HD, 4);

plothd(VAR.HD, time, printvars, folder, {'Monetary policy','Demand','Cost-push'}, [])


% Organize
clear s bb hh Rest Start End shock shockpos Qb IRFb lo up P
clear Acomp Icomp P_big HDshock_big HDshock j eps_big i hd
save(strcat(folder, 'out.mat'), 'VAR')




