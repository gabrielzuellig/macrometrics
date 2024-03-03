function [a]=testIRFsign(irf_vec,sign_to_check,rev_cond)
% 
% this functions checks whether the sign restriction on the kk variable is 
% verified 
% INPUT: 
% - irf_vec = hx1 vector of impulse response to the ith shock of the kth
%           variable (h = EndMAT(kk,ii)-StartMAT(kk,ii)); 
% - sign_to_check = +1 Positive 
%                   -1 Negative 
% - rev_cond = 1 (checks that conditions are verified as shown) 
%            = 0 (check whether the reverse are verified)
% 
% OUTPUT: 
% - a = 1 when sign restriction is verified (for the all h periods) 
%       0 otherwise

a=0;

if rev_cond==1
    if sign_to_check==-1;
    ggg=find(irf_vec<=0);
    elseif sign_to_check==1;
    ggg=find(irf_vec>=0);
    end
elseif rev_cond==0
    if sign_to_check==1;
    ggg=find(irf_vec<=0);
    elseif sign_to_check==-1;
    ggg=find(irf_vec>=0);
    end
end

if size(ggg,1)==size(irf_vec,1); 
    a=1; 
end