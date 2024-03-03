function C = dyn_multipliers(VAR, h)

% inputs from VAR structure needed:
% A: estimated coefficients
% shockpos: position of identified shock in Cholesky ordering

c_case = VAR.c_case;
p = VAR.p;
n = VAR.n;
A = VAR.A;

AR = A(c_case+1:end,:)';            % AR coefficient [A1,A2,...,Ap]

AR_3d=NaN*ones(n,n,p);
for pp=1:1:p
    AR_3d(:,:,pp)=AR(:,(pp-1)*n+1:pp*n);
end


for i=1:p
    C(:,:,i) = zeros(n); 
end
C(:,:,p+1) = eye(n);

for t=p+2:p+h
    acc = 0;
    for j=1:p
        acc = acc + AR_3d(:,:,j)*C(:,:,t-j);
    end
    C(:,:,t) = acc;
end

% reindicize
C = C(:,:,p+1:end);

VAR = VAR;
