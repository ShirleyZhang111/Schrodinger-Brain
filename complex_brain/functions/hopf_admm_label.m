function [x,Z,res,flag0] = hopf_admm_label(A,B,C,mu,maxit,tol,label)
%HOPF_ADMM solves
%   min sum_j 0.5*norm(A(:,:,j)*x + Z*B(:,j) - C(:,j))^2 + mu*norm(Z,ell_1)
%Input:
%   A: space * polynomial basis * time
%   B: space * time
%   C: space * time
%   mu: penalty for sparsity, larger mu implies sparser Z
%   maxit: maximum number of iterations
%   tol: tolerance in stopping rule
%Output:
%   x: coefficients of polynomial basis
%   C: coupling matrix
%   res: residuals

% constants 
% space:n polynomial basis:m time:t
l = max(label(:,2));
[n,m,t] = size(A);
sqrtnt = sqrt(n*n);
sqrtm = sqrt(m);
% initialization
x = zeros(m*n,1);
Y = zeros(n,n);
Z = zeros(n,n);
lambda = zeros(n,n);
rho = min(max(mu*10,1),1000);
tau = 1e-3;
A = permute(A,[3,1,2]); % A:time*space*m
F = B*B' + eye(n).*tau;
res = zeros(2,maxit);
flag = 1;
flag0 = 1;
% ADMM iteration
for iter = 1:maxit
    x0 = x;
    Z0 = Z;
    if flag
        dF = decomposition(F + eye(n).*rho,'chol');
        rhosigma = rho*1.618;
        rhotau = 1/(rho+tau);
        murhotau = mu*rhotau;
    end
    U = (C-Y*B).'; % time*space
    W = zeros(size(U));
    for k = 1:l
        % update x
        s1 = find(label(:,2)==k);
        A1 = reshape(A(:,s1,:),[t*length(s1),m]); % (time*label)*m
        U1 = reshape(U(:,s1),[t*length(s1),1]); % (time*label)*1
        G = A1'*A1 + eye(m).*tau;
        dG = decomposition(G,'chol');
        x(1+(k-1)*m:k*m,1) = dG\(x(1+(k-1)*m:k*m,1).*tau + A1'*U1);
        % update W
        W(:,s1) = reshape(A1*x(1+(k-1)*m:k*m,1),[t,length(s1)]);
    end
    % update Y
    Y = (Y.*tau + Z.*rho - lambda + (C-W.')*B')/dF;
    % update Z
    Z = (Y.*rho + Z.*tau + lambda).*rhotau;
    Z = sign(Z).*max(abs(Z)-murhotau,0);
    % update lambda
    lambda = lambda + (Y-Z).*rhosigma;
    % check convergence
    res(1,iter) = norm(Y-Z,'fro')/sqrtnt;
    res(2,iter) = norm(x-x0,2)/sqrtm;
    res(3,iter) = norm(Z-Z0,'fro')/sqrtnt;
    if max(res(:,iter)) < tol
        flag0 = 0;
        break
    end
    
    % adapt penalty
    if res(1,iter) < mean(res(2:3,iter))*0.1
        rho = rho*0.618;
        flag = 1;
    elseif res(1,iter) > mean(res(2:3,iter))*10
        rho = rho/0.618;
        flag = 1;
    else
        flag = 0;
    end
end

end

