function [x,Z,res,flag0] = hopf_admm(A,B,C,mu,maxit,tol)
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
%by Weiyang Ding @Fudan November 22, 2023

% constants
[n,m,t] = size(A);
sqrtnt = sqrt(n*t);
sqrtm = sqrt(m);
% initialization
x = zeros(m,1);
Y = zeros(n,n);
Z = zeros(n,n);
lambda = zeros(n,n);
rho = min(max(mu*10,1),1000);
tau = 1e-3;
A = reshape(permute(A,[1,3,2]),[n*t,m]);
G = A'*A + eye(m).*tau;
dG = decomposition(G,'chol');
% F = B*B' + eye(n).*tau;
F = B'*B + eye(t).*tau;% modified
res = zeros(2,maxit);
flag = 1;
flag0 = 1;
% ADMM iteration
for iter = 1:maxit
    %disp(['iteration = ' num2str(iter)]);
    x0 = x;
    Z0 = Z;
    if flag
        dF = decomposition(F + eye(t).*rho,'chol');
        rhosigma = rho*1.618;
        rhotau = 1/(rho+tau);
        murhotau = mu*rhotau;
    end

    % update coefficient x
    x = dG\(x.*tau + A'*reshape(C-Y*B,[],1));

    % update Y
%     Y = (Y.*tau + Z.*rho - lambda + (C-reshape(A*x,[n,t]))*B')/dF;
    Y = (Y.*tau + Z.*rho - lambda + (C-reshape(A*x,[n,t]))*B')./(tau+rho);% modified
    Y = Y - Y*B/dF*B';% modified
    
    % update Z
    Z = (Y.*rho + Z.*tau + lambda).*rhotau;
    Z = sign(Z).*max(abs(Z)-murhotau,0);
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

