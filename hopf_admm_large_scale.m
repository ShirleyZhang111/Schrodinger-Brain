function [x,Z,res,isconverge] = hopf_admm_large_scale(A,B,C,mu,maxit,tol)
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
%by Weiyang Ding @Fudan December 15, 2024

% constants
[n,m,t] = size(A);
% initialization
x = zeros(m,1,'single');
Y = zeros(n,n,'single');
Z = zeros(n,n,'single');
lambda = zeros(n,n,'single');
rho = min(max(mu*10,1),1000);
tau = single(1e-3);
A = reshape(permute(A,[1,3,2]),[n*t,m]);
G = A'*A + eye(m,'single').*tau;
dG = decomposition(G,'chol');
F = B'*B + eye(t,'single').*tau;
res = zeros(3,maxit,'single');
flag = 1;
isconverge = 0;
% ADMM iteration
for iter = 1:maxit
    tic
    x0 = x;
    Z0 = Z;
    if flag
        R = chol(F + eye(t,'single').*rho);
        BR = B/R;
        rhosigma = rho*1.618;
        rhotau = 1/(rho+tau);
        murhotau = mu*rhotau;
    end
    x = dG\(x.*tau + A'*reshape(C-Y*B,[],1));
    Y = (Y.*tau + Z.*rho - lambda + (C-reshape(A*x,[n,t]))*B').*rhotau;
    Y = Y - (Y*BR)*BR';
    Z = (Y.*rho + Z.*tau + lambda).*rhotau;
    Z = sign(Z).*max(abs(Z)-murhotau,0);
    lambda = lambda + (Y-Z).*rhosigma;
    % check convergence
    nZ = norm(Z,'fro');
    res(1,iter) = norm(Y-Z,'fro')/nZ;
    res(2,iter) = norm(x-x0,2)/norm(x,2);
    res(3,iter) = norm(Z-Z0,'fro')/nZ;
    if max(res(:,iter)) < tol
        isconverge = 1;
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
    disp(['Iteration ',num2str(iter),', use time: ' num2str(toc/60) ' mins, residual: ',num2str(max(res(:,iter)))])
end

end

