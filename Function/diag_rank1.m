function [d,u,v] = diag_rank1(A0,maxit,tol)
%DIAG_RANK1 decomposes A0 = diag(d) + u*v'

A = A0;
d = diag(A0);
res = zeros(2,maxit);
for iter = 1:maxit
    dold = d;
    Aold = A;
    [U,S,V] = svd(A);
    u = U(:,1);
    v = V(:,1).*S(1,1);
    A = u*v';
    d = diag(A0-A);
    A = A0 - diag(d);
    res(1,iter) = norm(d-dold);
    res(2,iter) = norm(A-Aold,'fro');
    if norm(res(:,iter)) < tol
        break
    end
end

figure;
semilogy(res(1,:)); hold on
semilogy(res(2,:)); hold off
legend('diagonal part','rank-1 part');
end

