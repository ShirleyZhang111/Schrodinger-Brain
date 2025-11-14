function [eigenmodes,lambda] = calc_structural_eigenmodes(SC)
% compute eigenmodes
% Input: 
%     SC : structure connectivity matrix
% Output: 
%     eigenmode: all the eigenmodes
%     lambda: the eigenvalue of Laplacian
L1 = diag(sum(SC,2))-SC;
[~,D1] = eig(L1);
L = L1./max(diag(D1));
[E,D] = eig(L);
lambda = diag(D);

lambda = [(1:length(lambda))',lambda];
lambda = sortrows(lambda,2);
eigenmodes = E(:,lambda(:,1));
