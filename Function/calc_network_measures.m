function [Num,Amp_m,G,G0,density,r,Q] = calc_network_measures(H,list)
% Network measures calculation
% Input: 
%     H : Network connectivity matrix
%     list : Module distribution of the network
% Output:
%     Num：number of coupling nonzeros
%     Amp_m: Mean amplitude of Nonzeros
%     G : Network communicability
%     G0 : Communicability matrix
%     density : Network density
%     r : Assortativity
%     Q : Modularity Index

H = H - diag(diag(H));

% Number of Nonzeros
Num = nnz(H);

% Mean amplitude of Nonzeros
Amp_m = sum(sum(abs(H)))/Num;

% Communicability
d = sum(abs(H),2);
W = diag(d.^(-1/2))*H*diag(d.^(-1/2));
[E0,D0] = eig(W);
G0 = E0*diag(exp(diag(D0)))*E0';
G = mean(mean(abs(G0)));

[len,~] = size(H);
% Network Density
E_dir = sum(sum(abs(H)));
density = E_dir/len/(len-1);

% Assortativity
k = sum(abs(H),2);
r1 = 0; r2 = 0; r3 = 0; 
for i = 1:379
    for j = i:379
        r1 = r1 + k(i)*k(j)*A(i,j);
        r2 = r2 + 0.5*(k(i)+k(j))*A(i,j);
        r3 = r3 + 0.5*(k(i)^2+k(j)^2)*A(i,j);
    end
end
r = (r1/E_dir - (r2/E_dir)^2)/(r3/E_dir- (r2/E_dir)^2);

% Modularity
delta = zeros(len,len);
for i = 1:len
    for j = 1:len
        delta(i,j) = list(i,2) == list(j,2);
    end
end
delta = double(delta);
q = zeros(len,len);
for i = 1:len
    for j = 1:len
        q(i,j) = (A(i,j)-k(i)*k(j)/E_dir)*delta(i,j);
    end
end
Q = sum(sum(q))/E_dir/2;