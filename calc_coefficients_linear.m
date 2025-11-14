function [Q,G,lambda] = calc_coefficients_linear(inCfg)
% Individual-level linear model with orthogonal constraint
%     min_{Q*Q = I} sum_j norm(Qz(:,j) - z(:,j+1))^2
% Input: 
%     z : space * time
%     T : time points for fitting
%     field : complex or real field
% Output:
%     Q : the unitary matrix satisfies z(t+1) = Q z(t)
%     G : the eigenvectors of Q 
%     lambda : the eigenvalues of Q

% input configuration
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end 
if isfield(inCfg,'field');field = inCfg.field;else;field = 'complex';end 

TC = inCfg.TC;

% data z-score
if strcmp(field, 'complex') == 1
    % Complex-valued data
    hilbert_data = double(hilbert(zscore(TC'))).';
    data = hilbert_data(:,T); clearvars hilbert_data
elseif strcmp(field, 'real') == 1
    % Real-valued data
    zscore_data = double(zscore(TC')).';
    data = zscore_data(:,T); clearvars zscore_data
end

%  Compute model parameters with orthgonality  
T0 = length(T);
A = data(:,1:end-1);
B = data(:,2:end);
[m,~] = size(A);
[U_r,S_r,V_r] = svd([A,B],'econ');
V1 = V_r(1:T0-1,:);
V2 = V_r(T0:2*T0-2,:);
[U0,~,V0] = svd(S_r*V2'*V1*S_r');
Q = eye(m)-U_r*(eye(2*T0-2)-U0*V0')*U_r';
[G,D] = eig(Q);
lambda = diag(D);
% save([save_path file_name(4:end-4) '_linear_fitting_' field '.mat'], 'Q','G','lambda');