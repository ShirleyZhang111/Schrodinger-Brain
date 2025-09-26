function [Q,G,lambda] = calc_coefficients_linear(inCfg)
% Individual linear model with orthogonal constraint
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
data_path = inCfg.data_path;
file_name = inCfg.file_name;
save_path = inCfg.save_path;
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end 
if isfield(inCfg,'field');field = inCfg.field;else;field = 'complex';end 
if ~contains(file_name,'.mat');file_name = [file_name '.mat'];end

% start processing
load([data_path file_name],'TC');% N*T
TC = TC(:,6:end-5);
if isempty(T);T = 5:size(TC,2)-4;end

if strcmp(field, 'complex') == 1
    % Complex-valued linear model
    hilbert_data = double(hilbert(zscore(TC'))).'; clearvars TC
    data = hilbert_data(:,T); clearvars hilbert_data
elseif strcmp(field, 'real') == 1
    % Real-valued linear model
    zscore_data = double(zscore(TC')).';clearvars TC
    data = zscore_data(:,T); clearvars zscore_data
end

%  linear model with orthgonality  
A = data(:,1:end-1);
B = data(:,2:end);
[m,~] = size(A);
[U_r,S_r,V_r] = svd([A,B],'econ');
V1 = V_r(1:T-1,:);
V2 = V_r(T:2*(T-2),:);
[U0,~,V0] = svd(S_r*V2'*V1*S_r');
Q = eye(m)-U_r*(eye(2*(T-2))-U0*V0')*U_r';
[G,D] = eig(Q);
lambda = diag(D);
save([save_path file_name(4:end-4) '_linear_fitting_' field '.mat'], 'Q','G','lambda');