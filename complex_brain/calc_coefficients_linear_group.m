function [Q,G,lambda] = calc_coefficients_linear_group(inCfg)
% Group linear model with orthogonal constraint
%     min_{Q*Q = I} sum_s sum_j norm(Q reshape(z(:,j,s),[],1) - reshape(z(:,j+1,s),[],1))^2
% Input: 
%     z : space * time * subjects
%     T : time points for fitting
%     field : complex field or real field
% Output:
%     Q : the unitary matrix
%     G : the eigenvectors of Q 
%     lambda : the eigenvalues of Q


% input configuration
data_path = inCfg.data_path;
save_path = inCfg.save_path;
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end % Time points for fitting
if isfield(inCfg,'field');field = inCfg.field;else;field = 'complex';end 
if isfield(inCfg,'item');item = inCfg.item;else;item = 1;end % number of subjects
if isfield(inCfg,'idx');idx = inCfg.idx;else;idx = 1;end % space, number of node (region/voxel)

W = zeros(idx,idx);
for k = 1:item
    item = items(k).name;
    data = load(strcat(data_path, '/', item));
    TC = data.TC;
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
    W = W + B*A';
end

[U,~,V] = svd(W);
Q = U*V';
[G,D] = eig(Q);
lambda = diag(D);
save([save_path 'Group_rfMRI_resting_Q_G_lambda_' field '.mat'], 'Q','G','lambda');