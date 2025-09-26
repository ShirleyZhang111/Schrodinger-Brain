function [coefficient, coupling_mat, res, mse, mse_real, energy] = calc_coefficients_nonlinear(inCfg)
% Individual nonlinear model
% Input: 
%   TC : space * time
%   T : time points for fitting
%   field : complex field or real field
% Output:
%   coefficient: coefficients of polynomial basis
%   coupling_mat: coupling matrix
%   res: residuals
%   mse/mse_real : mean squared error of data fitting
%   energy : total energy, the absolute number of hamiltonian of each time point


% input configuration
data_path = inCfg.data_path;
file_name = inCfg.file_name;
save_path = inCfg.save_path;
is_largescale = inCfg.is_largescale;
if isfield(inCfg,'maxit');maxit = inCfg.maxit;else;maxit = 500;end
if isfield(inCfg,'tol');tol = inCfg.tol;else;tol = 1e-4;end
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end
if isfield(inCfg,'mu');mu = inCfg.mu;else;mu = single(0.12*length(T));end
if isfield(inCfg,'field');field = inCfg.field;else;field = 'complex';end % complex field: field = 1; real field: field = 0;
if ~contains(file_name,'.mat');file_name = [file_name '.mat'];end

% start processing
load([data_path file_name],'TC');% N*T
TC = TC(:,6:end-5);
if isempty(T);T = 5:size(TC,2)-4;end

if strcmp(field, 'complex') == 1
    % Complex-valued nonlinear data-driven model
    hilbert_data = double(hilbert(zscore(TC'))).';clearvars TC
    B = hilbert_data(:,T);t0 = round(size(B,2)/5*4);
    C = conv2(hilbert_data(:,[T(1:4)-4 T T(end-3:end)+4]),fliplr([1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280]),'valid');%.*-1i
    % C = conv2(hilbert_data(:,[T(1)-1 T T(end)+1]),[1/2 0 -1/2],'valid');
    hilbert_data = permute(hilbert_data(:,T),[1,3,2]);
    A = [ones(size(hilbert_data)),...
        hilbert_data,conj(hilbert_data),...
        hilbert_data.*hilbert_data,hilbert_data.*conj(hilbert_data),conj(hilbert_data.*hilbert_data),...
        hilbert_data.*hilbert_data.*hilbert_data,hilbert_data.*hilbert_data.*conj(hilbert_data),hilbert_data.*conj(hilbert_data.*hilbert_data),conj(hilbert_data.*hilbert_data.*hilbert_data)];
    clearvars hilbert_data
elseif strcmp(field, 'real') == 1
    % Real-valued nonlinear data-driven model
    zscore_data = double(zscore(TC')).';clearvars TC
    B = zscore_data(:,T); t0 = round(size(B,2)/5*4);
    C = conv2(zscore_data(:,[T(1:4)-4 T T(end-3:end)+4]),fliplr([1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280]),'valid');%.*-1i
    zscore_data = permute(zscore_data(:,T),[1,3,2]);
    A = [ones(size(zscore_data)),...
        zscore_data,...
        zscore_data.*zscore_data,...
        zscore_data.*zscore_data.*zscore_data];
    clearvars zscore_data
end

mse = []; mse_real = []; energy = zeros(t0,1);
tic;
if length(mu)==1 % calculate results with optimal sparsity parameter
    if is_largescale
        A = single(A);B = single(B);C = single(C);
        [coefficient, coupling_mat, res, isconverge] = hopf_admm_large_scale(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,maxit,tol);
    else
        [coefficient, coupling_mat, res, isconverge] = hopf_admm(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,maxit,tol);
    end
    for time = 1:t0
        energy(time,1) = B(:,time)'*construct_data(B(:,time),coefficient,coupling_mat)*1i; 
    end
    if isconverge;n_iter = length(find(res(1,:)>0));else;n_iter = maxit;end
    disp(['Subject ',file_name(4:end-4),' finished, use time ',num2str(n_iter*toc/60),' mins, residual: ' num2str(max(res(:,n_iter)))]);
    % Save the result to a file
    save([save_path file_name(4:end-4) '_nonlinear_fitting_' field '.mat'], 'coefficient', 'coupling_mat', 'res','isconverge','mu','tol','energy');
    
else % find optimal sparsity parameter
    for i_mu = length(mu):-1:1
        if is_largescale
            A = single(A);B = single(B);C = single(C);
            [coefficient, coupling_mat, res, ~] = hopf_admm_large_scale(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu(i_mu),maxit,tol);
        else
            [coefficient, coupling_mat, res, ~] = hopf_admm(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu(i_mu),maxit,tol);
        end
        delta_Z = C(:,t0+1:end);
        dZ = construct_data(B(:,t0+1:end),coefficient,coupling_mat);
        mse(1,i_mu) = norm(delta_Z-dZ,'fro');
        mse_real(1,i_mu) = norm(real(delta_Z-dZ),'fro');
    end   
end