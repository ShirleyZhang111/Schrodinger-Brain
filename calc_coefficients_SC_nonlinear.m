function [coefficient, g] = calc_coefficients_SC_nonlinear(inCfg)
% Individual nonlinear model
% Input: 
%   TC : space * time
%   T : time points for fitting
%   SC : space * space, the structural connectivity matrix
%   field : complex field or real field
% Output:
%   coefficient: coefficients of polynomial basis
%   g: the global coupling parameter

% input configuration
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end % Time
if isfield(inCfg,'SC');SC = inCfg.SC;else; SC = [];end % the structural connectivity matrix
if isfield(inCfg,'field');field = inCfg.field;else;field = complex;end 
TC = inCfg.TC;

if strcmp(field, 'complex') == 1
    % Complex-valued Structure-informed model
    hilbert_data = double(hilbert(zscore(TC'))).';clearvars TC    
    B = hilbert_data(:,T);
    C = conv2(hilbert_data(:,[T(1:4)-4 T T(end-3:end)+4]),fliplr([1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280]),'valid');%.*-1i
    hilbert_data = permute(hilbert_data(:,T),[1,3,2]);
    A = [ones(size(hilbert_data)),...
        hilbert_data,conj(hilbert_data),...
        hilbert_data.*hilbert_data,hilbert_data.*conj(hilbert_data),conj(hilbert_data.*hilbert_data),...
        hilbert_data.*hilbert_data.*hilbert_data,hilbert_data.*hilbert_data.*conj(hilbert_data),hilbert_data.*conj(hilbert_data.*hilbert_data),conj(hilbert_data.*hilbert_data.*hilbert_data)];
    clearvars hilbert_data
elseif strcmp(field, 'real') == 1
    % Real-valued Structure-informed model
    zscore_data = double(zscore(TC')).';clearvars TC
    B = zscore_data(:,T);
    C = conv2(zscore_data(:,[T(1:4)-4 T T(end-3:end)+4]),fliplr([1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280]),'valid');%.*-1i
    zscore_data = permute(zscore_data(:,T),[1,3,2]);
    A = [ones(size(zscore_data)),...
        zscore_data,...
        zscore_data.*zscore_data,...
        zscore_data.*zscore_data.*zscore_data];
    clearvars zscore_data
end

% calculate coefficient through least square method
[n,m,t] = size(A);
A1 = reshape(permute(A,[1,3,2]),[n*t,m]);
A2 = reshape(SC*B,[n*t,1]);
A0 = [A1,A2];
dG = decomposition(A0'*A0,'chol');
coe = dG\(A0'*reshape(C,[],1));
coefficient = coe(1:m);
g = coe(m+1);