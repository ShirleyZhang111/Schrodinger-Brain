function [data_predicted,data_obeserved,correlation] = calc_prediction(inCfg)
% Data prediction using linear schrodinger model/real counterpart
% Input: 
%     Initial_time : time of the initial state
%     signal_length : the predicted/reconstructed siganl length
%     T : time points for fitting
%     field : complex or real field
%     TC : space * time
% Output:
%     data_predicted : the predicted data
%     data_observed : the empirical data
%     corr : the peason's correlation of the predicted data and the
%     empirical data 


% input configuration
% data_path = inCfg.data_path;
% file_name = inCfg.file_name;
Q = inCfg.Q; 

if isfield(inCfg,'T');T = inCfg.T;else;T = [];end % Time points
if isfield(inCfg,'Initial_time');Initial_time = inCfg.Initial_time;else;Initial_time = 1;end % Initial state
if isfield(inCfg,'signal_length');signal_length = inCfg.signal_length;else;signal_length = 1;end % prediction length
if isfield(inCfg,'field');field = inCfg.field;else;field = 'complex';end 
% if ~contains(file_name,'.mat');file_name = [file_name '.mat'];end


TC = inCfg.TC;

if strcmp(field, 'complex') == 1
    % Complex-valued linear model
    hilbert_data = double(hilbert(zscore(TC'))).'; 
    data = hilbert_data(:,T); clearvars hilbert_data
elseif strcmp(field, 'real') == 1
    % Real-valued linear model
    zscore_data = double(zscore(TC')).';
    data = zscore_data(:,T); clearvars zscore_data
end

% data prediction with initial state 'Initial_time'
S0 = real(data);
z0 = data(:,Initial_time);
data_predicted = zeros(length(data(:,1)),signal_length);
data_predicted(:,1) = real(z0);
for kk = 2:signal_length
    z0 = Q*z0;
    data_predicted(:,kk) = real(z0)';
end      
correlation = zeros(length(data(:,1)),1);
data_obeserved = S0(:,Initial_time:Initial_time-1+signal_length);
for k = 1:length(data(:,1))
    [r,~] = corr(data_obeserved(k,:)',data_predicted(k,:)');
    correlation(k,1) = r;
end