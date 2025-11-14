%% Linear Schrodinger-like model
% This code implements the linear Schrodinger-like model for fMRI data analysis

% load voxel-level rs-fMRI data
load('./Data/TC_HCP_9982.mat');

% Set model configuration parameters
field = 'complex';
inCfg = struct('field',field, 'T', 101:400, 'TC',TC); % Configuration structure: use time points 101-400 to compute model parameters

% Calculate parameters for linear Schrodinger-like model
[Q,G,lambda] = calc_coefficients_linear(inCfg);

% Time series prediction
% Prediction length: 25 TRs, Initial state: 401st TR
cfg = struct('field',field, 'T', 101:1000, 'Q',Q,...
    'Initial_time', 401-101, 'signal_length',25,'TC',TC);

% Execute prediction and get results
[data_predicted,data_obeserved,correlation] = calc_prediction(cfg);

% Analyze 25 TR prediction results
figure; 
plot(data_obeserved(1,:)); hold on;
plot(data_predicted(1,:)); hold on;
legend('observed','predicted');
figure;
boxplot(correlation);