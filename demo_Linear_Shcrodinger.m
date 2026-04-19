%% Linear Schrodinger-like model
% This code implements the linear Schrodinger-like model for fMRI data analysis (prediction for voxel-level data)

% load voxel-level rs-fMRI data
load('./Data/TC_HCP_9982.mat');

% Set model configuration parameters
field = 'complex';
inCfg = struct('field',field, 'T', 101:500, 'TC',TC); % Configuration structure: use time points 101-500 to compute model parameters

% Calculate parameters for linear Schrodinger-like model
[Q,G,lambda] = calc_coefficients_linear(inCfg);

% Time series prediction
% Prediction length: 20 TRs, Initial state: 501st TR
cfg = struct('field',field, 'T', 101:1000, 'Q',Q,...
    'Initial_time', 501-101, 'signal_length',20,'TC',TC); % initial 501, for time points 101-500, initial = 501-101

% Execute prediction and get results
[data_predicted,data_obeserved,correlation] = calc_prediction(cfg);

% Analyze 20 TR prediction results
figure; 
plot(data_obeserved(1,:)); hold on;
plot(data_predicted(1,:)); hold on;
legend('observed','predicted');
figure;
boxplot(correlation);
