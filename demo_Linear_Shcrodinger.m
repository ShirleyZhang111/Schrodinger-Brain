%% Linear Schrodinger-like model
% This code implements the linear Schrodinger-like model for fMRI data analysis (prediction for voxel-level data)

% load voxel-level rs-fMRI data
load('./Data/TC_HCP_9982_time1.mat');
load('./Data/TC_HCP_9982_time2.mat');
TC = [TC1,TC2];

% Set model configuration parameters
field = 'complex';
inCfg = struct('field',field, 'T', 1:400, 'TC',TC); % Configuration structure: use time points 101-500 to compute model parameters

% Calculate parameters for linear Schrodinger-like model
Q = calc_coefficients_linear(inCfg);

% Time series prediction
% Prediction length: 10 TRs, Initial state: 400st TR
cfg = struct('field',field, 'T', 1:500, 'Q',Q,...
    'Initial_time', 400, 'signal_length',10,'TC',TC); % initial 400

% Execute prediction and get results
[data_predicted,data_obeserved,correlation] = calc_prediction(cfg);

% Analyze 10 TR prediction results
figure; 
plot(data_obeserved(1,:)); hold on;
plot(data_predicted(1,:)); hold on;
legend('observed','predicted');
figure;
boxplot(correlation);
