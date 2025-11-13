% data prediction
parameters_path = 'your_save_path_of_individual_linear_model_parameters'; % save path in calc_coefficients_linear.m
file_name = 'your_file_name'; % file name in calc_coefficients_linear.m
data_path = 'your_data_path'; % data path in calc_coefficients_linear.m
field = 'complex'; % or real
load([parameters_path file_name(4:end-4) '_linear_fitting_' field '.mat']);

% prediction length: 25TRs 
inCfg = struct('data_path', data_path, 'file_name', file_name, 'Q', Q,...
    'field',field, 'T', 100:1000 ,'Initial_time', 300,'signal_length',25);

[data_predicted,data_obeserved,corr] = calc_prediction(inCfg);

% Analysis of 25TRs prediction
figure; 
plot(data_obeserved(1,:)); hold on;
plot(data_predicted(1,:)); hold on;
legend('observed','predicted');
