% Group linear model parameters
load('./Results/Group_rfMRI_resting_Q_G_lambda_complex.mat');

% Dynamics of the idx-th spatiotemporal mode 
% idx = 1,2,3 in our study: three slowest modes; 
% t = 200 TRs : temporal duration of mode dynamics
inCfg = struct('idx', idx, 't', 200,'lambda',lambda,'G',G);
[mode_dynamics,frequency] = calc_spatiotemporal_modes_dynamics(inCfg);

% Data reconstuction using modes (structural eigenmodes or spatiotemporal modes)
data_path = 'your_data_path';
file_name = 'your_file_name';

% spatiotemporal modes
phase_lambda = angle(lambda)/(2*pi*0.72); % 1TR = 0.72s
phase_lambda = sortrows([(1:length(lambda))',phase_lambda],2);
spatiotemporal_modes = G(:,phase_lambda(:,1));
inCfg01 = struct('data_path', data_path, 'file_name', file_name,'mode',spatiotemporal_modes);
corr_spatiotemporal_mode  = calc_recon_with_modes(inCfg01);

% structural eigenmodes
SC_path = 'your_path_of_SC_matrix'; 
load(SC_path);
[eigenmodes,~] = calc_structural_eigenmodes(SC);
inCfg02 = struct('data_path', data_path, 'file_name', file_name,'mode',eigenmodes);
corr_structural_eigenmode = calc_recon_with_modes(inCfg02);