% linear model parameter: Q0 (voxel-level, unitary matrix)
data_path = 'your_save_path_of_linear_group_model'; % save_path from calc_coefficients_linear_group.m
load([date_path 'Group_rfMRI_resting_Q_G_lambda_complex.mat']);

node_items = length(lambda);

% stimulus to V1_L -> region-th (HCP360 atlas: region = 181)
% list contains the voxels belonging to the region V1_L.
% interval = 0.01TR, lambda0 = lambda.^(1/100)
pre_angle = angle(lambda) + 2*pi*(angle(lambda)<0);
complex_angle = pre_angle/100;
lambda0 = cos(complex_angle) + 1i*sin(complex_angle);

z0 = zeros(node_items,1);
z0(region,1) = 100; % stimulus = 100
S_response = zeros(node_items,2000); % temporal 20TRs
S_response(:,1) = z0;
y = G'*z0;
for kk = 2:2000
    y = diag(lambda0)*y;
    S_response(:,kk) = y;
end
S_response = G*S_response; % Response Signals

% Intrinsic timescale of each node
Intrinsic_timescale = zeros(node_items,1);
for k = 1:node_items
    [c,lags] = xcorr(abs(S_response(k,:)),'normalized');
    c = c(2000:end);
    [~,index] = find(0.6 <= c & c <= 0.62);
    Intrinsic_timescale(k,1) = index(1)/100;
end