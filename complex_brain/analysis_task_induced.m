% Resting-state coupling matrix (Group-level)
data_path_resting = 'your_save_path_of_resting_state_nonlinear_model'; % save_path in calc_coefficients_nonlinear.m
items = dir(data_path_resting);
items = items(3:end);
subjects_number_resting = length(items);
C_r = [];
for id = 1:subjects_number_resting
    item = items(id).name;
    hopf = load(strcat(data_path_resting, '/', item));    
    C_r = C_r + hopf.coupling_mat;
end
C_r = C_r./subjects_number_resting;

% Task-state coupling matrix (Group-level)
data_path_task = 'your_save_path_of_task_state_nonlinear_model'; % save_path in calc_coefficients_nonlinear.m
items = dir(data_path_task);
items = items(3:end);
subjects_number_task = length(items);
C_t = [];
for id = 1:subjects_number_task
    item = items(id).name;
    hopf = load(strcat(data_path2, '/', item));    
    C_t = C_t + hopf.coupling_mat;
end
C_t = C_t./subjects_number_task;

% Relative coupling matrix R
R = C_t/C_r;

% Decomposition of relative coupling matrix R
maxit = 100;
tol = 1e-12;
[d,u,v] = diag_rank1(abs(R),maxit,tol);
R_approximation = diag(D) + u*v';