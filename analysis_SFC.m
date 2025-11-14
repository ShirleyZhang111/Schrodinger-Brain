% calculate SFC 
load('./Data/SC_voxel.mat');
[node_items,~] = size(SC);

% SFC of FC (group level)
data_path = 'your_data_path';
items = dir(data_path);
items = items(3:end);
subjects_number = length(items);
FC = zeros(node_items,node_items,subjects_number);
for id = 1:subjects_number
    item = items(id).name;
    TC = load(strcat(data_path, '/', item));    
    FC(:,:,id) = corrcoef(TC');
end
FC = mean(FC,3); % group FC
SFC_FC = calc_SFC(FC,SC);

% SFC of H_complex
load('./Results/Group_rfMRI_resting_Q_G_lambda_complex.mat');
H_complex = abs(logm(G));
SFC_complex = calc_SFC(H_complex,SC); 

% SFC of H_real
load('./Results/Group_rfMRI_resting_Q_G_lambda_real.mat');
H_real = abs(logm(G));
SFC_real = calc_SFC(H_real,SC); 