function main_nonlinear_voxel_HCP_EP(i_batch,batch_size)
% clear;clc
n_subs = 180;sub_index = batch_size*(i_batch - 1) + 1:min(n_subs, batch_size*i_batch);
cfg.sub_index = sub_index;

addpath(genpath('functions'));

hasPrep = 0;
hasCalc = 1;
Project = 'HCP_EP';
Resolution = '2mm';
cfg.data_path = '/public/mig_old_storage/home1/ISTBI_data/HCP/HCP_EP_1.0/Preprocessed/fmriprep/';
cfg.save_path = ['results/' Project '_voxel/'];if ~exist(cfg.save_path,'dir'); mkdir(cfg.save_path); end
cfg.is_largescale = 1;

load(['data/templates/label_MMP360_BN_shen_' Resolution '.mat'],'label');cfg.label = label;
% label = load_nii('data/templates/shen_3mm_268_parcellation.nii');cfg.label = label.img;
% label = load_nii('data/templates/HCPex_v1.0/HCPex_2mm.nii');cfg.label = label.img;

if ~exist(cfg.save_path,'dir');mkdir(cfg.save_path);end

if ~hasPrep
    eval(['subjectList = generate_TC_' Project '(cfg);']);
end

if ~hasCalc
    cfg.data_path = cfg.save_path;
    cfg.save_path = [cfg.save_path 'nonlinear_coefficients/'];if ~exist(cfg.save_path,'dir');mkdir(cfg.save_path);end
    
    if ~exist('subjectList','var')
        Files = dir(fullfile(cfg.data_path, '*mat'));
        n_subs = length(Files);
        subjectList = cell(n_subs,1);
        for i_sub = 1:n_subs
            subjectList{i_sub} = Files(i_sub).name;
        end
    end

    mu = 10:5:60;%0.12*(410-8)/5*4;
    % sub_index = batch_size*(i_batch - 1) + 1:min(n_subs, batch_size*i_batch);

    cfg.maxit = 500;
    cfg.tol = 1e-4;
    cfg.mu = mu;
    cfg.field = 'complex';

    for i_sub = sub_index%n_subs:-1:1 % sub_index
        disp(['i_sub = ' num2str(i_sub)]);
        cfg.file_name = subjectList{i_sub};
        [~,~,~,~,mse_real(:,i_sub),~] = calc_coefficients_nonlinear(cfg);
        % calc_coefficients_nonlinear(cfg);
    end
    cfg.save_path = [cfg.save_path 'test_mu/'];if ~exist(cfg.save_path,'dir');mkdir(cfg.save_path);end
    save([cfg.save_path 'mse_real_' num2str(sub_index) '.mat'],'mse_real');
end

