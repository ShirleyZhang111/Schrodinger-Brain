function main_nonlinear_region_SCZ(i_batch,batch_size)
%clear;clc

addpath(genpath('functions'));



hasPrep = 1;
hasCalc = 0;

cfg.data_path = '/public/home/cyqian/SCZ/fMRI_Jingzong_cifti_new/fmriprep/';
cfg.save_path = 'results/SCZ_Jingzong/';
load('data/templates/label_MMP360_BN_shen_2mm.mat','label');
cfg.label = label;
if ~exist(cfg.save_path,'dir');mkdir(cfg.save_path);end

if ~hasPrep
    subjectList = generate_TC_SCZ(cfg);
end

if ~hasCalc
    cfg.data_path = cfg.save_path;
    cfg.save_path = [cfg.save_path 'nonlinear_coefficients/'];
    if ~exist(cfg.save_path,'dir');mkdir(cfg.save_path);end
    
    
    if ~exist('subjectList', 'var')
        Files = dir(fullfile(cfg.data_path, '*.mat'));
        n_subs = length(Files);
        subjectList = cell(n_subs,1);
        for i_sub = 1:n_subs
            subjectList{i_sub} = Files(i_sub).name;
        end
    end
    

    mu = 0.9; %[0:0.2:2 3:10];
    sub_index = batch_size*(i_batch - 1) + 1:min(n_subs, batch_size*i_batch);

    cfg.is_largescale = 0;
    cfg.maxit = 1000;
    cfg.tol = 1e-12;
    cfg.mu = mu;
    cfg.feild = 'complex';

    for i_sub = sub_index
        disp(['i_sub = ' num2str(i_sub)]);
        cfg.file_name = subjectList{i_sub};
        %[~,~,~,~,mse_real(:,i_sub),~] = calc_coefficients_nonlinear(cfg);
        calc_coefficients_nonlinear(cfg);
    end
    
   %save(['results/SCZ_Jingzong/test_mu/mse_real_' num2str(sub_index(1)) '_' num2str(sub_index(end)) '.mat'],'mse_real');
end
