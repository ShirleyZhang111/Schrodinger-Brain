function main_nonlinear_voxel(i_batch,batch_size)
% clear;clc

addpath(genpath('functions'));

cfg = [];
%cfg.Project = 'HCP-D';
%cfg.Resolution = '8004+sub4008';
%cfg.Type = 'rfMRI';
%cfg.Condition = 'REST1';
%cfg.PE = 'AP';
cfg.is_largescale = 1;
%cfg.Project = 'HCP-YA';
%cfg.Resolution = '8004+sub4008';
%cfg.Type = 'rfMRI';
%cfg.Condition = 'REST1';
%cfg.PE = 'LR';

cfg.Project = 'HCP-A';
cfg.Resolution = '8004+sub4008';
cfg.Type = 'rfMRI';
cfg.Condition = [];
cfg.PE = [];

hasPrep = 1;
hasCalc = 0;

if ~isempty(cfg.Condition) % Young Adults and Development
    cfg.Experiment = [cfg.Type '_' cfg.Condition '_' cfg.PE];
else % Aging
    cfg.Experiment = cfg.Type;
end

cfg.data_path = ['results/' cfg.Project '_' cfg.Experiment '_' cfg.Resolution '/'];

if ~exist(cfg.data_path,'dir');mkdir(cfg.data_path);end

if ~hasPrep
    subjectList = generate_TC(cfg);
end

mu = 39;

if ~hasCalc
    cfg.save_path = [cfg.data_path 'nonlinear_coefficients_mu_' num2str(mu) '/'];
    if ~exist(cfg.save_path,'dir');mkdir(cfg.save_path);end
    
    
    if ~exist('subjectList', 'var')
        Files = dir(fullfile(cfg.data_path, '*.mat'));
        n_subs = length(Files);
        subjectList = cell(n_subs,1);
        for i_sub = 1:n_subs
            subjectList{i_sub} = Files(i_sub).name;
        end
    end
   
    sub_index = batch_size*(i_batch - 1) + 1:min(n_subs, batch_size*i_batch);
    
    for i_sub = sub_index
        disp(['i_sub = ' num2str(i_sub)]);
        calc_coefficients_nonlinear(cfg.data_path, subjectList{i_sub}, cfg.save_path, mu, cfg.is_largescale);
    end
    
end
