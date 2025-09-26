function test_para_nonlinear_voxel(i_batch,batch_size,sub_index)
% clear;clc

addpath(genpath('functions'));

cfg = [];

% cfg.Project = 'HCP-D';
% cfg.Resolution = '8004+sub4008';
% cfg.Type = 'rfMRI';
% cfg.Condition = 'REST1';
% cfg.PE = 'AP';

% cfg.Project = 'HCP-YA';
% cfg.Resolution = '8004+sub4008';
% cfg.Type = 'rfMRI';
% cfg.Condition = 'REST1';
% cfg.PE = 'LR';

cfg.Project = 'HCP-A';
cfg.Resolution = '8004+sub4008';
cfg.Type = 'rfMRI';
cfg.Condition = [];
cfg.PE = [];
cfg.is_largescale = 1;

hasPrep = 1;
hasCalc = 0;

if ~isempty(cfg.Condition) % Young Adults and Development
    cfg.Experiment = [cfg.Type '_' cfg.Condition '_' cfg.PE];
else % Aging
    cfg.Experiment = cfg.Type;
end

if strcmp(cfg.Project,'HCP-YA')
    cfg.data_path = ['results/' cfg.Project '_' cfg.Experiment '_' cfg.Resolution '_T_400/'];
else
    cfg.data_path = ['results/' cfg.Project '_' cfg.Experiment '_' cfg.Resolution '/'];
end
if ~exist(cfg.data_path,'dir');mkdir(cfg.data_path);end

if ~hasPrep
    subjectList = generate_TC(cfg);
end

if ~hasCalc
    cfg.save_path = [cfg.data_path 'nonlinear_coefficients/'];
    if ~exist(cfg.save_path,'dir');mkdir(cfg.save_path);end
    
    
    if ~exist('subjectList', 'var')
        Files = dir(fullfile(cfg.data_path, '*.mat'));
        n_subs = length(Files);
        subjectList = cell(n_subs,1);
        for i_sub = 1:n_subs
            subjectList{i_sub} = Files(i_sub).name(1:end-4);
        end
    end
    
    if nargin<3
        sub_index = 1:5;%ceil(rand(1,10)*length(subjectList));
    end
    
    mu = batch_size*(i_batch - 1) + 1:min(100, batch_size*i_batch);
    if ~exist([cfg.save_path 'mu_test_mse/'],'dir');mkdir([cfg.save_path 'mu_test_mse/']);end
    
    for i_sub = sub_index
        disp(['i_sub = ' num2str(i_sub)]);
        if ~exist([cfg.save_path 'mu_test_mse/' subjectList{i_sub}],'dir');mkdir([cfg.save_path 'mu_test_mse/' subjectList{i_sub}]);end
        if ~exist([cfg.save_path 'mu_test_mse/' subjectList{i_sub} '/mu_' num2str(mu(1)) '_' num2str(mu(end)) '.mat'],'file')
            [~,~,~,mse,mse_real] = calc_coefficients_nonlinear(cfg.data_path, subjectList{i_sub}, cfg.save_path, mu,cfg.is_largescale);
            save([cfg.save_path 'mu_test_mse/' subjectList{i_sub} '/mu_' num2str(mu(1)) '_' num2str(mu(end)) '.mat'],'mse','mse_real');
        end
    end
    
    
end
