function extract_coupling_mat_prop(idx_project)
% clear;clc

Resolution = '360+sub19'; % voxel-level: '8004+sub4008';
Project = {'HCP-D','HCP-YA','HCP-A'};
Type = 'rfMRI';

load('label_12k_MMP360+sub19.mat','label');
n_regions = max(label,[],'all');
mask_withinRegion = single(sqrt(label*label'));
mask_withinRegion(mask_withinRegion~=round(mask_withinRegion)) = 0;
mask_withinRegion = mask_withinRegion - diag(diag(mask_withinRegion));
mask_withinRegion(mask_withinRegion==0) = nan;
label0 = mean(mask_withinRegion,'omitnan')';
idx_wrong = find(label0~=label);
for i_region = 1:length(idx_wrong)
    idx_i = idx_wrong(i_region);
    idx_j = find(~isnan(mask_withinRegion(idx_i,:)));
    idx_j(mask_withinRegion(idx_i,idx_j)==label(idx_j)')=[];
    mask_withinRegion(idx_i,idx_j)=nan;
end
mask_withinRegion(~isnan(mask_withinRegion))=1;
mask_withinRegion(isnan(mask_withinRegion))=0;

for i_project = idx_project
    % configuration
    switch Project{i_project}
        case 'HCP-YA'
            Condition = 'REST1';
            PE = 'LR';
            Experiment = [Type '_' Condition '_' PE];
        case 'HCP-D'
            Condition = 'REST1';
            PE = 'AP';
            Experiment = [Type '_' Condition '_' PE];
        case 'HCP-A'
            Experiment = Type;
    end

    % set path
    data_path = ['results/' Project{i_project} '_' Experiment '_' Resolution '/'];
    save_path = [data_path 'results/'];
    if ~exist(save_path,'dir');mkdir(save_path);end
    
    % subjects
    selected_path = dir(fullfile(data_path,'nonlinear*'));
    selected_path = selected_path(~contains({selected_path.name}, 'mu_47'));
    data_path = [data_path selected_path.name '/'];

    selected_files = dir(fullfile(data_path, '*.mat'));
    n_subs = length(selected_files);
    subjectList = cell(n_subs,1);
    for i_sub = 1:n_subs
        tmp_sub = split(selected_files(i_sub).name,'_');
        subjectList{i_sub} = tmp_sub{1};
    end

    for i_sub = n_subs:-1:1
        disp(['subject ' num2str(i_sub)]);
        load([data_path subjectList{i_sub} '_nonlinear_fitting.mat'],'coupling_mat');
        coupling_mat = coupling_mat - diag(diag(coupling_mat));
        
        Nonzeros_within(i_sub,1) = nnz(coupling_mat(mask_withinRegion==1));
        Coe = abs(coefficient(mask_withinRegion==1));
        Coe_within(i_sub,1) = Coe(2);
        MeanAmp_within(i_sub,1) = mean(abs(coupling_mat(coupling_mat~=0 & mask_withinRegion==1)));

      
        Nonzeros_between(i_sub,1) = nnz(coupling_mat(mask_withinRegion==0));
        Coe = abs(coefficient(mask_withinRegion==0));
        Coe_between(i_sub,1) = Coe(2);
        MeanAmp_between(i_sub,1) = mean(abs(coupling_mat(coupling_mat~=0 & mask_withinRegion==0)));
        
        % generate region-level coupling mat
        coupling_mat0 = coupling_mat;
        coupling_mat0(coupling_mat0==0)=nan;
        coupling_mat_region1 = zeros(n_regions,'single');
        coupling_mat_region2 = coupling_mat_region1;

        for i_region = 1:n_regions
            for j_region = 1:n_regions
                coupling_mat_region1(i_region,j_region) = mean(coupling_mat(label==i_region,label==j_region),"all");
                coupling_mat_region2(i_region,j_region) = mean(coupling_mat0(label==i_region,label==j_region),"all",'omitnan');
            end
        end

        coupling_mat_region1_all(:,:,i_sub) = coupling_mat_region1;
        coupling_mat_region2_all(:,:,i_sub) = coupling_mat_region2;

    end

    % save([save_path 'coupling_mat_properties.mat'],'Nonzeros','Fnorm','MeanAmp');
    save([save_path 'coupling_mat_properties_within_and_between.mat'],...
        'Nonzeros_within','Coe_within','MeanAmp_within',...
        'Nonzeros_between','Coe_between','MeanAmp_between');

    save([save_path 'coupling_mat_region.mat'],'coupling_mat_region1_all','coupling_mat_region2_all');
end