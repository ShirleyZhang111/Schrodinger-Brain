function subjectList = generate_TC(inCfg)

Project = inCfg.Project;
Experiment = inCfg.Experiment;
Resolution = inCfg.Resolution;
Condition = inCfg.Condition;
data_path = inCfg.data_path;

if strcmp(Project,'HCP-D')
    rawData_path = '/public/mig_old_storage/home1/ISTBI_data/HCP/HCP_Development/';
    File = dir(fullfile(rawData_path,'HCD*'));
    subjectList = cell(size(File));
    for m = 1:length(File)
        subjectList{m}=File(m).name(1:end-6);
    end
elseif strcmp(Project,'HCP-YA')
    rawData_path = ['/public/home/cyqian/HCP/' Experiment '/'];
    File = dir(fullfile(rawData_path));
    subjectList = cell(size(File));
    for m = 1:length(File)
        subjectList{m}=File(m).name;
    end
    subjectList(1:2)=[];
elseif strcmp(Project,'HCP-A')
    rawData_path ='/public/mig_old_storage/home1/ISTBI_data/HCP/HCP_Aging_2.0/';
    File = dir(fullfile(rawData_path,'HCA*'));
    subjectList = cell(size(File));
    for m = 1:length(File)
        subjectList{m}=File(m).name(1:end-6);
    end
end

if strcmp(Resolution,'8004+sub4008')
    load('idx_cifti_downsample_12k.mat','idx_downsample');
elseif strcmp(Resolution,'8004+sub2666')
    load('idx_cifti_downsample_10k.mat','idx_downsample');
end

for sub = length(subjectList):-1:1
    disp(['sub = ' num2str(sub)]);
    if strcmp(Project,'HCP-D')
        cifti_file = [rawData_path subjectList{sub} '_V1_MR/MNINonLinear/Results/' Experiment '/' Experiment '_Atlas_MSMAll_hp0_clean.dtseries.nii'];
        
    elseif strcmp(Project,'HCP-YA')
        switch Condition
            case 'rfMRI'
                cifti_file = [rawData_path subjectList{sub} '/' Experiment '_Atlas_hp2000_clean.dtseries.nii'];
                if ~exist(cifti_file,'file')
                    cifti_file = [rawData_path subjectList{sub} '/' Experiment '_Atlas.dtseries.nii'];
                end
            case 'tfMRI'
                cifti_file = [rawData_path subjectList{sub} '/' Experiment '_Atlas_MSMAll.dtseries.nii'];
                if ~exist(cifti_file,'file')
                    cifti_file = [rawData_path subjectList{sub} '/' Experiment '_Atlas.dtseries.nii'];
                end
        end
        
    elseif strcmp(Project,'HCP-A')
        Session = {'REST1','REST2'};
        PE = {'AP','PA'};
        for se = 1:2
            for pe = 1:2
                tmp_experiment = [Condition '_' Session{se} '_' PE{pe}];
                % cifti_file = [rawData_path subjectList{sub} '/MNINonLinear/Results/' tmp_experiment '/' tmp_experiment '_Atlas_MSMAll_hp0_clean.dtseries.nii'];
                cifti_file = [rawData_path subjectList{sub} '_V1_MR/MNINonLinear/Results/' tmp_experiment '/' tmp_experiment '_Atlas_MSMAll.dtseries.nii'];
                if exist(cifti_file,'file')~=0
                    break
                end
            end
            if exist(cifti_file,'file')~=0
                break
            end
        end
        
    end
    
    if ~exist(cifti_file,'file')
        subjectList{sub}=[];
        continue
    end
    
    fMRI = ft_read_cifti(cifti_file);
    BrainImg = zscore(fMRI.dtseries(idx_downsample,:),[],2);
    TC = single(BrainImg);
    
    save([data_path 'TC_' subjectList{sub} '.mat'],'TC');
end