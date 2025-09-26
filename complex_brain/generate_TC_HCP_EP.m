function subjectList = generate_TC_HCP_EP(inCfg)

rawData_path = inCfg.data_path;
save_path = inCfg.save_path;
label = inCfg.label;
is_largescale = inCfg.is_largescale;

File = dir(fullfile(rawData_path,'sub*'));
subjectList = cell(size(File));
for m = length(File):-1:1
    if File(m).isdir
        subjectList{m}=File(m).name;
    else
        subjectList(m) = [];
    end
end
if isfield(inCfg,'sub_index'); sub_index = inCfg.sub_index; else; sub_index = length(subjectList):-1:1; end
for sub = sub_index
    disp(['sub = ' num2str(sub)]);

    image_file = [rawData_path subjectList{sub} '/func/' subjectList{sub} '_task-rest_dir-AP_run-1_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz'];
    if ~exist(image_file,'file')
        subjectList{sub}=[];
        continue
    end
    
    fMRI = load_nii(image_file);
    fMRI = fMRI.img;
    fMRI = reshape(fMRI,[],size(fMRI,4));
    if ~is_largescale
        fMRI(label==0,:)=[];label(label==0)=[];

        % global signal regression
        GS = mean(fMRI,1,'omitnan');
        for node = 1:size(fMRI,1)
            beta = (fMRI(node,:)*GS')/(GS*GS');
            fMRI(node,:) = fMRI(node,:) - beta*GS;
        end

        for re = max(label,[],'all'):-1:1
            BrainImg(re,:) = mean(zscore(fMRI(label==re,:),[],2));
        end
    else
        mask = zeros(size(label));
        dsr = 3;
        mask(1:dsr:end,1:dsr:end,1:dsr:end) = 1;
        mask(label==0) = 0;
        fMRI(mask==0,:)=[];
        BrainImg = zscore(fMRI,[],2);
    end
    TC = single(BrainImg);clearvars BrainImg fMRI
    
    save([save_path 'TC_' subjectList{sub} '_GSR.mat'],'TC');
end