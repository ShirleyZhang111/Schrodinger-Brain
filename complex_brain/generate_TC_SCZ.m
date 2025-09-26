function subjectList = generate_TC_SCZ(inCfg)

rawData_path = inCfg.data_path;
save_path = inCfg.save_path;
label = inCfg.label;

File = dir(fullfile(rawData_path,'sub*'));
subjectList = cell(size(File));
for m = length(File):-1:1
    if File(m).isdir
        subjectList{m}=File(m).name;
    else
        subjectList{m} = [];
    end
end

for sub = length(subjectList):-1:1
    disp(['sub = ' num2str(sub)]);

    image_file = [rawData_path subjectList{sub} '/func/' subjectList{sub} '_task-rest_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz'];
    
    if ~exist(image_file,'file')
        subjectList{sub}=[];
        continue
    end
    
    fMRI = load_nii(image_file);
    fMRI = fMRI.img;fMRI = reshape(fMRI,[],size(fMRI,4));
    for re = 1:max(label,[],'all')
        BrainImg(re,:) = mean(zscore(fMRI(label==re,:),[],2));
    end
    TC = single(BrainImg);
    
    save([save_path 'TC_' subjectList{sub} '.mat'],'TC');
end