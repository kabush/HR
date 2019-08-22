%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Inter-subject LOOCV MVPA RGR GM Features -> HRs'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.hr_all]);
    eval(['! rm -rf ',proj.path.mvpa.hr_all]);
    disp(['Creating ',proj.path.mvpa.hr_all]);
    eval(['! mkdir ',proj.path.mvpa.hr_all]);
end

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);

%% extract only extrinsic stimuli
v_score = v_score(find(label_id==proj.param.trg.ex_id));

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% load group HR data
load([proj.path.physio.hr_bpm,'all_bpm.mat']);

%% Storage for MVPA inputs
all_ex_img = [];
all_v = [];
all_hr_bpm = [];
all_subj_i = [];
all_qlty_i = [];

%% Storage for MVPA output
out_bpm = [];
trg_bpm = [];
rho_bpm = [];
p_bpm = [];
out_v = [];
trg_v = [];
rho_v = [];
p_v = [];

%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% Load gray matter mask 
    gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
    mask = double(gm_nii.img);
    brain_size=size(mask);
    mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
    in_brain=find(mask==1);  

    %% Load beta-series
    base_nii = load_nii([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
    brain_size = size(base_nii.img);
    
    %% Vectorize the base image
    base_img = vec_img_2d_nii(base_nii);
    base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));

    %% Concatenate the MASKED base image
    all_img = base_img(in_brain,:)';
    
    %% Concatenate all label/subj identifiers
    subj_id = repmat(i,numel(all_bpm),1);
    
    %% Subselect extrinsic data
    ex_id = find(label_id==proj.param.trg.ex_id);
    ex_img = all_img(ex_id,:);

    %% Normalize within current subject
    ex_img = zscore(ex_img);

    %% Peform quality check of generated features
    qlty = check_gm_img_qlty(ex_img);

    if(qlty.ok)
        
        %% ----------------------------------------
        %% Build Inter-subjec structures
        all_ex_img = [all_ex_img;ex_img];
        all_v = [all_v;v_score];
        all_hr_bpm = [all_hr_bpm;all_bpm];
        all_subj_i = [all_subj_i;subj_id];
        all_qlty_i = [all_qlty_i;i];
        
    end

end

%% ----------------------------------------
%% Inter-subject modeling

for i = 1:numel(all_qlty_i)

    qlty_i = all_qlty_i(i);

    %% extract subject info
    subj_study = subjs{qlty_i}.study;
    name = subjs{qlty_i}.name;
    
    %% allocate result storage
    result = struct();

    %% predict HR (restricted data)
    [out,trg,~,mdl] = regress_inter_loocv(all_ex_img,all_hr_bpm, ...
                                          all_subj_i,qlty_i, ...
                                          proj.param.mvpa.kernel);
    
    result.bpm.out_bpm = out;
    result.bpm.trg_bpm = trg;
    result.bpm.p = mdl.stats.p;
    result.bpm.rho = mdl.stats.rho;
    result.bpm.beta = mdl.betas;
    result.bpm.ids = mdl.ids;

    %% predict Valence (restricted data)
    [out,trg,~,mdl] = regress_inter_loocv(all_ex_img,all_v,all_subj_i, ...
                                          qlty_i,proj.param.mvpa.kernel);

    %% assemble results
    result.v.out_bpm = out;
    result.v.trg_bpm = trg;
    result.v.p = mdl.stats.p;
    result.v.rho = mdl.stats.rho;
    result.v.beta = mdl.betas;
    result.v.ids = mdl.ids;

    %% save out prediction result
    save([proj.path.mvpa.hr_all,subj_study,'_',name,'_result.mat'],'result');

    %% log 
    logger([subj_study,'_',name,', bpm(rho)=',num2str(result.bpm.rho),[', ' ...
                        'v(rho)='],num2str(result.v.rho)],proj.path.logfile);

end
