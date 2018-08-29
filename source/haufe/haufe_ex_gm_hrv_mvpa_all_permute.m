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
logger(['***********************************************************'],proj.path.logfile);
logger([' Permutation tested Haufe transform of HRV BPM hyperplanes '],proj.path.logfile);
logger(['***********************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.haufe.hrv_permute_all]);
    eval(['! rm -rf ',proj.path.haufe.hrv_permute_all]);
    disp(['Creating ',proj.path.haufe.hrv_permute_all]);
    eval(['! mkdir ',proj.path.haufe.hrv_permute_all]);
end

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_scores = load([proj.path.trg.ex,'stim_v_scores.txt']);

%%  extract only extrinsic stimuli
ex_ids = find(label_id==proj.param.trg.ex_id);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% load group HRV data
load([proj.path.physio.hrv_bpm,'all_bpm.mat']);
 
%% ----------------------------------------
%% Haufe parameters
Nboot = proj.param.haufe.npermute;
Nchunk = proj.param.haufe.chunk;

%% Storage for MVPA inputs
all_ex_img = [];
all_hrv_bpm = [];
all_subj_i = [];
all_qlty_i = [];

brain_size = 0;
%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    disp([subj_study,':',name]);

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
        all_hrv_bpm = [all_hrv_bpm;all_bpm];
        all_subj_i = [all_subj_i;subj_id];
        all_qlty_i = [all_qlty_i;i];
        
    end

end

grp_haufe_hrv = [];

for j = 1:(Nboot+1)

    %% randomly permute the order of ex_ids the same way for
    %% all subjects in inter-subject prediction
    seed_rnd_ids = randsample(1:numel(ex_ids),numel(ex_ids))';
    seed_ids = 1:numel(ex_ids)';
    perm_ids = [];

    cnt = 1;
    for i =1:numel(all_qlty_i)

        %% First iteration is correct; others are random
        if j>1
            this_ids = numel(ex_ids)*(cnt-1)+seed_rnd_ids;            
        else
            this_ids = numel(ex_ids)*(cnt-1)+seed_ids;
        end
        perm_ids = [perm_ids;this_ids];
        cnt = cnt+1;

    end

    %%storage for group haufe 
    all_haufe_wts = []; 
    all_haufe_mask = [];

    parfor i = 1:numel(all_qlty_i)

        tic % start timer

        qlty_i = all_qlty_i(i);

        %% extract subject info
        subj_study = subjs{qlty_i}.study;
        name = subjs{qlty_i}.name;
        disp(['n=',num2str(j),':',subj_study,':',name]);

        %% predict HRV (restricted data)
        [out,trg,~,mdl] = regress_inter_loocv(all_ex_img(perm_ids,:), ...
                                              all_hrv_bpm(perm_ids,:), ...
                                              all_subj_i,qlty_i, ...
                                              proj.param.mvpa.kernel);

        %% HAUFE-TRANSFORM
        wts = zeros(1,size(all_ex_img,2));
        wts(mdl.ids) = mdl.betas;
        haufe_wts = zscore(fast_haufe(all_ex_img,wts',Nchunk));

        %% STORE
        tmp_wts = zeros(prod(brain_size(1:3)),1);
        tmp_mask = zeros(prod(brain_size(1:3)),1);

        tmp_wts(in_brain) = haufe_wts;
        tmp_mask(in_brain) = 1;

        all_haufe_wts = [all_haufe_wts,tmp_wts];
        all_haufe_mask = [all_haufe_mask,tmp_mask];

        toc 

    end

    %% Group hrv Haufe
    ahf_sum = sum(all_haufe_mask,2);
    row_ids = find(ahf_sum>numel(all_qlty_i)/2);

    %% ========================================
    %% Save out grp haufe info (AS WE GO)
    grp_haufe_hrv = [grp_haufe_hrv,mean(all_haufe_wts,2)];
    save([proj.path.haufe.hrv_permute_all,'grp_haufe_hrv_n', ...
          num2str(Nboot),'_j',num2str(j),'.mat'],'grp_haufe_hrv');

    %% delete the old file that is no longer relevant
    if(j>1)
        eval(['! rm ',proj.path.haufe.hrv_permute_all,'grp_haufe_hrv_n', ...
              num2str(Nboot),'_j',num2str(j-1),'.mat']);
    end

end


%% ========================================
%% Find Bootstrap Significance Voxels
alpha05 = 0.05;
alpha01 = 0.01;
alpha001 = 0.001;

%% ----------------------------------------
%% Valence
sig_ids_05 = [];
sig_ids_01 = [];
sig_ids_001 = [];

for j=1:numel(row_ids);

    % ----------------------------------------
    % Count extrem random samples
    Next = 0;
    if(grp_haufe_hrv(row_ids(j),1)>0)
        Next = numel(find(grp_haufe_hrv(row_ids(j),2:end)>grp_haufe_hrv(row_ids(j),1)));
    else
        Next = numel(find(grp_haufe_hrv(row_ids(j),2:end)<grp_haufe_hrv(row_ids(j),1)));
    end

    % ----------------------------------------
    % Do 2-sided tests
    if(Next<round((alpha05/2)*(Nboot)))
        sig_ids_05 = [sig_ids_05,row_ids(j)];
    end

    if(Next<round((alpha01/2)*(Nboot)))
        sig_ids_01 = [sig_ids_01,row_ids(j)];
    end

    if(Next<round((alpha001/2)*(Nboot-1)))
        sig_ids_001 = [sig_ids_001,row_ids(j)];
    end

end

% Save out: mean encoding of group gray-matter voxels
mu_hrv_haufe_nii = build_nii_from_gm_mask(grp_haufe_hrv(row_ids,1),gm_nii,row_ids);
save_nii(mu_hrv_haufe_nii,[proj.path.haufe.hrv_permute_all,'mu_hrv_haufe_n',num2str(Nboot),'.nii']);

% Save out: mean encoding of bootstrap sign. (p<0.05) group gray-matter voxels
mu_boot_hrv_haufe_nii = build_nii_from_gm_mask(grp_haufe_hrv(sig_ids_05,1),gm_nii,sig_ids_05);
save_nii(mu_boot_hrv_haufe_nii,[proj.path.haufe.hrv_permute_all,'mu_boot_hrv_haufe_n',num2str(Nboot),'_05.nii']);

% Save out: mean encoding of bootstrap sign. (p<0.01) group gray-matter voxels
mu_boot_hrv_haufe_nii = build_nii_from_gm_mask(grp_haufe_hrv(sig_ids_01,1),gm_nii,sig_ids_01);
save_nii(mu_boot_hrv_haufe_nii,[proj.path.haufe.hrv_permute_all,'mu_boot_hrv_haufe_n',num2str(Nboot),'_01.nii']);
% Save out: mean encoding of bootstrap sign. (p<0.001) group gray-matter voxels
mu_boot_hrv_haufe_nii = build_nii_from_gm_mask(grp_haufe_hrv(sig_ids_001,1),gm_nii,sig_ids_001);
save_nii(mu_boot_hrv_haufe_nii,[proj.path.haufe.hrv_permute_all,'mu_boot_hrv_haufe_n',num2str(Nboot),'_001.nii']);
