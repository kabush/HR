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
logger(['Intra-subject LOOCV MVPA RGR GM Features -> BPM'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.bpm_all]);
    eval(['! rm -rf ',proj.path.mvpa.bpm_all]);
    disp(['Creating ',proj.path.mvpa.bpm_all]);
    eval(['! mkdir ',proj.path.mvpa.bpm_all]);
end

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);

%% ----------------------------------------
%% Load BPM index (calculated earlier)
load([proj.path.physio.hr_bpm,'min_traj_idx.mat']);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj); 

%% ----------------------------------------
%% iterate over study subjects
measures = [];
predictors = [];
subjects = [];

cnt = 0;

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    
    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    try
        
        %% Load BPM data
        load([proj.path.physio.hr_bpm,subj_study,'_',name,'_trajs.mat']);
        bpm = zscore(trajs(:,min_traj_idx));
        
        %% Load gray matter mask 
        gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
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
        
        %% Subselect extrinsic data
        ex_id = find(label_id==proj.param.trg.ex_id);
        ex_img = zscore(all_img(ex_id,:));
        
        %% Peform quality check of generated features
        qlty = check_gm_img_qlty(ex_img);
        
        if(qlty.ok)
            
            %% ----------------------------------------
            %% MVPA of Brain state onto BPM
            
            %% Fit model
            [out,trg,mdl,stats] = regress_intra_loocv(ex_img,bpm,proj.param.mvpa.kernel);
            disp(num2str(stats.rho));
            
            prds = struct();
            prds.out = out;
            prds.trg = trg;
            prds.rho = stats.rho;
            
            save([proj.path.mvpa.bpm_all,subj_study,'_',name,'_prds.mat'],'prds');
            
        else
            disp(['   quality check failed']);
            
        end
        
    catch
        disp(['   MVPA Error: possible missing beta series']);
    end
    
end
