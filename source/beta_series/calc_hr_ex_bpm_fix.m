%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% This script systematically constructs
%% a threshold to filter out neutrally
%% valenced stimuli from the HR analysis
%% in line with 
%%
%% Load in path data
load('proj.mat');

%% ----------------------------------------
%% Set-up Directory Structure for HR
if(proj.flag.clean_build)
    logger(['Removing ',proj.path.physio.hr_bpm],proj.path.logfile);
    eval(['! rm -rf ',proj.path.physio.hr_bpm]);
    logger(['Creating ',proj.path.physio.hr_bpm],proj.path.logfile);
    eval(['! mkdir ',proj.path.physio.hr_bpm]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Load labels;
v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
v_score = v_score(find(label_id==proj.param.trg.ex_id));

%% allocate storage
grp_trajs = zeros(numel(v_score),numel(proj.param.physio.hr.intrv));
grp_cnts = zeros(numel(v_score),1);

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    try
        load([proj.path.physio.hr_beta,subj_study,'_',name,'_ex_betas.mat']);
    catch
        logger([subj_study,'_',name],proj.path.logfile);
        logger(['    Could not find hr beta file for processing.'],proj.path.logfile);
    end

    trajs = [ex_betas.trajs1;ex_betas.trajs2];
    intrvs = [ex_betas.t_intrvs1]; 

    %%Change name to handle missing HR
    if(~isempty(trajs))

        size(trajs)

        %% process trajs
        trajs = proj.param.physio.hr.convert_bpm*trajs; %convert to bpm
        grp_trajs = grp_trajs + trajs;
        grp_cnts = grp_cnts + 1;
        
        %% save out
        save([proj.path.physio.hr_bpm,subj_study,'_',name,'_trajs.mat'],'trajs');
        
    end

end

%% Compute stimulus wise trajectories (mean over group)
for i = 1:size(grp_trajs,1)
    grp_trajs(i,:) = grp_trajs(i,:)./grp_cnts(i);
end

%% Identify pos./neg. classes
mu = median(v_score);
pos_ids = find(v_score>mu);
neg_ids = find(v_score<mu);
mu_neg_traj = mean(grp_trajs(neg_ids,:));
min_traj_idx = find(mu_neg_traj==min(mu_neg_traj));

%% Save out index of trajectory minimum;
save([proj.path.physio.hr_bpm,'min_traj_idx.mat'],'min_traj_idx');