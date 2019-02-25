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
logger(['Prediction effect sizes of HRV and valence      '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% -------------------------------------------------
%% Analyze state predictions using  thresholded data
%% -------------------------------------------------

rho_bpm_thresh = [];
rho_v_thresh = [];
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    try
        load([proj.path.mvpa.hrv_thresh,subj_study,'_',name,'_result.mat']);
    catch
        disp('    Could not find HRV beta file for processing.');
    end

    rho_bpm_thresh = [rho_bpm_thresh,result.bpm.rho];
    rho_v_thresh = [rho_v_thresh,result.v.rho];

end

s_p_bpm_thresh = signrank(rho_bpm_thresh);
s_p_v_thresh = signrank(rho_v_thresh);

logger(['State effect hrv (thresh): ',...
      num2str(median(rho_bpm_thresh)),', p=',num2str(s_p_bpm_thresh)],proj.path.logfile);
logger(['State effect v (thresh): ',...
      num2str(median(rho_v_thresh)),', p=',num2str(s_p_v_thresh)],proj.path.logfile);

%% ----------------------------------------
%% Analyze state predictions using all data
%% ----------------------------------------

rho_bpm_all = [];
rho_v_all = [];
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    try
        load([proj.path.mvpa.hrv_all,subj_study,'_',name,'_result.mat']);
    catch
        disp('    Could not find HRV beta file for processing.');
    end

    rho_bpm_all = [rho_bpm_all,result.bpm.rho];
    rho_v_all = [rho_v_all,result.v.rho];

end

s_p_bpm_all = signrank(rho_bpm_all);
s_p_v_all = signrank(rho_v_all);

logger(['----------------------------------------'],proj.path.logfile);
logger(['State effect hrv (all): ',...
      num2str(median(rho_bpm_all)),', p=',num2str(s_p_bpm_all)],proj.path.logfile);
logger(['State effect v (all): ',...
      num2str(median(rho_v_all)),', p=',num2str(s_p_v_all)],proj.path.logfile);


%% -------------------------------------------------
%% Analyze HRV predictions
%% -------------------------------------------------
load([proj.path.physio.hrv_bpm,'cv_rho_all.mat']);
load([proj.path.physio.hrv_bpm,'cv_rho_thresh.mat']);

hrv_p_v_all = signrank(cv_rho_all);
hrv_p_v_thresh = signrank(cv_rho_thresh);

logger(['----------------------------------------'],proj.path.logfile);
logger(['HRV effect thresh: ', ...
        num2str(median(cv_rho_thresh)),', p=',num2str(hrv_p_v_thresh)],proj.path.logfile);
logger(['HRV effect all: ', ...
        num2str(median(cv_rho_all)),', p=',num2str(hrv_p_v_all)],proj.path.logfile);

logger(['----------------------------------------'],proj.path.logfile);
cmp_p_thresh = ranksum(cv_rho_thresh,rho_v_thresh);
logger(['state diff from hrv (thresh): p=',num2str(cmp_p_thresh)],proj.path.logfile);
cmp_p_all = ranksum(cv_rho_all,rho_v_all);
logger(['state diff from hrv (all): p=',num2str(cmp_p_all)],proj.path.logfile);

logger(['----------------------------------------'],proj.path.logfile);
cmp_p_hrv = ranksum(cv_rho_thresh,cv_rho_all);
logger(['hrv (all) from hrv (thresh): p=',num2str(cmp_p_hrv)],proj.path.logfile);
cmp_p_state = ranksum(rho_v_thresh,rho_v_all);
logger(['state (all) from state (thresh): p=',num2str(cmp_p_state)],proj.path.logfile);
cmp_p_state_bpm = ranksum(rho_bpm_thresh,rho_bpm_all);
logger(['state_bpm (all) from state_bpm (thresh): p=',num2str(cmp_p_state_bpm)],proj.path.logfile);
