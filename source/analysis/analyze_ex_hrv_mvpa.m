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

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Analyze thresholded data
%% ----------------------------------------

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

figure(1)
[b stat] = robustfit(rho_v_thresh,rho_bpm_thresh);
stat.p(2)
stat.p(1)
scatter(rho_v_thresh,rho_bpm_thresh);
hold on;
plot(sort(rho_v_thresh),sort(rho_v_thresh)*b(2)+b(1));

%% ----------------------------------------
%% Analyze all data
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

figure(1)
[b stat] = robustfit(rho_v_all,rho_bpm_all);
stat.p(2)
stat.p(1)
scatter(rho_v_all,rho_bpm_all);
hold on;
plot(sort(rho_v_all),sort(rho_v_all)*b(2)+b(1));

%% ----------------------------------------
load([proj.path.physio.hrv_bpm,'cv_rho_all.mat']);
load([proj.path.physio.hrv_bpm,'cv_rho_thresh.mat']);
