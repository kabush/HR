%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% ------------------------------------------------------------
%% Clean up matlab environment
matlab_reset; % debug

tic

%% ------------------------------------------------------------
%% Link all source code
addpath(genpath('./source/'));

%% ------------------------------------------------------------
%% STEP 1: Initialize the projects directories and parameters.
init_project;

% %% ------------------------------------------------------------
% %% STEP 2: Clear and reconstruct the project data folder
% clean_project;
% 
% %% ------------------------------------------------------------
% %% STEP 3: Preprocess raw data (wrangling, filtering, formatting)
% 
% % fMRI data
% preprocess_fmri;
% preprocess_mask;
% 
% % HRV data (convert to csv for use in Kubios software)
% preprocess_hrv_csv;
% 
% %% ******************************
% %% ******************************
% %% THERE IS A MANUAL KUBIOS 
% %% PROCESSING STEP HERE
% %% ******************************
% %% ******************************

% convert kubios *.mat format to *.csv format for python code
preprocess_hrv_kubios_reformat; 

%% ------------------------------------------------------------
%% STEP 4: Format Extrinsic Stimuli Design
format_ex_3dlss; 

%% ------------------------------------------------------------
%% STEP 5: Calculate Extrinsic (EX) Stimuli Beta-Series

% Physio betas (these are the targets)
calc_hrv_ex_beta; % bpm trajectories | neutral filtering
calc_hrv_ex_bpm;  % bpm targets for mvpa

%% fMRI betas (these are the features)
calc_fmri_ex_beta;

%% ------------------------------------------------------------
%% STEP 7: Conduct MVPA for Secondary Measures
mvpa_fmri_ex_gm_rgr_hrv_inter_thresh; % inter-subj Gray-matter
                                      % regression.
mvpa_fmri_ex_gm_rgr_hrv_inter_all;    % inter-subj Gray-matter
                                      % regression.

%% ------------------------------------------------------------ 
%% STEP 12: Secondary Analysis of EX physiology (Neuropsychologia Paper)
analyze_ex_hrv_mvpa;
% 
% %% ------------------------------------------------------------ 
% %% STEP 13: Hyperplane analysis
% haufe_ex_gm_hrv_mvpa_all_permute;


toc
