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
% % HR data (convert to csv for use in Kubios software)
% preprocess_hr_csv;
% 
% %% ******************************
% %% ******************************
% %% THERE IS A MANUAL KUBIOS 
% %% PROCESSING STEP HERE
% %% ******************************
% %% ******************************
% 
% % convert kubios *.mat format to *.csv format for python code
% preprocess_hr_kubios_reformat; 
% 
% %% ------------------------------------------------------------
% %% STEP 4: Format Extrinsic Stimuli Design
% format_ex_3dlss; 
% 
% %% ------------------------------------------------------------
% %% STEP 5: Calculate Extrinsic (EX) Stimuli Beta-Series
% 
% % Physio betas (these are the targets)
% calc_hr_ex_beta; % bpm trajectories | neutral filtering
% calc_hr_ex_bpm;  % bpm targets for mvpa
% 
% %% fMRI betas (these are the features)
% calc_fmri_ex_beta;
% 
% %% ------------------------------------------------------------
% %% STEP 7: Conduct MVPA for Secondary Measures
% mvpa_fmri_ex_gm_rgr_hr_inter_thresh; % inter-subj Gray-matter
%                                      % regression.
mvpa_fmri_ex_gm_rgr_hr_inter_all;    % inter-subj Gray-matter
                                     % regression.


%% ------------------------------------------------------------
%% STEP 8: Compare MVPA hyperplanes (V vs HR)
analyze_ex_v_vs_hr_mvpa;

%% ------------------------------------------------------------ 
%% STEP 9: Secondary Analysis of EX physiology (Neuropsychologia Paper)
analyze_ex_hr_mvpa;

% %% ------------------------------------------------------------ 
% %% STEP 10: Hyperplane analysis
% haufe_ex_gm_hr_mvpa_all_permute;

toc
