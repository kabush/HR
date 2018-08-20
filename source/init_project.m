%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% ----------------------------------------
%% Seed random number generator
rng(1,'twister');

%% ----------------------------------------
%% Initialize project param structure
proj = struct;

%% ----------------------------------------
%% Link tools
proj.path.kablab = '/home/kabush/lib/kablab/';
addpath(genpath(proj.path.kablab));

proj.path.export_fig = '/home/kabush/lib/export_fig/'
addpath(genpath(proj.path.export_fig));

%% ----------------------------------------
%% Project Flag Definitions
proj.flag.clean_build = 1;

%% ----------------------------------------
%% Project Path Definitions

%% Raw data
proj.path.raw_data = '/raw/bush/';
proj.path.raw_physio = 'physio';
proj.path.raw_logs = 'logfiles';
proj.path.raw_tabs = 'tabs';

%% Workspace
proj.path.home = '/home/kabush/workspace/';
proj.path.name = 'HRV';
proj.path.code = [proj.path.home,'code/',proj.path.name,'/'];
proj.path.data = [proj.path.home,'data/',proj.path.name,'/'];
proj.path.log =[proj.path.code,'log/'];
proj.path.fig = [proj.path.code,'fig/'];

%% ----------------------------------------
%% TICKET: Refactor output paths by data
%% mri[raw, i.e. clean and gm], hrv,
%% ----------------------------------------

%% Subject Lists
proj.path.subj_list = [proj.path.code,'subj_lists/'];

%% Design path (this is a meta source file)
proj.path.design = [proj.path.code,'design/'];

%% HRV paths
proj.path.hrv_beta = [proj.path.data,'hrv_beta/']; %ex and in put in
                                                   %the same
                                                   %directory
proj.path.hrv_bpm = [proj.path.data,'hrv_bpm/'];  
proj.path.hrv_mvpa_thresh = [proj.path.data,'hrv_mvpa_thresh/'];  
proj.path.hrv_haufe_thresh = [proj.path.data,'hrv_haufe_thresh/']; 
proj.path.hrv_mvpa_all = [proj.path.data,'hrv_mvpa_all/'];  
proj.path.hrv_haufe_all = [proj.path.data,'hrv_haufe_all/'];  
proj.path.hrv_haufe_permute_all = [proj.path.data,'hrv_haufe_permute_all/'];  
proj.path.hrv_haufe_permute_thresh = [proj.path.data,'hrv_haufe_permute_thresh/'];  

%% fMRI path
proj.path.fmri_clean = [proj.path.data,'fmri_clean/'];
proj.path.fmri_ex_beta = [proj.path.data,'fmri_ex_beta/'];
proj.path.fmri_in_beta = [proj.path.data,'fmri_in_beta/'];
proj.path.gm_mask = [proj.path.data,'gm_mask/'];

%% Target path
proj.path.trg = [proj.path.data,'target_ex/'];
proj.path.trg_in = [proj.path.data,'target_in/'];

%% MVPA path
proj.path.mvpa_fmri_ex_gs_cls = [proj.path.data,'mvpa_fmri_ex_gs_cls/'];
proj.path.mvpa_fmri_ex_gm_cls = [proj.path.data,'mvpa_fmri_ex_gm_cls/'];

%% Results logging file
proj.path.logfile = [proj.path.log,'logfile.txt'];
eval(['! rm ',proj.path.logfile]); % clear at initialization

%% Task file nomenclature
proj.path.name_id1 = 'Identify_run_1';
proj.path.name_id2 = 'Identify_run_2';
proj.path.name_rest = 'Rest';

%% ----------------------------------------
%% Project Parameter Definitions

%% Data source
proj.param.studies = {'CTM','INCA'};

%% fMRI Processing param
proj.param.TR = 2.0;
proj.param.slices = 37;
proj.param.slice_pattern = 'seq+z';
proj.param.do_anat = 'yes';
proj.param.do_epi = 'yes';
proj.param.tasks = 'identify'; %rest modulate 
proj.param.scans = 'run1 run2';
proj.param.rest_scans = 'run1';

%% *** Annoying extra parameter (silently swear at Philips software
%% engineers) ***  This shift is due to manner in which the design was
%% orginally constructed to accomodate the real-time
%% processing pipeline.  Prior to the Philips R5 upgrade
%% we were dropping 4 inital TRs, so the design built this in.
%% After the R5 upgrade we were dropping zero TRs but the
%% first TR is processed strangely and so is skipped. To
%% adjust for this we shift the design earlier in time by 3*TRs (TR=2s).
%% Basic problem is that the design assumed an 18 transient period
%% at the start of the identification runs which changed to 12 s
%% following R5 upgrades (shift was introduced to keep original
%% design files intact (possibly bad decision in the long run)
proj.param.r5_shift = -6;

%% Supervised learning labels of stimuli
proj.param.ex_id = 1;
proj.param.in_id = 2;
proj.param.feel_id = 3;

%% Cognitive dynamics labels of stimuli
proj.param.cogdyn_in_id = 1;
proj.param.cogdyn_cue_id = 2;
proj.param.cogdyn_feel_id = 3;
proj.param.cogdyn_rest_id = 4;

%% Likert scores adjustment parameters
proj.param.dummy_score = -1;
proj.param.mid_score = 5.0; % used to binarize classes

%% values representing binarized valence/arousal classes
proj.param.pos_class = 1;
proj.param.neg_class = -1;

%% Design construction fidelity (20 hz) 
%% all designs are manufactured at hi-fidelity
%% before downsampling to match fMRI acquisition
%% rate to minimize noise caused by slow TR
proj.param.hirez = 20;

%% Length of the tasks (in units of TR)
proj.param.n_trs_id1 = 282;
proj.param.n_trs_id2 = 282;
proj.param.n_trs_rest = 225;
proj.param.n_trs_mod1 = 310;
proj.param.n_trs_mod2 = 310;

%% Length of stimulus (in seconds)
proj.param.stim_t = 2;

%% Biopac channels
proj.param.chan_hr = 1;
proj.param.chan_rsp = 2;
proj.param.chan_scr = 3;
proj.param.chan_emg_zygo = 4;
proj.param.chan_emg_corr = 5;

%% Biopac recording freqs
proj.param.hz_scr = 2000;
proj.param.hz_emg = 2000;
proj.param.hz_hr = 2000;

%% HR analysis parameters
proj.param.hrv.intrv = 0.5:0.5:4.0;
proj.param.hrv.n_resamp = 30;
proj.param.hrv.thresh_seq = 0.1:0.1:3.4;
proj.param.hrv.convert_bpm = 60;

%% MVPA parameters
proj.param.mvpa_kernel = 'linear';
proj.param.mvpa_n_resamp = 1; % should be >= 30

%% Haufe parameters
proj.param.permute_nperm = 480;
proj.param.haufe_chunk = 10;

%% Plotting parameters
proj.param.plot.axisLabelFontSize = 18;
proj.param.plot.circleSize = 10;
proj.param.plot.white = [1,1,1];
proj.param.plot.light_grey = [.8,.8,.8];
proj.param.plot.dark_grey = [.6,.6,.6];
proj.param.plot.axis_nudge = 0.1;

%% ----------------------------------------
%% Write out initialized project structure
save('proj.mat','proj');