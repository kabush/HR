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
rng(11056,'twister');

%% ----------------------------------------
%% Initialize project param structure
proj = struct;

%% ----------------------------------------
%% Link tools
proj.path.kablab = '/home/kabush/lib/kablab/';
addpath(genpath(proj.path.kablab));

proj.path.export_fig = '/home/kabush/lib/export_fig/'
addpath(genpath(proj.path.export_fig));

proj.path.nifti = '/home/kabush/lib/nifti/'
addpath(genpath(proj.path.nifti));


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
proj.path.name = 'HR';
proj.path.code = [proj.path.home,'code/',proj.path.name,'/'];
proj.path.data = [proj.path.home,'data/',proj.path.name,'/'];
proj.path.log =[proj.path.code,'log/'];
proj.path.fig = [proj.path.code,'fig/'];

%% ----------------------------------------
%% TICKET: Refactor output paths by data
%% mri[raw, i.e. clean and gm], hr,
%% ----------------------------------------

%% Subject Lists
proj.path.subj_list = [proj.path.code,'subj_lists/'];

%% Design path (this is a meta source file)
proj.path.design = [proj.path.code,'design/'];

%% Results logging file
formatOut = 'yyyy_mm_dd_HH:MM:SS';
t = datetime('now');
ds = datestr(t,formatOut);
proj.path.logfile = [proj.path.log,'logfile_',ds,'.txt'];

% eval(['! rm ',proj.path.logfile]); % clear at initialization

%% ----------------------------------------
%% Data Output Directory (All top-level names)
proj.path.mri.name = 'mri/';
proj.path.physio.name = 'physio/';
proj.path.betas.name = 'beta_series/';
proj.path.trg.name = 'target/';
proj.path.mvpa.name = 'mvpa/';
proj.path.haufe.name = 'haufe/';

%% ----------------------------------------
%% Specific Output Paths

%% fMRI paths (all)
proj.path.mri.mri_clean = [proj.path.data,proj.path.mri.name,'mri_clean/'];
proj.path.mri.gm_mask = [proj.path.data,proj.path.mri.name,'gm_mask/'];
proj.path.betas.fmri_ex_beta = [proj.path.data,proj.path.betas.name,'fmri_ex_beta/'];
proj.path.betas.fmri_in_beta = [proj.path.data,proj.path.betas.name,'fmri_in_beta/'];

%% HR paths (all)
proj.path.physio.hr_beta = [proj.path.data,proj.path.betas.name,'hr_beta/']; 
proj.path.physio.hr_bpm = [proj.path.data,proj.path.betas.name,'hr_bpm/'];  
proj.path.physio.hr_kubios = [proj.path.data,proj.path.physio.name,'hr_kubios/']; 
proj.path.physio.hr_kubios_output = [proj.path.raw_data,'kubios/']; ...

% master level repository of all kubios processed HR files
proj.path.physio.hr_kubios_reformat = [proj.path.data,proj.path.physio.name,'kubios_reformat/']; 

%% Target paths (all)
proj.path.trg.ex = [proj.path.data,proj.path.trg.name,'target_ex/'];
proj.path.trg.in = [proj.path.data,proj.path.trg.name,'target_in/'];

%% MVPA paths (all)
proj.path.mvpa.v_all = [proj.path.data,proj.path.mvpa.name,'v_all/'];  
proj.path.mvpa.bpm_all = [proj.path.data,proj.path.mvpa.name,'bpm_all/'];  

%% Haufe paths (all)
proj.path.haufe.v_permute_all = [proj.path.data,proj.path.haufe.name,'v_permute_all/'];  

%% Task file nomenclature
proj.path.task.name_id1 = 'Identify_run_1';
proj.path.task.name_id2 = 'Identify_run_2';
proj.path.task.name_rest = 'Rest';

%% ----------------------------------------
%% Project Parameter Definitions

%% Data source
proj.param.studies = {'CTM','INCA'};

%% fMRI Processing param
proj.param.mri.TR = 2.0;
proj.param.mri.slices = 37;
proj.param.mri.slice_pattern = 'seq+z';
proj.param.mri.do_anat = 'yes';
proj.param.mri.do_epi = 'yes';
proj.param.mri.tasks = 'identify'; %rest modulate 
proj.param.mri.scans = 'run1 run2';
proj.param.mri.rest_scans = 'run1';

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
proj.param.trg.r5_shift = -6;

%% Supervised learning labels of stimuli
proj.param.trg.ex_id = 1;
proj.param.trg.in_id = 2;
proj.param.trg.feel_id = 3;


%% values representing binarized valence/arousal classes
proj.param.trg.pos_class = 1;
proj.param.trg.neg_class = -1;

%% Likert scores adjustment parameters
proj.param.trg.dummy_score = -1;
proj.param.trg.mid_score = 5.0; % used to binarize classes

%% Cognitive dynamics labels of stimuli
proj.param.trg.cogdyn.in_id = 1;
proj.param.trg.cogdyn.cue_id = 2;
proj.param.trg.cogdyn.feel_id = 3;
proj.param.trg.cogdyn.rest_id = 4;

%% Start times of feel TRs relative to IN stimulus times
proj.param.trg.feel_times = 4:proj.param.mri.TR:10;
proj.param.trg.cue_times = 2; 
proj.param.trg.post_in_rest_times = 12;

%% Length of stimulus (in seconds)
proj.param.trg.stim_t = 2;

%% Length of the tasks (in units of TR)
proj.param.mri.n_trs_id1 = 282;
proj.param.mri.n_trs_id2 = 282;
proj.param.mri.n_trs_rest = 225;
proj.param.mri.n_trs_mod1 = 310;
proj.param.mri.n_trs_mod2 = 310;

%% Design construction fidelity (20 hz) 
%% all designs are manufactured at hi-fidelity
%% before downsampling to match fMRI acquisition
%% rate to minimize noise caused by slow TR
proj.param.betas.hirez = 20;

%% Biopac channels
proj.param.physio.chan_hr = 1;
proj.param.physio.chan_rsp = 2;
proj.param.physio.chan_scr = 3;
proj.param.physio.chan_emg_zygo = 4;
proj.param.physio.chan_emg_corr = 5;

%% Biopac recording freqs
proj.param.physio.hz_scr = 2000;
proj.param.physio.hz_emg = 2000;
proj.param.physio.hz_hr = 2000;

%% HR analysis parameters
proj.param.physio.hr.intrv = 0.5:0.5:4.0;
proj.param.physio.hr.n_resamp = 30;
proj.param.physio.hr.thresh_seq = 0.1:0.1:3.4;
proj.param.physio.hr.convert_bpm = 60;

%% MVPA parameters
proj.param.mvpa.kernel = 'linear';

%% Haufe parameters
proj.param.haufe.npermute = 1200;
proj.param.haufe.chunk = 10;

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