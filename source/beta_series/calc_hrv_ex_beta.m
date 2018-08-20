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
%% Load in path data
load('proj.mat');

%% ----------------------------------------
%% Set-up Directory Structure for HRV
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.hrv_beta]);
    eval(['! rm -rf ',proj.path.hrv_beta]);
    disp(['Creating ',proj.path.hrv_beta]);
    eval(['! mkdir ',proj.path.hrv_beta]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Split Designs into Run 1 and Run 2
id_path = [proj.path.trg,'stim_ids.txt'];
time_path = [proj.path.trg,'stim_times.1D'];

ids = load(id_path);
N = numel(ids)/2;
ids1 = ids(1:N);
ids2 = ids((N+1):end);

save([proj.path.code,'tmp/ids1.txt'],'ids1','-ascii');
save([proj.path.code,'tmp/ids2.txt'],'ids2','-ascii');

times = load(time_path);
N = numel(times)/2;
times1 =times(1:N);
times2 = times((N+1):end)-(proj.param.TR*proj.param.n_trs_id1);
save([proj.path.code,'tmp/times1.txt'],'times1','-ascii');
save([proj.path.code,'tmp/times2.txt'],'times2','-ascii');

%% ----------------------------------------
%% Compute HRVs over each subject individually
for i=1:numel(subjs)
   
    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    disp(['***********************************']);
    disp([subj_study,':',name]);

    %% Initialize hrv beta structure
    ex_betas = struct();

    %% ----------------------------------------
    %% Define run invariant paths
    out_path = [proj.path.code,'tmp/',subj_study,'_',name];
    rest_path = [proj.path.raw_data,subj_study,'/', ...
               proj.path.raw_physio,'/',subj_study,'_',name,'/', ...
               subj_study,'_',name,'_Rest.mat'];

    %% ----------------------------------------
    %% Process Identify Run 1

    %% Define input/outputs paths
    in_path = [proj.path.raw_data,subj_study,'/', ...
               proj.path.raw_physio,'/',subj_study,'_',name,'/', ...
               subj_study,'_',name,'_Identify_run_1.mat'];
    id_path = [proj.path.code,'tmp/ids1.txt'];
    time_path = [proj.path.code,'tmp/times1.txt'];

    disp('****MATLAB***');
    disp(out_path)

    %% RUN Kayla's HRV python code (will save out to) 
    eval(['! /usr/local/miniconda/bin/python ',proj.path.code,...
          'source/beta_series/hrv_analysis.py ',in_path,' ',id_path,' ',time_path,' ',out_path]);

    ex_betas.t_intrvs1 = [];
    ex_betas.trajs1 = [];
    try
        %% NEW FILES
        ex_betas.t_intrvs1 = load([proj.path.code,'tmp/',subj_study,'_',name,'_t_intrvs.txt']);
        ex_betas.trajs1 = load([proj.path.code,'tmp/',subj_study,'_',name,'_trajs.txt']);
    catch
        disp('Could not load HRV files for Identify run 1');
    end

    %% ----------------------------------------
    %% Process Identify Run 2

    %% Define input/outputs paths
    in_path = [proj.path.raw_data,subj_study,'/', ...
               proj.path.raw_physio,'/',subj_study,'_',name,'/', ...
               subj_study,'_',name,'_Identify_run_2.mat'];
    id_path = [proj.path.code,'tmp/ids2.txt'];
    time_path = [proj.path.code,'tmp/times2.txt'];

    out_path = [proj.path.code,'tmp/',subj_study,'_',name];
    rest_path = [proj.path.raw_data,subj_study,'/', ...
               proj.path.raw_physio,'/',subj_study,'_',name,'/', ...
               subj_study,'_',name,'_Rest.mat'];

    %% RUN Kayla's HRV python code (will save out to) 
    eval(['! /usr/local/miniconda/bin/python ',proj.path.code, ...
          'source/beta_series/hrv_analysis.py ',in_path,' ',id_path,' ',time_path,' ',out_path]);


    ex_betas.t_intrvs2 = [];
    ex_betas.trajs2 = [];
    try
        %% NEW FILES
        ex_betas.t_intrvs2 = load([proj.path.code,'tmp/',subj_study,'_',name,'_t_intrvs.txt']);
        ex_betas.trajs2 = load([proj.path.code,'tmp/',subj_study,'_',name,'_trajs.txt']);
    catch
        disp('Could not load HRV files for Identify run 1');
    end

    %% ----------------------------------------
    %% Save individual HRV structures
    save([proj.path.hrv_beta,subj_study,'_',name,'_ex_betas.mat'],'ex_betas');

end

%% ----------------------------------------
%% Clean-up
eval(['! rm ',proj.path.code,'tmp/*']);

