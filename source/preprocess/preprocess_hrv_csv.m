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

%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.physio.hrv_kubios]);
    eval(['! rm -rf ',proj.path.physio.hrv_kubios]);
    disp(['Creating ',proj.path.physio.hrv_kubios]);
    eval(['! mkdir ',proj.path.physio.hrv_kubios]);
end

%% Create the subjects to be analyzed (possible multiple studies)
subjs = load_subjs(proj);
disp(['Processing fMRI of ',num2str(numel(subjs)),' subjects']);

%% Preprocess fMRI of each subject in subjects list 
for i=1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;

    %% debug
    disp([subj_study,':',name]);

    try
        path_id_1 = [proj.path.raw_data,subj_study,'/physio/', ...
                     subj_study,'_',name,'/',subj_study,'_',name, ...
                     '_Identify_run_1.mat'];
        load(path_id_1);
        hrv = data(:,1);
        csvwrite([proj.path.physio.hrv_kubios,subj_study,'_',name,'_Identify_run_1.csv'],hrv);

    catch
        disp('could not load HRV files for Identify run 1');
    end


    try
        path_id_2 = [proj.path.raw_data,subj_study,'/physio/', ...
                     subj_study,'_',name,'/',subj_study,'_',name, ...
                     '_Identify_run_2.mat'];
        load(path_id_2);
        hrv = data(:,1);
        csvwrite([proj.path.physio.hrv_kubios,subj_study,'_',name,'_Identify_run_2.csv'],hrv);
    catch
        disp('could not load HRV files for Identify run 2');

    end






end
