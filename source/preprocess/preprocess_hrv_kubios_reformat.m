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
        path_id_1 = [proj.path.physio.hrv_kubios_output,subj_study, ...
                     '_',name,'_Identify_run_1_kubios.mat'];
        load(path_id_1);
        rr = Res.HRV.Data.T_RR;
        csvwrite([proj.path.physio.hrv_kubios_reformat,subj_study, ...
                  '_',name,'_Identify_run_1_kubios_reformat.csv'],rr);
    catch
        disp('could not reformat Identify 1');
    end


    try
        path_id_2 = [proj.path.physio.hrv_kubios_output,subj_study, ...
                     '_',name,'_Identify_run_2_kubios.mat'];
        load(path_id_2);
        rr = Res.HRV.Data.T_RR;
        csvwrite([proj.path.physio.hrv_kubios_reformat,subj_study, ...
                  '_',name,'_Identify_run_2_kubios_reformat.csv'],rr);
    catch
        disp('could not reformat Identify 2');
    end

end
