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

%% Compute group mean trajectories
mu_pos_traj = mean(grp_trajs(pos_ids,:));
mu_neg_traj = mean(grp_trajs(neg_ids,:));

neg_ci_hi = [];
neg_ci_lo = [];
for i=1:size(mu_neg_traj,2)
    [h p ci stat] = ttest(grp_trajs(neg_ids,i));
    neg_ci_hi = [neg_ci_hi,ci(2)];
    neg_ci_lo = [neg_ci_lo,ci(1)];
end

pos_ci_hi = [];
pos_ci_lo = [];
for i=1:size(mu_pos_traj,2)
    [h p ci stat] = ttest(grp_trajs(pos_ids,i));
    pos_ci_hi = [pos_ci_hi,ci(2)];
    pos_ci_lo = [pos_ci_lo,ci(1)];
end


figure(1)
set(gcf,'color','w');

xseq = 0.5:0.5:4;
yseq = -.6:0.1:.6;

% Plot figure references
plot(xseq,0*xseq,'k-','LineWidth',2);
hold on;
plot(repmat(2,1,length(yseq)),yseq,'k:','LineWidth',2);

% Plot trajectory
plot(xseq,mu_pos_traj,'r-','LineWidth',3);
% plot(xseq,pos_ci_hi,'r--','LineWidth',1);
% plot(xseq,pos_ci_lo,'r--','LineWidth',1);

hold on;
plot(xseq,mu_neg_traj,'b-','LineWidth',3);
plot(xseq,neg_ci_hi,'b--','LineWidth',1);
plot(xseq,neg_ci_lo,'b--','LineWidth',1);
hold off;

xlim([0.5,4.0]);

% Manipulate labels
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Time since Stimulus onset (s)');
ylabel('HR change from pre-stimulus (bpm)');

% Export the image
export_fig 'EX_bpm_mean_trajectory.png' -r300
eval(['! mv ',proj.path.code,'EX_bpm_mean_trajectory.png ',proj.path.fig]);
