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
%% valenced stimuli from the HRV analysis
%% in line with 
%%
%%


%% Load in path data
load('proj.mat');

%% ----------------------------------------
%% Set-up Directory Structure for HRV
if(proj.flag.clean_build)
    logger(['Removing ',proj.path.hrv_bpm],proj.path.logfile);
    eval(['! rm -rf ',proj.path.hrv_bpm]);
    logger(['Creating ',proj.path.hrv_bpm],proj.path.logfile);
    eval(['! mkdir ',proj.path.hrv_bpm]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Load labels;
v_label = load([proj.path.trg,'stim_v_labs.txt']);
a_label = load([proj.path.trg,'stim_a_labs.txt']);
label_id = load([proj.path.trg,'stim_ids.txt']);
v_score = load([proj.path.trg,'stim_v_scores.txt']);
a_score = load([proj.path.trg,'stim_a_scores.txt']);

v_score = v_score(find(label_id==proj.param.ex_id));
a_score = a_score(find(label_id==proj.param.ex_id));

n_stim = numel(v_score);
seq_all = 1:n_stim;
seq_id1 = 1:(n_stim/2);
seq_id2 = ((n_stim/2)+1):n_stim;

%% allocate storage
grp_trajs = zeros(n_stim,numel(proj.param.hrv.intrv));
grp_cnts = zeros(n_stim,1);
grp_intrvs = [];


all_trajs = [];
all_v_scores = [];
all_i = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;
    try
        load([proj.path.hrv_beta,subj_study,'_',name,'_ex_betas.mat']);
    catch
        logger([subj_study,'_',name],proj.path.logfile);
        logger(['    Could not find hrv beta file for processing.'],proj.path.logfile);
    end

    trajs = [ex_betas.trajs1;ex_betas.trajs2];
    intrvs = [ex_betas.t_intrvs1]; 

    %%Change name to handle missing HRV
    if(~isempty(trajs))

        hrv_v_score = v_score;
        hrv_seq = seq_all;

        %%Adjust training data to handle mising HRV
        if(isempty(ex_betas.trajs1))
            trajs = ex_betas.trajs2;
            intrvs = ex_betas.t_intrvs2; 
            grp_intrvs = ex_betas.t_intrvs2; 
            hrv_seq = seq_id2;
            hrv_v_score = v_score(seq_id2);
        end
        
        if(isempty(ex_betas.trajs2))
            trajs = ex_betas.trajs1;
            intrvs = ex_betas.t_intrvs1; 
            grp_intrvs = ex_betas.t_intrvs2; 
            hrv_seq = seq_id1;
            hrv_v_score = v_score(seq_id1);
        end
        
        trajs = proj.param.hrv.convert_bpm*trajs; %convert to bpm
        grp_trajs(hrv_seq,:) = grp_trajs(hrv_seq,:)+trajs;
        grp_cnts(hrv_seq) = grp_cnts(hrv_seq)+1;
        
        %% Collect data for inter-subject LOOCV
        all_trajs = [all_trajs;trajs];
        all_v_scores = [all_v_scores;hrv_v_score];
        all_i = [all_i;repmat(i,size(trajs,1),1)];
        
    end

end

for i = 1:size(grp_trajs,1)
    grp_trajs(i,:) = grp_trajs(i,:)./grp_cnts(i);
end

% ----------------------------------------
% Creating and visualizing the measure
% ----------------------------------------
%
% 1) We will generate the group mean trajectory
%    for each stimulus image.  
% 
% 2) We will use a threshold about Likert=5 to
%    categorize images as pos/neg
%
% 3) We will increment the threshold from 0-3
%    and calculate the difference between
%    the mean category trajectory
%
% 4) We will plot the relationship between
%    distance between cateogry mean trajectories
%    showing that there exists a threshold at
%    which the difference starts to grow (this
%    demarks the edge of neutral. When the difference
%    of the differences is significant we can mark
%    this as the point at which categories are "PURE"
%    wrt the measure

Nsample = proj.param.hrv.n_resamp;;
not_found = 1;
thresh_seq = proj.param.hrv.thresh_seq;
thresh_vec = [];
all_b = [];
all_p = [];

for i=1:numel(thresh_seq)

    thresh = thresh_seq(i);
    thresh_vec = [thresh_vec, thresh];
    
    %% Identify pos./neg. classes
    mu = median(v_score);
    pos_ids = find((v_score-mu)>thresh);
    neg_ids = find((v_score-mu)<-thresh);

    v_pos = v_score(pos_ids);
    v_neg = v_score(neg_ids);
    
    grp_pos_trajs = grp_trajs(pos_ids,:);
    grp_neg_trajs = grp_trajs(neg_ids,:);


    all_mu_grp_rnd_pos_trajs = [];
    all_mu_grp_rnd_neg_trajs = [];

    for j=1:Nsample

        % randomly subsample
        pos_rnd_ids = randsample(1:numel(pos_ids),numel(pos_ids)-1);
        neg_rnd_ids = randsample(1:numel(neg_ids),numel(neg_ids)-1);
        
        grp_rnd_pos_trajs = grp_pos_trajs(pos_rnd_ids,:);
        grp_rnd_neg_trajs = grp_neg_trajs(neg_rnd_ids,:);
       
        all_mu_grp_rnd_pos_trajs = [all_mu_grp_rnd_pos_trajs;mean(grp_rnd_pos_trajs)];
        all_mu_grp_rnd_neg_trajs = [all_mu_grp_rnd_neg_trajs;mean(grp_rnd_neg_trajs)];
        
    end

    %% calcuate mean trajs over resampling
    mu_rnd_pos_trajs = mean(all_mu_grp_rnd_pos_trajs);
    mu_rnd_neg_trajs = mean(all_mu_grp_rnd_neg_trajs);

    %% find min of neg
    min_id = find(mu_rnd_neg_trajs==min(mu_rnd_neg_trajs));

    %% extract bpms at min
    pos_set = grp_pos_trajs(:,min_id);
    neg_set = grp_neg_trajs(:,min_id);

    xpred = [pos_set;neg_set];
    ypred = [v_pos;v_neg];

    [b stat] = robustfit(xpred,ypred);
    all_b = [all_b,b(2)];
    all_p = [all_p,stat.p(2)];
    
end

%% Find first significant prediction
%% This maximize kept values
sig_ids = find(all_p<0.05);
best_thresh = thresh_seq(sig_ids(1));

%% Compute POS/NEG sets
mu = median(v_score);
pos_ids = find((v_score-mu)>best_thresh);
neg_ids = find((v_score-mu)<-best_thresh);
grp_pos_trajs = grp_trajs(pos_ids,:);
grp_neg_trajs = grp_trajs(neg_ids,:);

%% Find minimum of negative stim deceleration
%% This is the new end point
mu_grp_neg_trajs = mean(grp_neg_trajs);
min_id = find(mu_grp_neg_trajs==min(mu_grp_neg_trajs));
fit_seq = 1:min_id;

%% Save out minimum values
pos_bpm = grp_pos_trajs(:,min_id);
neg_bpm = grp_neg_trajs(:,min_id);
all_bpm = grp_trajs(:,min_id);

save([proj.path.hrv_bpm,'best_thresh.mat'],'best_thresh');
save([proj.path.hrv_bpm,'all_bpm.mat'],'all_bpm');
save([proj.path.hrv_bpm,'pos_ids.mat'],'pos_ids');
save([proj.path.hrv_bpm,'pos_bpm.mat'],'pos_bpm');
save([proj.path.hrv_bpm,'neg_ids.mat'],'neg_ids');
save([proj.path.hrv_bpm,'neg_bpm.mat'],'neg_bpm');

%% ----------------------------------------
%% Plot intermediate results
%% ----------------------------------------

% ----------------------------------------
% This figure is likely for supplemental.
% It shows that the thresholding decreases
% the regression p-value 
% to convey numerous peices of methodology
% on a single plot.


figure(1);
set(gcf,'color','w');
scatter(thresh_vec,all_p,40,'MarkerFaceColor','b', ...
        'MarkerEdgeColor','b','MarkerFaceAlpha',0.2);

hold off;
xlabel('Threshold for Neutral Stimulus');
ylabel('Significance of HRV Model Prediction');
xlim([min(thresh_vec),max(thresh_vec)]);
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

% export hi-resolution figure
export_fig 'plot_hrv_thresh_significance.png' -r300  
eval(['! mv ',proj.path.code,'plot_hrv_thresh_significance.png ',proj.path.fig]);


% ----------------------------------------
% This is an extremely complex figure designed
% to convey numerous peices of methodology
% on a single plot.

figure(2)
set(gcf,'color','w');

% Plot zero acceleration line
x = linspace(proj.param.hrv.intrv(1),proj.param.hrv.intrv(end));
plot(x,0*x,'-k','LineWidth',2);
hold on;


% Plot grp mean trajectories for reference
mu_pos_trajs = mean(grp_pos_trajs);
mu_neg_trajs = mean(grp_neg_trajs);
std_pos_trajs = std(grp_pos_trajs);
std_neg_trajs = std(grp_neg_trajs);
ci_pos_trajs = std_pos_trajs/sqrt(size(grp_pos_trajs,1))*1.96;
ci_neg_trajs = std_neg_trajs/sqrt(size(grp_neg_trajs,1))*1.96;

plot(grp_intrvs,mu_pos_trajs,'-b','LineWidth',3);
plot(grp_intrvs,mu_neg_trajs,'-r','LineWidth',3);

plot1 = plot(grp_intrvs,mu_pos_trajs+ci_pos_trajs,'--b','LineWidth',2);
plot2 = plot(grp_intrvs,mu_pos_trajs-ci_pos_trajs,'--b','LineWidth',2);
plot3 = plot(grp_intrvs,mu_neg_trajs+ci_neg_trajs,'--r','LineWidth',2);
plot4 = plot(grp_intrvs,mu_neg_trajs-ci_neg_trajs,'--r','LineWidth',2);

line_alpha_value = 0.4;
plot1.Color(4) = line_alpha_value;
plot2.Color(4) = line_alpha_value;
plot3.Color(4) = line_alpha_value;
plot4.Color(4) = line_alpha_value;

% Plot timepoint when cue turns off
nudge = proj.param.plot.axis_nudge;;

y = linspace(min(mu_neg_trajs-ci_neg_trajs)-nudge,max(mu_pos_trajs+ci_pos_trajs)+nudge);
plot(0*y+2,y,':k','LineWidth',2);

% Figure out axis limits
cmb_trajs = [mu_pos_trajs+ci_pos_trajs,mu_neg_trajs-ci_neg_trajs];
xlim([min(grp_intrvs),max(grp_intrvs)]);
ylim([min(mu_neg_trajs-ci_neg_trajs)-nudge,max(mu_pos_trajs+ci_pos_trajs)+nudge]);

% Mark-up plot
xlabel('Times since stimulus onset (s)');
ylabel('Heartrate change from pre-stimulus (bpm)');
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'hrv_thresh_traj.png' -r300  
eval(['! mv ',proj.path.code,'hrv_thresh_traj.png ', ...
      proj.path.fig]);

%% ----------------------------------------
%% Perform Inter-subject LOOCV analysis

subj_i = unique(all_i);
all_bpm = all_trajs(:,min_id);


%% threshold out stimuli
mu = median(all_v_scores);
thresh_ids = find(abs(all_v_scores-mu)>best_thresh);

thresh_v_scores = all_v_scores(thresh_ids);
thresh_bpm = all_bpm(thresh_ids);
thresh_i = all_i(thresh_ids);

cv_rho_all = [];
cv_rho_thresh = [];
for subj=1:numel(subj_i);

    i = subj_i(subj);

    %% ----------------------------------------
    %% Calculate CV for all data
    tst_ids = find(all_i==i);
    trn_ids = setdiff(1:numel(all_v_scores),tst_ids);

    %% fit trn data
    [b stat] = robustfit(all_bpm(trn_ids),all_v_scores(trn_ids));

    %% predict tst data
    pred_v = all_bpm(tst_ids)*b(2)+b(1);

    %% Check correlation
    cv_rho_all = [cv_rho_all,corr(pred_v,all_v_scores(tst_ids))];

    %% ----------------------------------------
    %% Calculate CV for thresholded data
    tst_ids = find(thresh_i==i);
    trn_ids = setdiff(1:numel(thresh_v_scores),tst_ids);

    %% fit trn data
    [b stat] = robustfit(thresh_bpm(trn_ids),thresh_v_scores(trn_ids));

    %% predict tst data
    pred_v = thresh_bpm(tst_ids)*b(2)+b(1);

    %% Check correlation
    cv_rho_thresh = [cv_rho_thresh,corr(pred_v,thresh_v_scores(tst_ids))];


end

save([proj.path.hrv_bpm,'cv_rho_all.mat'],'cv_rho_all');
save([proj.path.hrv_bpm,'cv_rho_thresh.mat'],'cv_rho_thresh');

disp('Group Valence values for Pos/Neg');
disp(['Mu Pos V: ',num2str(mean(v_score(pos_ids)))]);
disp(['SD Pos V: ',num2str(std(v_score(pos_ids)))]);
disp(['Mu Neg V: ',num2str(mean(v_score(neg_ids)))]);
disp(['SD Pos V: ',num2str(std(v_score(neg_ids)))]);

disp('Group Arousal values for Pos/Neg');
disp(['Mu Pos A: ',num2str(mean(a_score(pos_ids)))]);
disp(['SD Pos A: ',num2str(std(a_score(pos_ids)))]);
disp(['Mu Neg A: ',num2str(mean(a_score(neg_ids)))]);
disp(['SD Pos A: ',num2str(std(a_score(neg_ids)))]);
