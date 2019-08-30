%%========================================
%%========================================
%%
%% Keith Bush, PhD (2018)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

% %% Load in path data
% load('proj.mat');
% 
%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Compare Predicted VAL and BPM in prediting VAL  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Load Normative Valence scores;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
v_score = v_score(find(label_id==proj.param.trg.ex_id));

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj); 

%% Load trajectory param
load([proj.path.physio.hr_bpm,'min_traj_idx.mat']);


% search space 
thresh_seq = 0.0:0.2:3;
mu_likert = 5.0;

% effect size
pval_Rsqr = [];
bpm_Rsqr = [];
cmb_Rsqr = [];

% significance of fixed effect
pval_p = [];
bpm_p = [];
cmb_p1 = [];
cmb_p2 = [];

% dataset statistics
Nids = [];
mu_pos_set = [];
mu_neg_set = [];


for j = 1:numel(thresh_seq);

    thresh = thresh_seq(j);
    good_ids = find(abs(v_score-mu_likert)>thresh);
    
    %% ----------------------------------------
    %% iterate over study subjects
    val = [];
    bpm = [];
    pval = [];
    sids = [];
    cnt = 0;

    for i = 1:numel(subjs)
        
        %% extract subject info
        subj_study = subjs{i}.study;
        name = subjs{i}.name;
        
        %% debug
        logger([subj_study,':',name],proj.path.logfile);
        
        traj_exist = 0;
        svm_exist = 0;
        
        try
            % Load BPM data
            load([proj.path.physio.hr_bpm,subj_study,'_',name,'_trajs.mat']);
            traj_exist = 1;
        catch
            logger(['   Missing traj'],proj.path.logfile);
        end
        
        try
            % Load SVM predictions
            load([proj.path.mvpa.v_all,subj_study,'_',name,'_prds.mat']);
            svm_exist = 1;
        catch
            logger(['   Missing svm'],proj.path.logfile);
        end
        
        if(traj_exist & svm_exist)
            
            cnt = cnt + 1;
            
            % Load everythin
            bpm = [bpm;zscore(trajs(good_ids,min_traj_idx))];
            pval = [pval;zscore(prds.out(good_ids))];
            val = [val;zscore(v_score(good_ids))];
            sids = [sids;repmat(cnt,numel(good_ids),1)];
            
        end
        
    end

    %% ----------------------------------------
    %% Collect statistics
    Nids = [Nids,numel(good_ids)];


    mu_pos_set = [mu_pos_set,mean(v_score(good_ids(find(v_score(good_ids)>=mu_likert))))];
    mu_neg_set = [mu_neg_set,mean(v_score(good_ids(find(v_score(good_ids)<mu_likert))))];
    
    
    %% ----------------------------------------
    %% Group GLMM fit
    
    %Variables
    m_val = double(val);
    m_pval = double(pval);
    m_bpm = double(bpm);
    m_sids = double(sids);
    
    % MVPA Fixed-effect
    tbl = table(m_val,m_pval,m_sids,'VariableNames',{'trg','pred1','subj'});
    mdl = fitlme(tbl,['trg ~ 1 + pred1']);
    [~,~,FE] = fixedEffects(mdl);
    logger(['  FE(2)=',num2str(FE.Estimate(2))],proj.path.logfile);
    logger(['     p=',num2str(FE.pValue(2))],proj.path.logfile);
    logger(['  Rsqr=',num2str(mdl.Rsquared.Ordinary)],proj.path.logfile);
    pval_Rsqr = [pval_Rsqr,mdl.Rsquared.Ordinary];
    pval_p = [pval_p,FE.pValue(2)];

    logger(' ',proj.path.logfile)
    
    % BPM Fixed-effect
    tbl = table(m_val,m_bpm,m_sids,'VariableNames',{'trg','pred1','subj'});
    mdl = fitlme(tbl,['trg ~ 1 + pred1']);
    [~,~,FE] = fixedEffects(mdl);
    logger(['  FE(2)=',num2str(FE.Estimate(2))],proj.path.logfile);
    logger(['     p=',num2str(FE.pValue(2))],proj.path.logfile);
    logger(['  Rsqr=',num2str(mdl.Rsquared.Ordinary)],proj.path.logfile);
    bpm_Rsqr = [bpm_Rsqr,mdl.Rsquared.Ordinary];
    bpm_p = [bpm_p,FE.pValue(2)];

    logger(' ',proj.path.logfile)

    % Comb Fixed-effect
    tbl = table(m_val,m_pval,m_bpm,m_sids,'VariableNames',{'trg','pred1','pred2','subj'});
    mdl = fitlme(tbl,['trg ~ 1 + pred1 + pred2']);
    [~,~,FE] = fixedEffects(mdl);
    logger(['  FE(2)=',num2str(FE.Estimate(2))],proj.path.logfile);
    logger(['     p=',num2str(FE.pValue(2))],proj.path.logfile);
    logger(['  Rsqr=',num2str(mdl.Rsquared.Ordinary)],proj.path.logfile);
    cmb_Rsqr = [cmb_Rsqr,mdl.Rsquared.Ordinary];
    cmb_p1 = [cmb_p1,FE.pValue(2)];
    cmb_p2 = [cmb_p2,FE.pValue(3)];

end

%% ----------------------------------------
%% plot effects

%% ----------------------------------------
figure(2)
set(gcf,'color','w');

plot(thresh_seq,mu_pos_set,'b-','LineWidth',2);
hold on;
plot(thresh_seq,mu_neg_set,'r-','LineWidth',2);

%Bradley 2001
gids = find(mu_pos_set>7.1);
lids = find(mu_pos_set<7.3);
brad2001_ids = intersect(gids,lids);
scatter(mean(thresh_seq(brad2001_ids)),mean(mu_pos_set(brad2001_ids)),80,'k','filled');

gids = find(mu_neg_set>2.5);
lids = find(mu_neg_set<2.7);
brad2001_ids = intersect(gids,lids);
scatter(mean(thresh_seq(brad2001_ids)),mean(mu_neg_set(brad2001_ids)),80,'k','filled');

% Postive means
plot(thresh_seq,repmat(7.21,1,numel(thresh_seq)),'k:'); % Bradley
plot(thresh_seq,repmat(7.37,1,numel(thresh_seq)),'k:');
plot(thresh_seq,repmat(6.57,1,numel(thresh_seq)),'k:');
plot(thresh_seq,repmat(7.15,1,numel(thresh_seq)),'k:');
plot(thresh_seq,repmat(6.90,1,numel(thresh_seq)),'k:');

% Negative means
plot(thresh_seq,repmat(2.63,1,numel(thresh_seq)),'k:'); % Bradley
plot(thresh_seq,repmat(2.70,1,numel(thresh_seq)),'k:');
plot(thresh_seq,repmat(1.45,1,numel(thresh_seq)),'k:');
plot(thresh_seq,repmat(3.0,1,numel(thresh_seq)),'k:');
plot(thresh_seq,repmat(2.22,1,numel(thresh_seq)),'k:');

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Neutral Valence Threshold');
ylabel('Mean Valence Score');

export_fig 'EX_threshold_mu_v.png' -r300  
eval(['! mv ',proj.path.code,'EX_threshold_mu_v.png ',proj.path.fig]);

%% ----------------------------------------
figure(1)
set(gcf,'color','w');

plot(thresh_seq,cmb_Rsqr,'k:','LineWidth',2);
hold on;
plot(thresh_seq,pval_Rsqr,'b-','LineWidth',2);
plot(thresh_seq,bpm_Rsqr,'Color',[.85,.325,.098],'LineWidth',2);
plot(thresh_seq,bpm_Rsqr,'r-','LineWidth',2);

%% Denote non-significant BPM effects
sig_pts = 0*cmb_p2;
sig_pts(find(cmb_p2>0.05))=1;
non_sig_ids = find(sig_pts==1);
scatter(thresh_seq(non_sig_ids),bpm_Rsqr(non_sig_ids)+0.002,'k*');

%% Denote approximation of Bradley 2001
scatter(mean(thresh_seq(brad2001_ids)),mean(pval_Rsqr(brad2001_ids)),80,'k','filled');
scatter(mean(thresh_seq(brad2001_ids)),mean(bpm_Rsqr(brad2001_ids)),80,'k','filled');

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Neutral Valence Threshold');
ylabel('Prediction R-squared');

export_fig 'EX_threshold_Rsqr.png' -r300  
eval(['! mv ',proj.path.code,'EX_threshold_Rsqr.png ',proj.path.fig]);


%% ----------------------------------------
figure(3)
set(gcf,'color','w');

plot(thresh_seq,Nids/numel(v_score),'k-','LineWidth',2);
hold on;
scatter(mean(thresh_seq(brad2001_ids)),mean(Nids(brad2001_ids)/numel(v_score)),80,'k','filled');

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Neutral Valence Threshold');
ylabel('Fraction of Stimuli Surviving Threshold');

export_fig 'EX_threshold_fraction.png' -r300  
eval(['! mv ',proj.path.code,'EX_threshold_fraction.png ',proj.path.fig]);

