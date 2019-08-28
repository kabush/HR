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
% %% Initialize log section
% logger(['************************************************'],proj.path.logfile);
% logger(['Intra-subject LOOCV MVPA RGR GM Features -> Valence'],proj.path.logfile);
% logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Load Normative Valence scores;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
v_score = v_score(find(label_id==proj.param.trg.ex_id));

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj); 

%% ----------------------------------------
%% iterate over study subjects

val = [];
bpm = [];
pval = [];
sids = [];

load([proj.path.physio.hr_bpm,'min_traj_idx.mat']);

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
        disp(['   Missing traj']);
    end

    try
        % Load SVM predictions
        load([proj.path.mvpa.v_all,subj_study,'_',name,'_prds.mat']);
        svm_exist = 1;
    catch
        disp(['   Missing svm']);
    end

    if(traj_exist & svm_exist)
        
        cnt = cnt + 1;

        % Load everythin
        bpm = [bpm;zscore(trajs(:,min_traj_idx))];
        pval = [pval;zscore(prds.out)];
        val = [val;zscore(v_score)];
        sids = [sids;repmat(cnt,numel(v_score),1)];
        
    end

end



%% ----------------------------------------
%% Group GLMM fit

%Variables
m_pval = double(pval);
m_bpm = double(bpm);
m_sids = double(sids);

% Primary Test:  Model SVM & BPM combined
tbl = table(m_pval,m_bpm,m_sids,'VariableNames',{'trg','pred1','subj'});
mdl_fe = fitlme(tbl,['trg ~ 1 + pred1']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred1 + (1+pred1|subj)']);

%%Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_v_re.pValue<0.05);
    disp('   random effects matter');
    mdl = mdl_re;
else
    disp('   random effects DO NOT matter');
end

disp(' ');

%% ----------------------------------------
%% Examine Main Effect
[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    disp('Fixed Effects are significant');
    disp(['  p=',num2str(FE.pValue(2))]);
end

%% ----------------------------------------
%% compute effect size
Rsqr = mdl.Rsquared.Ordinary;
Fsqr = Rsqr/(1-Rsqr);
logger(['  All Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);

disp(' ');


%% ----------------------------------------
%% format figure
ymin = -3;
ymax = 3;
xmin = -3;
xmax = 3;


%% ----------------------------------------
figure(1)
set(gcf,'color','w');

%% plot all the datapoints
scatter(m_bpm,m_pval,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);
hold on;

%% overlay the SVM group effeect
[y,idx] = sort(m_bpm);
sx1 = y;
sy1 = FE.Estimate(1) + FE.Estimate(2)*sx1; 
plot(sx1,sy1,'r-','LineWidth',3);


xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Predicted Valence Scores');
ylabel('Valence Scores');

% %% explot hi-resolution figure
% export_fig 'EX_predicted_v_summary.png' -r300  
% eval(['! mv ',proj.path.code,'EX_predicted_v_summary.png ',proj.path.fig]);
