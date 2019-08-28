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
logger(['Compare VAL predicting VAL prediction error     '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

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
        % Load SVM predictions
        load([proj.path.mvpa.v_all,subj_study,'_',name,'_prds.mat']);
        svm_exist = 1;
    catch
        logger(['   Missing svm'],proj.path.logfile);
    end

    if(svm_exist)
        
        cnt = cnt + 1;

        % Load everythin
        pval = [pval;zscore(prds.out)];
        val = [val;zscore(v_score)];
        sids = [sids;repmat(cnt,numel(v_score),1)];
        
    end

end


%% ----------------------------------------
%% Group GLMM fit

%Variables
m_err = double(zscore(val-pval));
m_val = double(val);
m_sids = double(sids);

% Primary Test: 
tbl = table(m_err,m_val,m_sids,'VariableNames',{'trg', ...
                    'pred1','subj'});
mdl_fe = fitlme(tbl,['trg ~ 1 + pred1']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred1 + (1+pred1|subj)']);

% Explore random effects across model types
fe_v_re = compare(mdl_fe,mdl_re);

mdl = mdl_fe;
if(fe_v_re.pValue<0.05);
    logger('   random effects matter',proj.path.logfile);
    mdl = mdl_re;
else
    logger('   random effects DO NOT matter',proj.path.logfile);
end

logger(' ',proj.path.logfile);

%% ----------------------------------------
%% Examine Main Effect

[~,~,FE] = fixedEffects(mdl);
if(FE.pValue(2)<0.05)
    logger('Fixed Effects are significant',proj.path.logfile);
    logger(['  FE(1)=',num2str(FE.Estimate(1))],proj.path.logfile);
    logger(['     p=',num2str(FE.pValue(1))],proj.path.logfile);
    logger(['  FE(2)=',num2str(FE.Estimate(2))],proj.path.logfile);
    logger(['     p=',num2str(FE.pValue(2))],proj.path.logfile);
else
    logger('Fixed Effects are **NOT** significant',proj.path.logfile);
    logger(['  FE(1)=',num2str(FE.Estimate(1))],proj.path.logfile);
    logger(['     p=',num2str(FE.pValue(1))],proj.path.logfile);
    logger(['  FE(2)=',num2str(FE.Estimate(2))],proj.path.logfile);
    logger(['     p=',num2str(FE.pValue(2))],proj.path.logfile);
end

%% ----------------------------------------
%% compute effect size

Rsqr = mdl.Rsquared.Ordinary;
Fsqr = Rsqr/(1-Rsqr);
logger(['  All Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['  Fsqr=',num2str(Fsqr)],proj.path.logfile);
logger(' ',proj.path.logfile);

%% ----------------------------------------
%% plot effects

figure(1)
set(gcf,'color','w');

%% format figure
ymin = -3;
ymax = 3;
xmin = -2;
xmax = 2;

%% plot all the datapoints
scatter(m_val,m_err,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);
hold on;

%% overlay the group effect
sx1 = linspace(min(m_val),max(m_val));
sy1 = FE.Estimate(1) + FE.Estimate(2)*sx1; 
plot(sx1,sy1,'r-','LineWidth',3);

xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Valence Scores');
ylabel('Valence Prediciton Error');

%% explot hi-resolution figure
export_fig 'EX_v_err_wrt_v_summary.png' -r300  
eval(['! mv ',proj.path.code,'EX_v_err_wrt_v_summary.png ',proj.path.fig]);

