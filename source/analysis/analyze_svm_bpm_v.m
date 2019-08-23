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
v_score = zscore(v_score);

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj); 

%% ----------------------------------------
%% iterate over study subjects

valence = [];
bpm = [];
svm = [];
subjects = [];

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
        svm = [svm;prds];
        valence = [valence;v_score];
        subjects = [subjects;repmat(cnt,numel(v_score),1)];
        
    end

end



%% ----------------------------------------
%% Group GLMM fit

%TEST
measures = double(valence);
pred1  = double(svm-mean(svm));
pred2 = double(bpm-mean(bpm));
subjects = double(subjects);

% % Model SVM only
% tbl = table(measures,pred1,subjects,'VariableNames',{'trg', ...
%                     'pred1','subj'});
% mdl_fe = fitlme(tbl,['trg ~ 1 + pred1']);
% mdl_re= fitlme(tbl,['trg ~ 1 + pred1 + (1+pred1|subj)']);
% 
% % Model BPM only
% tbl = table(measures,pred2,subjects,'VariableNames',{'trg', ...
%                     'pred2','subj'});
% mdl_fe = fitlme(tbl,['trg ~ 1 + pred2']);
% mdl_re= fitlme(tbl,['trg ~ 1 + pred2 + (1+pred2|subj)']);
% 
% Model SVM & BPM combined
tbl = table(measures,pred1,pred2,subjects,'VariableNames',{'trg', ...
                    'pred1','pred2','subj'});
mdl_fe = fitlme(tbl,['trg ~ 1 + pred1 + pred2']);
mdl_re= fitlme(tbl,['trg ~ 1 + pred1 + pred2 + (1+pred1|subj) + (1+pred2|subj)']);

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
scatter(pred1,measures,10,'MarkerFaceColor', ...
        proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);
hold on;

%% overlay the group VR skill plot
[y,idx] = sort(pred1);
spred1 = y;
y_hat = FE.Estimate(1) + FE.Estimate(2)*spred1; 
plot(spred1,y_hat,'r-','LineWidth',3);

%% overlay the estimate
[y,idx] = sort(pred2);
spred2 = y;
y_hat = FE.Estimate(1) + FE.Estimate(3)*spred2; 
plot(spred2,y_hat,'b-','LineWidth',3);


xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

xlabel('Predicted Valence Scores');
ylabel('Valence Scores');

%% explot hi-resolution figure
export_fig 'EX_predicted_v_summary.png' -r300  
eval(['! mv ',proj.path.code,'EX_predicted_v_summary.png ',proj.path.fig]);


%% ----------------------------------------
figure(2)
set(gcf,'color','w');

%% plot all the datapoints
scatter(pred2*FE.Estimate(3),measures-(FE.Estimate(1)+FE.Estimate(2)*spred1),10, ...
        'MarkerFaceColor', proj.param.plot.white,'MarkerEdgeColor', ...
        proj.param.plot.light_grey);
hold on;

%% overlay the estimate
[y,idx] = sort(pred2);
spred2 = y;
y_hat = FE.Estimate(3)*spred2; 
plot(spred2,y_hat,'b-','LineWidth',3);

xlim([-.1,.1]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;

ylabel('Valence Scores');
ylabel('Main effect residuals');
