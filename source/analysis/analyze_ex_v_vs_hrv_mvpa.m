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

% %% Set-up Directory Structure for fMRI betas
% if(proj.flag.clean_build)
%     disp(['Removing ',proj.path.haufe.cosine]);
%     eval(['! rm -rf ',proj.path.haufe.cosine]);
%     disp(['Creating ',proj.path.haufe.cosine]);
%     eval(['! mkdir ',proj.path.haufe.cosine]);
% end
% 
% % %% Initialize log section
% % logger(['************************************************'],proj.path.logfile);
% % logger([' Compare hyperplane angles (V vs HRV)           '],proj.path.logfile);
% % logger(['************************************************'],proj.path.logfile);
% 
% %% ----------------------------------------
% %% Extrinsic stimuli
% label_id = load([proj.path.trg.ex,'stim_ids.txt']);
% ex_ids = find(label_id==proj.param.trg.ex_id);
% 
% %% ----------------------------------------
% %% load group HRV data
% load([proj.path.physio.hrv_bpm,'all_bpm.mat']);
% 
% %% ----------------------------------------
% %% load subjs
% subjs = load_subjs(proj);
% 
% %% Storage for MVPA inputs
% all_ex_img = [];
% all_hrv_bpm = [];
% all_subj_i = [];
% all_qlty_i = [];
% 
% %% -------------------------------------------------
% %% Analyze state predictions using  thresholded data
% %% -------------------------------------------------
% brain_size = 0;
% for i = 1:numel(subjs)
% 
%     %% extract subject info
%     subj_study = subjs{i}.study;
%     name = subjs{i}.name;
%     id = subjs{i}.id;
%     disp([subj_study,':',name]);
% 
%     %% Load gray matter mask 
%     gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
%     mask = double(gm_nii.img);
%     brain_size=size(mask);
%     mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
%     in_brain=find(mask==1);  
% 
%     %% Load beta-series
%     base_nii = load_nii([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
%     brain_size = size(base_nii.img);
%     
%     %% Vectorize the base image
%     base_img = vec_img_2d_nii(base_nii);
%     base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
% 
%     %% Concatenate the MASKED base image
%     all_img = base_img(in_brain,:)';
%     
%     %% Concatenate all label/subj identifiers
%     subj_id = repmat(i,numel(all_bpm),1);
%     
%     %% Subselect extrinsic data
%     ex_id = find(label_id==proj.param.trg.ex_id);
%     ex_img = all_img(ex_id,:);
% 
%     %% Normalize within current subject
%     ex_img = zscore(ex_img);
% 
%     %% Peform quality check of generated features
%     qlty = check_gm_img_qlty(ex_img);
% 
%     if(qlty.ok)
%         
%         %% ----------------------------------------
%         %% Build Inter-subjec structures
%         all_ex_img = [all_ex_img;ex_img];
%         all_hrv_bpm = [all_hrv_bpm;all_bpm];
%         all_qlty_i = [all_qlty_i;i];
%         
%     end
% 
% end
% 
% v_haufe_wts = [];
% bpm_haufe_wts = [];
%   
% for j=1:numel(all_qlty_i)
% 
%     qlty_i = all_qlty_i(j);
%     
%     %% extract subject info
%     subj_study = subjs{qlty_i}.study;
%     name = subjs{qlty_i}.name;
%     disp([subj_study,':',name]);
%     
%     %% load valencea and hrv hyperplanes
%     load([proj.path.mvpa.hrv_all,subj_study,'_',name,'_result.mat']);
% 
%     %% HAUFE-TRANSFORM (valence)
%     wts = zeros(1,size(all_ex_img,2));
%     wts(result.v.ids) = result.v.beta;
%     haufe_wts = zscore(fast_haufe(all_ex_img,wts',proj.param.haufe.chunk));
%     v_haufe_wts = [v_haufe_wts,haufe_wts];
% 
%     %% HAUFE-TRANSFORM (bpm)
%     wts = zeros(1,size(all_ex_img,2));
%     wts(result.v.ids) = result.bpm.beta;
%     haufe_wts = zscore(fast_haufe(all_ex_img,wts',proj.param.haufe.chunk));
%     bpm_haufe_wts = [bpm_haufe_wts,haufe_wts];
%     
% end
% 
% N = numel(all_qlty_i);
% vv_cosTheta = zeros(N,N);
% bpmbpm_cosTheta = zeros(N,N);
% vbpm_cosTheta = zeros(N,N);
% 
% %% Compute Valence-Valence Cosine Similiarities
% for i=1:N
%     for j=1:N
%         vv_cosTheta(i,j) = dot(v_haufe_wts(:,i),v_haufe_wts(:,j))/(norm(v_haufe_wts(:,i))*norm(v_haufe_wts(:,j)));
%     end    
% end
% save([proj.path.haufe.cosine,'vv_cosTheta.mat'],'vv_cosTheta');
% 
% 
% %% Compute BPM-BPM Cosine Similiarities
% for i=1:N
%     for j=1:N
%         bpmbpm_cosTheta(i,j) = dot(bpm_haufe_wts(:,i),bpm_haufe_wts(:,j))/(norm(bpm_haufe_wts(:,i))*norm(bpm_haufe_wts(:,j)));
%     end    
% end
% save([proj.path.haufe.cosine,'bpmbpm_cosTheta.mat'],'bpmbpm_cosTheta');
% 
% %% Compute Valence-BPMe Cosine Similiarities
% for i=1:N
%     for j=1:N
%         vbpm_cosTheta(i,j) = dot(v_haufe_wts(:,i),bpm_haufe_wts(:,j))/(norm(v_haufe_wts(:,i))*norm(bpm_haufe_wts(:,j)));
%     end    
% end
% save([proj.path.haufe.cosine,'vbpm_cosTheta.mat'],'vbpm_cosTheta');


%% Gathering all the pre-computed similiarities
load([proj.path.haufe.cosine,'vv_cosTheta.mat']);
load([proj.path.haufe.cosine,'vbpm_cosTheta.mat']);
load([proj.path.haufe.cosine,'bpmbpm_cosTheta.mat']);

mu_vv = mean(vv_cosTheta,2);
mu_vbpm = mean(vbpm_cosTheta,2);
mu_bpmbpm = mean(bpmbpm_cosTheta,2);


cosim_all = [mu_vv',mu_bpmbpm',mu_vbpm'];
cosim_types = [repmat(0.5,1,numel(mu_vv)),repmat(2,1,numel(mu_bpmbpm)),repmat(3.5,1,numel(mu_vbpm))];

%% PLOTTING
figure(1)
set(gcf,'color','w');

h = notBoxPlot(cosim_all,cosim_types,'jitter',0.8); %,'MarkerSize',5);
xticklabels({'v vs v','\Delta vs \Delta','v vs \Delta'});

fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
ylabel('Encoding Cosine Similarity');
xlabel('Encoding Comparison');

export_fig 'hyperplane_boxplot.png' -r600
eval(['! mv ',proj.path.code,'hyperplane_boxplot.png ',proj.path.fig]);

%% Statistics output

%% V vs V
disp(['median V vs V: ',num2str(median(mu_vv))]);
disp(['p=',num2str(signrank(mu_vv))]);
disp('');

%% V vs HR change
disp(['median V vs delta HR: ',num2str(median(mu_vbpm))]);
disp(['p=',num2str(signrank(mu_vbpm))]);
disp('');

%% HR change vs HR change
disp(['median delta HR vs delta HR: ',num2str(median(mu_bpmbpm))]);
disp(['p=',num2str(signrank(mu_bpmbpm))]);
disp('');

%% V vs V (vs) HR vs HR
disp(['median V vs V (vs) HR vs HR, p=',num2str(ranksum(mu_vv,mu_bpmbpm))]);
disp('');

%% V vs V (vs) V vs HR
disp(['median V vs V (vs) V vs HR, p=',num2str(ranksum(mu_vv,mu_vbpm))]);
disp('');

%% V vs HR (vs) HR vs HR
disp(['median V vs HR (vs) HR vs HR, p=',num2str(ranksum(mu_vbpm,mu_bpmbpm))]);
disp('');



