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

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger([' Global Permutation Testing of MVPA Hyperplanes '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.haufe.v_permute_all]);
    eval(['! rm -rf ',proj.path.haufe.v_permute_all]);
    disp(['Creating ',proj.path.haufe.v_permute_all]);
    eval(['! mkdir ',proj.path.haufe.v_permute_all]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Load labels;
label_id = load([proj.path.trg.ex,'stim_ids.txt']);
v_label = load([proj.path.trg.ex,'stim_v_labs.txt']);
v_score = load([proj.path.trg.ex,'stim_v_scores.txt']);
ex_id = find(label_id==proj.param.trg.ex_id);

%% ----------------------------------------
%% iterate over permuations
Nperm = proj.param.haufe.npermute;
Nloop = Nperm + 1; %% first loop is true model structure
Nchunk = proj.param.haufe.chunk;

%% storage for group Haufe 
grp_haufe_v = zeros(172800,Nloop);  %%*** TICKET hardcoded ***

%% permutation significance levels
alpha05 = 0.05;
alpha01 = 0.01;
alpha001 = 0.001;

for i = 1:Nloop

    tic

    %%storage for group haufe 
    all_haufe_v_wts = zeros(172800,numel(subjs));
    all_haufe_v_mask = zeros(172800,numel(subjs));

    qlty_vec = zeros(1,numel(subjs));

    %% ----------------------------------------
    %% iterate over study subjects
    for j = 1:numel(subjs)
        
        %% extract subject info
        subj_study = subjs{j}.study;
        name = subjs{j}.name;
        id = subjs{j}.id;
        
        %% debug
        disp([subj_study,':',name,':i=',num2str(i)]);
        
        try

            %% Load gray matter mask 
            gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
            mask = double(gm_nii.img);
            brain_size=size(mask);
            mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
            in_brain=find(mask==1);  
            
            %% Load beta-series
            base_nii = load_nii([proj.path.betas.fmri_ex_beta,subj_study,'_',name,'_lss.nii']);
            brain_size = size(base_nii.img);
            
            %% Vectorize the base image
            base_img = vec_img_2d_nii(base_nii);
            base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
            
            %% Concatenate the MASKED base image
            all_img = base_img(in_brain,:)';
            
            %% Extract the extrinsic betas
            ex_img = zscore(all_img(ex_id,:));

            %% Extract the extrinsic valence scores;
            ex_v_score = zscore(v_score(ex_id));

            %% Peform quality check of generated features
            qlty = check_gm_img_qlty(ex_img);

            if(qlty.ok)

                %% Grab all labels in proper order
                label_ids = 1:numel(ex_v_score);
                
                %% Only first iteration is structure (remaining
                %% loops are permutations, therefore randomize labels
                if(i>1)
                    label_ids = randsample(label_ids,numel(label_ids));
                end

                %% Fit classifier
                mdl = fitrsvm(ex_img,ex_v_score(label_ids),'KernelFunction',proj.param.mvpa.kernel);

                %% Construct Valence Haufe tranform
                wts = mdl.Beta;
                haufe_v_wts = zscore(fast_haufe(ex_img,wts,Nchunk));
                all_haufe_v_wts(in_brain,j) = haufe_v_wts;
                all_haufe_v_mask(in_brain,j) = 1;
                
                qlty_vec(j)=1;
                

            else
                logger(['  -Failed Quality Check'],proj.path.logfile);
            end
            
        catch
            logger(['  -Haufe Error'],proj.path.logfile);
        end

    end

    %% ----------------------------------------
    %% Extract Quality fits        
    qlty_ids = find(qlty_vec==1);
    qlty_n = numel(qlty_ids);
    
    qlty_haufe_v_wts = all_haufe_v_wts(:,qlty_ids);
    qlty_haufe_v_mask = all_haufe_v_mask(:,qlty_ids);
    
    %% ----------------------------------------
    %% Group Mean Haufe transforms (1 saved per permutation)
    
    %% Group valence Haufe
    ahf_sum = sum(qlty_haufe_v_mask,2);
    row_ids_v = find(ahf_sum>(qlty_n/2));

    % Distribution version
    grp_haufe_v_dist = zeros(size(grp_haufe_v,1),qlty_n);
    grp_haufe_v_dist(row_ids_v,:) = qlty_haufe_v_wts(row_ids_v,:);

    % Mean version
    grp_haufe_v(row_ids_v,i) = mean(qlty_haufe_v_wts(row_ids_v,:),2);

    % T-score version
    if(i==1)
        grp_haufe_v_tstat = 0*grp_haufe_v(:,1);
        for k=1:numel(row_ids_v)
            [h p ci stat] = ttest(qlty_haufe_v_wts(row_ids_v(k),:));
            grp_haufe_v_tstat(row_ids_v(k),1)=stat.tstat;
        end
    end

    %% ----------------------------------------
    %% Do permutation test given samples available

    
    if(i>1)

        save([proj.path.haufe.v_permute_all,'grp_haufe_v_n=',num2str(i-1),'_of_N=',num2str(Nperm),'.mat'],'grp_haufe_v');

        if(i>2)
            eval(['! rm ',proj.path.haufe.v_permute_all,'grp_haufe_v_n=',num2str(i-2),'_of_N=',num2str(Nperm),'.mat']);
        end
    
        %% ----------------------------------------
        %% Valence
        sig_ids_05_v = [];
        sig_ids_01_v = [];
        sig_ids_001_v = [];
        
        for k=1:numel(row_ids_v)
            
            % Count extrem random samples
            Next = 0;
            if(grp_haufe_v(row_ids_v(k),1)>0)
                Next = numel(find(grp_haufe_v(row_ids_v(k),2:i)>grp_haufe_v(row_ids_v(k),1)));
            else
                Next = numel(find(grp_haufe_v(row_ids_v(k),2:i)<grp_haufe_v(row_ids_v(k),1)));
            end
            
            % Do 2-sided tests
            if(Next<round((alpha05/2)*i))
                sig_ids_05_v = [sig_ids_05_v,row_ids_v(k)];
            end
            
            if(Next<round((alpha01/2)*i))
                sig_ids_01_v = [sig_ids_01_v,row_ids_v(k)];
            end
            
            if(Next<round((alpha001/2)*i))
                sig_ids_001_v = [sig_ids_001_v,row_ids_v(k)];
            end
            
        end

        % ----------------------------------------
        % Save out: mean encoding of group gray-matter voxels
        if(numel(row_ids_v)>0)
            % mu_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v_tstat(row_ids_v,1),gm_nii,row_ids_v);
            mu_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(row_ids_v,1),gm_nii,row_ids_v);
            save_nii(mu_v_haufe_nii,[proj.path.haufe.v_permute_all,'mu_haufe_v_N=',num2str(Nperm),'.nii']);
        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.05) group
        % gray-matter voxels
        if(numel(sig_ids_05_v)>0)
            % mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v_tstat(sig_ids_05_v,1),gm_nii,sig_ids_05_v);
            mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_05_v,1),gm_nii,sig_ids_05_v);
            save_nii(mu_perm_v_haufe_nii,[proj.path.haufe.v_permute_all,'mu_perm_haufe_v_N=',num2str(Nperm),'_05.nii']);
        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.01) group
        % gray-matter voxels
        if(numel(sig_ids_01_v)>0)
            % mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_01_v,1),gm_nii,sig_ids_01_v);
            mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_01_v,1),gm_nii,sig_ids_01_v);
            save_nii(mu_perm_v_haufe_nii,[proj.path.haufe.v_permute_all,'mu_perm_haufe_v_N=',num2str(Nperm),'_01.nii']);
        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.001) group gray-matter voxels
        if(numel(sig_ids_001_v)>0)
            % mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_001_v,1),gm_nii,sig_ids_001_v);
             mu_perm_v_haufe_nii = build_nii_from_gm_mask(grp_haufe_v(sig_ids_001_v,1),gm_nii,sig_ids_001_v);
             save_nii(mu_perm_v_haufe_nii,[proj.path.haufe.v_permute_all,'mu_perm_haufe_v_N=',num2str(Nperm),'_001.nii']);
        end

    else
        
        % ----------------------------------------
        % Save out: wts of encodign for power analysis
        disp('saving distribution***************')
        save([proj.path.haufe.v_permute_all,'grp_haufe_distr.mat'],'grp_haufe_v_dist');

    end

        
    %% log completion of permuation testing
    logger([subj_study,':',name,':i=',num2str(i)],proj.path.logfile);

    toc        
end
