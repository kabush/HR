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
logger([' Compare Encoding (This Paper to SciReport 2018) '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

try

    %% Load encoding
    Nperm = proj.param.haufe.npermute;
    this_nii = load_nii([proj.path.haufe.v_permute_all,'mu_perm_haufe_v_N=',num2str(Nperm),'_05.nii']);
    this_img = vec_img_2d_nii(this_nii);
    in_this_img = find(abs(this_img)>0);

    %% *** TICKET *** Hardcoded path to prior atlas
    sci_nii = load_nii('/home/kabush/atlas/Bush_SciRep_2018/INCA_mu_boot_v_haufe_n1200_05.nii');
    sci_img = vec_img_2d_nii(sci_nii);    
    in_sci_img = find(abs(sci_img)>0);

    %% Find joint space
    in_img = intersect(in_this_img,in_sci_img);

    %% Group GLMM fit
    
    %Variables
    m_sci = double(sci_img(in_img));
    m_this = double(this_img(in_img));
    
    % Primary Test: 
    tbl = table(m_sci,m_this,'VariableNames',{'trg', ...
                        'pred1'});
    mdl_fe = fitlme(tbl,['trg ~ 1 + pred1']);
    mdl = mdl_fe;
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
    
    figure(1)
    set(gcf,'color','w');
    
    %% Plot effects
    scatter(m_this,m_sci,10,'MarkerFaceColor', ...
            proj.param.plot.white,'MarkerEdgeColor', ...
            proj.param.plot.light_grey);
    hold on;
    plot(sort(m_this),sort(m_this)*FE.Estimate(2)+FE.Estimate(1),'r-','LineWidth',3);

    hold off;
    fig = gcf;
    ax = fig.CurrentAxes;
    ax.FontSize = proj.param.plot.axisLabelFontSize;
    
    xlabel('Encoding Hyperplane');
    ylabel('Encoding Hyperplane (Bush et al., 2018)');
    
    %% explot hi-resolution figure
    export_fig 'EX_encoding_v_compare.png' -r300  
    eval(['! mv ',proj.path.code,'EX_encoding_v_compare.png ',proj.path.fig]);


    
    
catch
    logger(['  -Encoding Load Error'],proj.path.logfile);
end
