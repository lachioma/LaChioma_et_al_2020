%% DI distribution (Fig. 2B)
% You just need to load variable Disp:
% load('Vars\Disp_DGD_V1_LM_RL.mat');

Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
SFs   = {'0.01'};
nSFs = numel(SFs);

thr_std = 4;
nbins = 7;
CC = load_DGD_colors;

for ax = 1 : nAreas
    Area = Areas{ax};
    
    for i = 1 : nSFs
        SF = SFs{i};
        s  = find(strcmp({Disp.(Area).SF}, SF));
        ts = find(Disp.(Area)(s).thr_std==thr_std);

        nExps_tmp = length(Disp.(Area)(s).StatsDI_perExp{ts}.DI_eachExp);
        DI_distrib_eachExp = nan(nExps_tmp, nbins);
        for e = 1 : nExps_tmp
            DI_thisExp = Disp.(Area)(s).StatsDI_perExp{ts}.DI_eachExp{e};
            [DI_distrib_eachExp(e,:), edges] = histcounts(DI_thisExp, 0:1/nbins:1, 'Normalization','Probability');
        end
        DI_distrib_eachExp = DI_distrib_eachExp*100;
        
        DI_distrib_meanAcrossExps = mean(DI_distrib_eachExp, 1);
        DI_distrib_stdAcrossExps  =  std(DI_distrib_eachExp, 0,1);
        DI_distrib_semAcrossExps  = DI_distrib_stdAcrossExps/sqrt(nExps_tmp);


        hdi = figure;
        dx = (edges(2)-edges(1));
        xticks = edges(2:end)-dx/2;
        superbar( xticks, DI_distrib_meanAcrossExps,...
            'E', DI_distrib_semAcrossExps,...
            'BarWidth',0.9*dx, 'BarFaceColor', CC.(Areas{ax}){i}, 'BarEdgeColor','none',...
            'ErrorbarRelativeWidth',0.0, 'ErrorbarLineWidth',3);
        ylabel('Proportion of cells (%)' );
        xlim([-0.04 1.04])
        if strcmp(SF,'0.10')
            ylim([0 52])
        else
            ylim([0 52])
        end
        set(gca, 'YTick',[0:10:100])
        set(gca, 'XTick',[0 1],'XTickLabel',[0 1])
        xlabel('Disparity selectivity Index')
        title_str = ['Area ' Area ' - SF ' SF ' cpd'];
        title( title_str, 'Interpreter','none' ); 

    end
end


%% Nearest Neighbor analysis (appended exps, systematically across Areas and SFs)
Areas = {'RL'};%, 'LM', 'RL'};
SFs = {'0.01'};
thr_type = 1; thr_std  = 4;
thr_di   = 0.3;
plot_figs = [];
Edges = [0,10,20,40:40:200,300,400,Inf];
PhasesToExclude = [ ];
nAreas = numel(Areas);
nSFs = numel(SFs);
% NNexp=[]; NNall=[]; 
NNexp45=[]; NNall45=[]; NNexp90=[]; NNall90=[];

for i = 1 : nSFs
    SF = SFs{i};
    NNexp(i).SF = SF; NNall(i).SF = SF;
    for ax = 1 : nAreas
        Area = Areas{ax};
        if     strcmp(SF,'0.01')
            aDGDvarname = 'aDGD3';
        elseif strcmp(SF,'0.05')
            aDGDvarname = 'aDGD';
        elseif strcmp(SF,'0.10')
            aDGDvarname = 'aDGD2';
        end
        filename = ['J:\Alessandro La Chioma\from_I_drive\AnalyzedData\Alessandro\DGD\2017-01-06\'...
                        'vars_' Area '.mat'];
        aDGDap = load( filename, aDGDvarname );
        aDGDap = aDGDap.(aDGDvarname);
        
        [NNexp90(i).(Area)] = Map_IOPhase_NearestNeighbor_appended_exclude3( aDGDx, [90,270]       , thr_type, thr_std, thr_di, SF, Edges, PhasesToExclude, plot_figs);

    end
end





%%

ExpIDsall = [aDGDap.ROIs.ExpID];
SpatFreq = str2double(SF);
SFRoiNrs = find([aDGDap.ROIs.SpatFreq] == SpatFreq);
ExpIDs = unique(ExpIDsall(SFRoiNrs));
nExps  = numel(ExpIDs);


for e = 1:nExps
    
    NNexp(e).BinEdges = Edges;
    deltaIOPhaseNN = [];
    D = [];
    PrefIOPhase = [];
    PrefIOPhaseIx = [];
    DIOSFrespRoiNrs_excl = [];
    PrefIOPhase_excl = [];
%     deltaIOPhaseNN_mean = [];
%     deltaIOPhaseNN_eachExp = [];
%     Dist_um_eachExp        = [];
    
    expid = ExpIDs(e);
    Info = aDGDap.Info{expid};
    exptag = Info.Exp(1:5);




    save_dir  = ['J:\Alessandro La Chioma\from_I_drive\AnalyzedData\Alessandro\' Info.Mouse '\' Info.Spot '\' Info.Date '\'];
    FP = load([save_dir 'FieldParameters-' exptag '.mat']);
    load([save_dir 'ROI-' exptag '.mat']);
    aDGDap.FP{e} = FP;
    aDGDap.ROI{e} = ROI;
end

%% Plot RDC tuning curves (as in Fig. 7A)

% load('Fit__aRDS__V1_LM_RL__lsqcurvefit_asymGauss__TCmax__2020-06-30.mat')



%%% Select area, exp nr., ROIs to plot
    Areas = {'LM'}; ...{'V1','LM','RL'};
    eee = 7;  %%%%%%%%%%%%%%%%%%%%%%%%%%
    ROIsToPlot = 1 : 10; % indixes of selROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Responsiveness criteria:
thr_std = 4;
thr_type = 1;
thr_rep_fraction_rds = 0.5;
thr_Rsquared = 0.5;%0.5;

RDStype  = 'aRDS';
nAreas = numel(Areas);
CC = load_DGD_colors(1); cx = 1;

for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    nExps = length(Fit(ai).(RDStype));

    for e = eee%1 : 3%%%%%nExps
        
%         respRDS_1  = above_stdThreshold_DGD(aD_RDS(ai).(RDStype )(e) , thr_std, thr_type, [], thr_rep_fraction_rds);
%         respRDS_2  = above_stdThreshold_DGD(aD_RDS(ai).(RDStype2)(e) , thr_std, thr_type, [], thr_rep_fraction_rds);
%         respRDS    = intersect(respRDS_1, respRDS_2);
%         rsqRDS_1   = above_RsquaredThreshold(Fit(ai).(RDStag )(e).TC, thr_Rsquared, 'Rsquared_meanTC_1');
%         rsqRDS_2   = above_RsquaredThreshold(Fit(ai).(RDStag )(e).TC, thr_Rsquared, 'Rsquared_meanTC_2');
%         selROIs    = mintersect(respRDS, rsqRDS_1, rsqRDS_2);
        selROIs = [];
        for i = 1 : length(Fit(ai).(RDStype)(e).TC)
            if isfield(Fit(ai).(RDStype)(e).TC(i).GoF, 'Rsquared_meanTC') && Fit(ai).(RDStype)(e).TC(i).GoF.Rsquared_meanTC > thr_Rsquared
                selROIs = [selROIs; i]; %#ok<AGROW>
            end
        end
        
        disp([' Selected ROIs: ' num2str(length(selROIs)) ' / ' num2str(length(Fit(ai).(RDStype )(e).TC)) ]);
        
        for i = ROIsToPlot % 1 : 10%20%length(selROIs)
            
            nr = selROIs(i);

            hf=figure(100000000+ax*10000000+e*10000+nr); hold on;
                fig_pos = get(gcf, 'position'); fig_pos(3) = fig_pos(3)*1.3;
                set(gcf, 'position',fig_pos); % [left, bottom, width, height]

                x_val = Fit(ai).(RDStype)(e).StimSettings.bin_centers_deg;

            hsa=subaxis(1,4,1,1,3,1); hold on;
                supererr(x_val, Fit(ai).(RDStype)(e).TC(nr).TuningCurve_mean, [],Fit(ai).(RDStype)(e).TC(nr).TuningCurve_sem,...
                    '|','Color',CC.(Area){cx})
                plot( x_val, Fit(ai).(RDStype)(e).TC(nr).TuningCurve_mean, 'Color',CC.(Area){cx},...
                    'Marker','o','LineWidth',2, 'MarkerFaceColor',CC.(Area){cx}, 'MarkerEdgeColor','none');
%                 plot( Fit(ai).(RDStype)(e).TC(nr).FittedData_x,...
%                     Fit(ai).(RDStype)(e).TC(nr).FittedData, 'r', 'LineWidth',2);
                
                
                xlim([min(x_val) max(x_val)])
                set(gca, 'XTick',[-20:10:20])
                set(gca,'Color','none')
                box off

            subaxis(1,4,4,1,1,1)
            Info = Fit(ai).(RDStype)(e).Info;
                str_info= {...
                    [ Info.Mouse ];...
                    [ 'Area ' Info.Area ];...
                    [ Info.Spot ];...
                    [ 'Date ' Info.Date ,...
                      'exp' Info.Exp ];...' - ' 'IDexp ' num2str(Fit(ai).IDexp(e))];...
                    [ ' ' ];...
                    [ 'ROI # ' num2str(nr) ]; ...
                    [ ' ' ];
                    [ 'thr_std = ' num2str(thr_std) ];...
                    [ 'thr_type = ' num2str(thr_type) ];...
                    [ 'thr_rep_fraction_rds = ' num2str(thr_rep_fraction_rds)];...
                    [ 'thr_Rsquared = ' num2str(thr_Rsquared)];...
                    };

                FitParam = Fit(ai).(RDStype)(e).TC(nr).FitParam; str_param=cell(length(FitParam),1);
                for p = 1 : length(FitParam)
                    str_param{p} = [FitParam{p} ' = ' num2str(Fit(ai).(RDStype)(e).TC(nr).BestFitEsts(p)) ];
                end
%                     str_param{p+1} = ['phase shift' ' = ' num2str((Fit(ai).(RDStag)(e).TC(nr).BestFitEsts(6)-Fit(ai).(RDStag)(e).TC(nr).BestFitEsts(8))) ];

                str_param(end+1:end+4) = {...
                    [ ' ' ];...
                    [ 'SSE-meanTC_1 = '  num2str(Fit(ai).(RDStype)(e).TC(nr).GoF.SSE_meanTC) ];...
%                     [ 'SSE-meanTC_2 = '  num2str(Fit(ai).(RDStype)(e).TC(nr).GoF.SSE_meanTC_2) ];...
                    [ 'RMSE-meanTC_1 = ' num2str(Fit(ai).(RDStype)(e).TC(nr).GoF.RMSE_meanTC) ];...
%                     [ 'RMSE-meanTC_2 = ' num2str(Fit(ai).(RDStype)(e).TC(nr).GoF.RMSE_meanTC_2) ];...
                    [ 'R2-meanTC_1 = '   num2str(Fit(ai).(RDStype)(e).TC(nr).GoF.Rsquared_meanTC)];...
%                     [ 'R2-meanTC_2 = '   num2str(Fit(ai).(RDStype)(e).TC(nr).GoF.Rsquared_meanTC_2)];...
                    };
                text( 0, 0.5, [str_info; str_param],...
                'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','none');
                axis off

        end
    end
end


%% Plot RDC tuning width and asymmetry (Fig. 7C,D)
%
% load('Fit__aRDS__V1_LM_RL__lsqcurvefit_asymGauss__TCmax__2020-06-30.mat')

clearvars -except  db aD aD_RDS aD_RPS Fit Fit_RDS Fit_RPS

% Responsiveness criteria:
thr_std = 4;
thr_type = 1;
thr_Rsquared = 0.5;%0.5;

Areas = {'V1','LM','RL'};
% Areas = {'LM'};

RDStype1  = 'aRDS';
% RDStype2  = 'aRDS2';
RDStype = RDStype1; select_rois = 'intersect(respRDS_1, rsqRDS_1)';
% RDStype = RDStype2; select_rois = 'intersect(respRDS_2, rsqRDS_2)';

nAreas = numel(Areas);
nExps_max = 9;
Sigma_eachExp      = cell(nExps_max, nAreas, 2);
SigmaDelta_eachExp = cell(nExps_max, nAreas);
SigmaSum_eachExp   = cell(nExps_max, nAreas);
Sigma_allEachArea      = cell(nAreas, 2);
SigmaDelta_allEachArea = cell(1, nAreas);
SigmaSum_allEachArea   = cell(1, nAreas);
SigmaDelta_medEachExp  = cell(1, nAreas);
SigmaSum_medEachExp  = cell(1, nAreas);
selROIs_eachExp  = cell(nExps_max, nAreas);

for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    nExps = length(Fit(ai).(RDStype));

    for e = 1 : nExps
        
        respRDS_1  = above_RsquaredThreshold(Fit(ai).(RDStype1 )(e).TC, 0, 'Rsquared_meanTC');
%         respRDS_2  = above_RsquaredThreshold(Fit(ai).(RDStype2 )(e).TC, 0, 'Rsquared_meanTC');
        rsqRDS_1   = above_RsquaredThreshold(Fit(ai).(RDStype1 )(e).TC, thr_Rsquared, 'Rsquared_meanTC');
%         rsqRDS_2   = above_RsquaredThreshold(Fit(ai).(RDStype2 )(e).TC, thr_Rsquared, 'Rsquared_meanTC');
        selROIs    = eval(select_rois);
        
        disp([' Selected ROIs: ' num2str(length(selROIs)) ' / ' num2str(length(Fit(ai).(RDStype )(e).TC)) ]);
        
        Sigma      = nan(length(selROIs),2);
        SigmaDelta = nan(length(selROIs),1);
        SigmaSum   = nan(length(selROIs),1);
        
        for i = 1 : length(selROIs)
            
            nr = selROIs(i);
            
            Sigma(i,1)    = Fit(ai).(RDStype)(e).TC(nr).BestFitEsts(2);
            Sigma(i,2)    = Fit(ai).(RDStype)(e).TC(nr).BestFitEsts(5);
            SigmaDelta(i) = Sigma(i,2) - Sigma(i,1);
            SigmaSum(i)   = Sigma(i,2) + Sigma(i,1);
        end
        
        if length(selROIs) > 0
            Sigma_eachExp{e,ax,1}     = Sigma(:,1);
            Sigma_eachExp{e,ax,2}     = Sigma(:,2);
            SigmaDelta_eachExp{e,ax}  = SigmaDelta;
            SigmaSum_eachExp{e,ax}    = SigmaSum;
            selROIs_eachExp{e,ax}     = selROIs;
        end
        
    end
    Sigma_allEachArea{ax,1}    = cell2mat(Sigma_eachExp(:,ax,1));
    Sigma_allEachArea{ax,2}    = cell2mat(Sigma_eachExp(:,ax,2));
    SigmaDelta_allEachArea{ax} = cell2mat(SigmaDelta_eachExp(:,ax));
    SigmaSum_allEachArea{ax}   = cell2mat(SigmaSum_eachExp(:,ax));
    
    SigmaDelta_medEachExp{ax}  = cellfun(@(x) median(abs(x)), SigmaDelta_eachExp(:,ax),'UniformOutput',true);
    SigmaSum_medEachExp{ax}    = cellfun(@median, SigmaSum_eachExp(:,ax),'UniformOutput',true);
end

%%% Figures

CC = load_DGD_colors(1); cx = 1; 
CC2 = load_DGD_colors(2); cx = 1; 

% % Tuning asymmetry (Sigma2 - Sigma1) - cell-wise
% figure;
% for ax = 1 : nAreas
%     Area = Areas{ax};
%     h(ax) = subplot(nAreas,1,ax);
%     hold on;
%     histogram( (SigmaDelta_allEachArea{ax}),[-20:5:20],'Normalization','probability',...
%         'FaceColor',CC.(Area){cx});
% %     set(gca,'XTick',[-180:30:180]);
%     title(['mean = ' num2str(mean(SigmaDelta_allEachArea{ax})) ' - median = ' num2str(median(SigmaDelta_allEachArea{ax}))])
% end
% linkaxes(h);

% % Tuning asymmetry abs(Sigma2 - Sigma1) - cell-wise
% figure;
% for ax = 1 : nAreas
%     Area = Areas{ax};
%     h(ax) = subplot(nAreas,1,ax);
%     hold on;
%     histogram( abs(SigmaDelta_allEachArea{ax}),[-20:5:20],'Normalization','probability',...
%         'FaceColor',CC.(Area){cx});
% %     set(gca,'XTick',[-180:30:180]);
%     title(['mean = ' num2str(mean(abs(SigmaDelta_allEachArea{ax}))) ' - median = ' num2str(median(abs(SigmaDelta_allEachArea{ax})))])
% end
% linkaxes(h);

% Tuning asymmetry abs(Sigma2 - Sigma1) - exp-wise
figure;
sigma_delta_edges = [0:2.5:20, Inf];
n_bins = length(sigma_delta_edges) - 1;
for ax = 1 : nAreas
    Area = Areas{ax};
    h(ax) = subplot(nAreas,1,ax);
    hold on
    SigmaDelta_counts_eachExp = nan(length(SigmaSum_eachExp(:,ax)), n_bins);
    for e = 1 : length(SigmaSum_eachExp(:,ax))
        SigmaDelta_counts_eachExp(e,:) = histcounts(abs(SigmaDelta_eachExp{e,ax}),sigma_delta_edges,'Normalization','probability');
    end
    SigmaDelta_counts_meanEachExp = nanmean(SigmaDelta_counts_eachExp, 1);
    SigmaDelta_counts_semEachExp = nanstd(SigmaDelta_counts_eachExp,0,1) ./ sqrt(sum(~isnan(SigmaDelta_counts_eachExp),1));
    hb=superbar(1:n_bins, SigmaDelta_counts_meanEachExp, 'E',SigmaDelta_counts_semEachExp,...
        'BarFaceColor',CC.(Area){cx}, 'ErrorbarColor',CC2.(Area){cx}, 'ErrorbarStyle','|');
    set(hb, 'ShowBaseLine', 'off')
    if ax==nAreas
        set(gca, 'XTick',[0.5:n_bins-0.5,n_bins], 'XTickLabel',[num2cellstr(sigma_delta_edges(1:end-1)),{['>' num2str(sigma_delta_edges(end-1))]} ] )
        xlabel('Tuning asymmetry |\sigma2 - \sigma1| (deg)')
    else
        set(gca, 'XTick',[])
    end
end
linkaxes(h); ylim([-0.02 0.50]);

% % Tuning width cell-wise (Sigma1+Sigma2, as in La Chioma et al. 2019)
% figure;
% for ax = 1 : nAreas
%     Area = Areas{ax};
%     h(ax) = subplot(nAreas,1,ax);
%     hold on;
%     histogram( SigmaSum_allEachArea{ax},[0:3:42],'Normalization','probability',...
%         'FaceColor',CC.(Area){cx});
% %     set(gca,'XTick',[-180:30:180]);
%     title(['mean = ' num2str(mean(SigmaSum_allEachArea{ax})) ' - median = ' num2str(median(SigmaSum_allEachArea{ax}))])
% end
% linkaxes(h);

% Tuning width exp-wise (Sigma1+Sigma2, as in La Chioma et al. 2019)
figure;
sigma_sum_edges = [0:3:27, Inf];
n_bins = length(sigma_sum_edges) - 1;

for ax = 1 : nAreas
    Area = Areas{ax};
    h(ax) = subplot(nAreas,1,ax);
    hold on
    SigmaSum_counts_eachExp = nan(length(SigmaSum_eachExp(:,ax)), n_bins);
    for e = 1 : length(SigmaSum_eachExp(:,ax))
        SigmaSum_counts_eachExp(e,:) = histcounts(SigmaSum_eachExp{e,ax},sigma_sum_edges,'Normalization','probability');
    end
    SigmaSum_counts_meanEachExp = nanmean(SigmaSum_counts_eachExp, 1);
    SigmaSum_counts_semEachExp = nanstd(SigmaSum_counts_eachExp,0,1) ./ sqrt(sum(~isnan(SigmaSum_counts_eachExp),1));
    hb=superbar(1:n_bins, SigmaSum_counts_meanEachExp, 'E',SigmaSum_counts_semEachExp,...
        'BarFaceColor',CC.(Area){cx}, 'ErrorbarColor',CC2.(Area){cx}, 'ErrorbarStyle','|');
    set(hb, 'ShowBaseLine', 'off')
    if ax==nAreas
        set(gca, 'XTick',[0.5:n_bins-0.5,n_bins], 'XTickLabel',[num2cellstr(sigma_sum_edges(1:end-1)),{['>' num2str(sigma_sum_edges(end-1))]} ] )
        xlabel('Tuning width (\sigma1 + \sigma2) (deg)')
    else
        set(gca, 'XTick',[])
    end
    
end
linkaxes(h); ylim([-0.01 0.3]);


%%% Statistics

SigmaSum_counts_meanOfMediansAcrossExps = nan(nAreas, 1);
SigmaSum_counts_semOfMediansAcrossExps = nan(nAreas, 1);
SigmaDelta_meanOfMediansAcrossExps = nan(nAreas, 1);
SigmaDelta_semOfMediansAcrossExps = nan(nAreas, 1);
for ax = 1 : nAreas
    SigmaSum_counts_meanOfMediansAcrossExps(ax) = nanmean(SigmaSum_medEachExp{ax});
    SigmaSum_counts_semOfMediansAcrossExps(ax) = nanstd(SigmaSum_medEachExp{ax}) / sqrt(sum(~isnan(SigmaSum_medEachExp{ax})));
    
    SigmaDelta_meanOfMediansAcrossExps(ax) = nanmean(SigmaDelta_medEachExp{ax});
    SigmaDelta_semOfMediansAcrossExps(ax) = nanstd(SigmaDelta_medEachExp{ax}) / sqrt(sum(~isnan(SigmaDelta_medEachExp{ax})));
end

p_val = nan(2,2);

TestType = 'KW';
CType = ['bonferroni'];
% Tuning asymmetry abs(Sigma2 - Sigma1)
fprintf('\n Tuning asymmetry abs(Sigma2 - Sigma1) (cell-wise)\n\n')
[p_val(1,1), p_pairwise, ~, ~] = perform_anova1(cellfun(@abs,SigmaDelta_allEachArea,'UniformOutput',false), Areas, TestType, CType);

TestType = 'KW';
CType = ['bonferroni'];
% Tuning width (Sigma1+Sigma2)
fprintf('\n Tuning width (Sigma1 + Sigma2) (cell-wise)\n\n')
[p_val(2,1), p_pairwise, ~, ~] = perform_anova1(SigmaSum_allEachArea, Areas, TestType, CType);

% Tuning asymmetry abs(Sigma2 - Sigma1)
TestType = 'KW';
CType = ['bonferroni'];
fprintf('\n Tuning asymmetry abs(Sigma2 - Sigma1) (exp-wise)\n\n')
[p_val(1,2), p_pairwise, ~, ~] = perform_anova1(SigmaDelta_medEachExp, Areas, TestType, CType);
% Multiple comparisons with Mann-Whitney test:
pw_p = nan(3,1);
[pw_p(1)] = ranksum(SigmaDelta_medEachExp{1},SigmaDelta_medEachExp{2});
[pw_p(2)] = ranksum(SigmaDelta_medEachExp{1},SigmaDelta_medEachExp{3});
[pw_p(3)] = ranksum(SigmaDelta_medEachExp{2},SigmaDelta_medEachExp{3});
% bonferroni correction:
pw_p_bonf = pw_p * 3;

% Tuning width (Sigma1+Sigma2)
TestType = 'KW';
CType = ['bonferroni'];
fprintf('\n Tuning width (Sigma1 + Sigma2) (exp-wise)\n\n')
[p_val(2,2), p_pairwise, ~, ~] = perform_anova1(SigmaSum_medEachExp, Areas, TestType, CType);
% Multiple comparisons with Mann-Whitney test:
pw_p = nan(3,1);
[pw_p(1)] = ranksum(SigmaSum_medEachExp{1},SigmaSum_medEachExp{2});
[pw_p(2)] = ranksum(SigmaSum_medEachExp{1},SigmaSum_medEachExp{3});
[pw_p(3)] = ranksum(SigmaSum_medEachExp{2},SigmaSum_medEachExp{3});
% bonferroni correction:
pw_p_bonf = pw_p * 3;


%% Get number of tuned/untuned cells (to run for Figure 8C)

% load('Fit__aRDS__aRDS2__aRDS_aRDS2__V1_LM_RL__lsqcurvefit_gabor__TCmax__2019-10-16.mat')

Areas = {'V1','LM','RL'};
thr_Rsquared = 0.5;%0.5;
RDStype1  = 'aRDS';
RDStype2 = 'aRDS2';
RDStype12   = 'aRDS_aRDS2';

% clear all variables containing 'resp':
namesWorkspace = who; outStr = regexpi(namesWorkspace, 'resp'); ind = ~cellfun('isempty',outStr); varlist = namesWorkspace(ind);
if ~isempty(varlist), clear(varlist{:}); end

N_resp = get_number_tuned_cells(Fit, Areas, RDStype1, RDStype2, RDStype12, ...
                    thr_Rsquared);

%% Fraction tuned/untuned cells with RDC (Fig. 8C)
% Before this section, run the section:
% Get number of tuned/untuned cells (to run for Fig. 8C)

% load('Fit__aRDS__aRDS2__aRDS_aRDS2__V1_LM_RL__lsqcurvefit_gabor__TCmax__2019-10-16.mat')

nAreas = numel(Areas);
n_respRDS_any_untuned_meanAcrossExps = nan(nAreas, 1);
n_respRDS_any_tuned_1_only_meanAcrossExps = nan(nAreas, 1);
n_respRDS_any_tuned_2_only_meanAcrossExps = nan(nAreas, 1);
n_respRDS_12_tuned_12_meanAcrossExps = nan(nAreas, 1);
n_respRDS_any_untuned_semAcrossExps = nan(nAreas, 1);
n_respRDS_any_tuned_1_only_semAcrossExps = nan(nAreas, 1);
n_respRDS_any_tuned_2_only_semAcrossExps = nan(nAreas, 1);
n_respRDS_12_tuned_12_semAcrossExps = nan(nAreas, 1);


for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    n_respRDS_any_untuned_meanAcrossExps(ax) = mean( N_resp.n_respRDS_any_untuned_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} );
    n_respRDS_any_tuned_1_only_meanAcrossExps(ax) = mean( N_resp.n_respRDS_any_tuned_1_only_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} );
    n_respRDS_any_tuned_2_only_meanAcrossExps(ax) = mean( N_resp.n_respRDS_any_tuned_2_only_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} );
    n_respRDS_12_tuned_12_meanAcrossExps(ax) = mean( N_resp.n_respRDS_12_tuned_12_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} );
    
    n_respRDS_any_untuned_semAcrossExps(ax) = std( N_resp.n_respRDS_any_untuned_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} ) / sqrt(length(N_resp.n_respRDS_any_eachExp{ax}));
    n_respRDS_any_tuned_1_only_semAcrossExps(ax) = std( N_resp.n_respRDS_any_tuned_1_only_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} ) / sqrt(length(N_resp.n_respRDS_any_eachExp{ax}));
    n_respRDS_any_tuned_2_only_semAcrossExps(ax) = std( N_resp.n_respRDS_any_tuned_2_only_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} ) / sqrt(length(N_resp.n_respRDS_any_eachExp{ax}));
    n_respRDS_12_tuned_12_semAcrossExps(ax) = std( N_resp.n_respRDS_12_tuned_12_eachExp{ax} ...
        ./ N_resp.n_respRDS_any_eachExp{ax} ) / sqrt(length(N_resp.n_respRDS_any_eachExp{ax}));
    
end

CC = load_DGD_colors(1);
label_str = {'untuned: ';'only RDS-C: ';'only RDS-A: ';'both: '};

figure;
Y  = [ n_respRDS_any_untuned_meanAcrossExps';
      n_respRDS_any_tuned_1_only_meanAcrossExps';
      n_respRDS_any_tuned_2_only_meanAcrossExps';
      n_respRDS_12_tuned_12_meanAcrossExps' ] * 100;
EY = [ n_respRDS_any_untuned_semAcrossExps';
      n_respRDS_any_tuned_1_only_semAcrossExps';
      n_respRDS_any_tuned_2_only_semAcrossExps';
      n_respRDS_12_tuned_12_semAcrossExps' ] * 100;
C = CC.sf1';
superbar(Y, 'E', EY, 'BarFaceColor', C); 

set(gca, 'XTick',[1:size(Y,1)], 'XTickLabel',label_str, 'XTickLabelRotation',45)
set(gca, 'YTick',[0:10:100], 'YLim',[0 100])


%% Get phase offset and amplitude ratio for corr-RDS and anticorr RDS (to run for Figure 8E-G plots)

% load('Fit__aRDS__aRDS2__aRDS_aRDS2__V1_LM_RL__lsqcurvefit_gabor__TCmax__2019-10-16.mat')

clearvars -except  db aD aD_RDS aD_RPS Fit Fit_RDS Fit_RPS

Areas = {'V1','LM','RL'};

% Responsiveness criteria:
thr_std = 4;
thr_type = 1;
% thr_rep_fraction_rds = 0.5;
% thr_anova = 0.01;%0.01;
thr_Rsquared = 0.5;%0.5;

RDStype1  = 'aRDS';
RDStype2  = 'aRDS2';
RDStype12 = 'aRDS_aRDS2';

nAreas = numel(Areas);

phase_gabor_eachArea = cell(nAreas, 1);
amp_ratio_eachArea   = cell(nAreas, 1);
Rsquared_eachArea    = cell(nAreas, 1);
nExps_max = 9;
phase_gabor_offset_eachExp    = cell(nExps_max, nAreas);
phase_gabor_offset_medEachExp =  nan(nExps_max, nAreas);
amp_ratio_eachExp             = cell(nExps_max, nAreas);
amp_ratio_medEachExp          =  nan(nExps_max, nAreas);

for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    phase_gabor = [];
    amp_ratio   = [];
    Rsquared    = [];
    cnt = 0;
    
    nExps = length(Fit(ai).(RDStype12));

    for e = 1 : nExps
        
        respRDS_1  = above_RsquaredThreshold(Fit(ai).(RDStype12 )(e).TC, 0, 'Rsquared_meanTC_1');
        respRDS_2  = above_RsquaredThreshold(Fit(ai).(RDStype12 )(e).TC, 0, 'Rsquared_meanTC_2');
        respRDS    = intersect(respRDS_1, respRDS_2);
        rsqRDS_1   = above_RsquaredThreshold(Fit(ai).(RDStype12 )(e).TC, thr_Rsquared, 'Rsquared_meanTC_1');
        rsqRDS_2   = above_RsquaredThreshold(Fit(ai).(RDStype12 )(e).TC, thr_Rsquared, 'Rsquared_meanTC_2');
        selROIs    = mintersect(respRDS, rsqRDS_1, rsqRDS_2);
                
        disp([' Selected ROIs: ' num2str(length(selROIs)) ' / ' num2str(length(Fit(ai).(RDStype12 )(e).TC)) ]);
        
        phase_gabor_thisExp = [];
        amp_ratio_thisExp   = [];
        cnt2 = 0;
        
        for i = 1 : length(selROIs)
            
            nr = selROIs(i);
            cnt  = cnt  +1 ;
            cnt2 = cnt2 +1 ;
            
            phase_gabor(cnt,1) = ( Fit(ai).(RDStype12)(e).TC(nr).BestFitEsts(6) );
            phase_gabor(cnt,2) = ( Fit(ai).(RDStype12)(e).TC(nr).BestFitEsts(8) );
            amp_ratio(cnt,1)   = ( Fit(ai).(RDStype12)(e).TC(nr).BestFitEsts(7) / Fit(ai).(RDStype12)(e).TC(nr).BestFitEsts(2) );
            Rsquared(cnt,1)    = Fit(ai).(RDStype12)(e).TC(nr).GoF.Rsquared_meanTC_1;
            Rsquared(cnt,2)    = Fit(ai).(RDStype12)(e).TC(nr).GoF.Rsquared_meanTC_2;
            
            phase_gabor_thisExp(cnt2,1) = phase_gabor(cnt,1);
            phase_gabor_thisExp(cnt2,2) = phase_gabor(cnt,2);
            amp_ratio_thisExp(cnt2,1)   = amp_ratio(cnt,1);
        end
        
        if length(selROIs) > 0
            phase_gabor_offset_eachExp{e,ax}    = phase_gabor_thisExp(:,1) - phase_gabor_thisExp(:,2);
            phase_gabor_offset_medEachExp(e,ax) = median(phase_gabor_thisExp(:,1) - phase_gabor_thisExp(:,2));
            amp_ratio_eachExp{e,ax}    = amp_ratio_thisExp;
            amp_ratio_medEachExp(e,ax) = median(amp_ratio_thisExp);
        else
            phase_gabor_offset_eachExp{e,ax}    = [];
            phase_gabor_offset_medEachExp(e,ax) = nan;
            amp_ratio_eachExp{e,ax}    = [];
            amp_ratio_medEachExp(e,ax) = nan;
        end
    end
    
    phase_gabor_eachArea{ax} = phase_gabor;
    amp_ratio_eachArea{ax}   = amp_ratio;
    Rsquared_eachArea{ax}    = Rsquared;
    
end


%% Scatter plot phase offset vs. amplitude ratio (Fig. 8E)
% Before this section, run the section:
% Get phase offset and amplitude ratio for corr-RDS and anticorr RDS (to run for Figure 8E-G plots)

PhaseStep = 30;
phase_bin_centers = -90:PhaseStep:270;
phase_bin_edges = [phase_bin_centers,phase_bin_centers(end)+PhaseStep]-PhaseStep/2;
CC = load_DGD_colors(1); cx = 1;
marker_size = 31;
% yrange = [0 50]; % percentage

figure; hold on;
for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    phase_offset_thisArea = phase_gabor_eachArea{ax}(:,1) - phase_gabor_eachArea{ax}(:,2) ;
    phase_offset_thisArea (phase_offset_thisArea < -90) = phase_offset_thisArea (phase_offset_thisArea < -90) + 360;

    amp_ratio_thisArea = amp_ratio_eachArea{ax} ;
    
    scatter(phase_offset_thisArea, amp_ratio_thisArea, marker_size, CC.(Areas{ax}){cx}, 'filled');
    
end

xlim([phase_bin_edges(1) phase_bin_edges(end)])
set(gca, 'XTick',phase_bin_centers)
set(gca, 'YScale','log')
set(gca, 'YTick',[1:10], 'YTickLabel',[1:10])
ylabel('Amplitude ratio' );
xlabel('Phase offset (phase deg)');


%% Histogram phase offset frequency distribution (Fig. 8F)
% Before this section, run the section:
% Get phase offset and amplitude ratio for corr-RDS and anticorr RDS (to run for Figure 8E-G plots)

PhaseStep = 30;
phase_bin_centers = -90:PhaseStep:270;
phase_bin_edges = [phase_bin_centers,phase_bin_centers(end)+PhaseStep]-PhaseStep/2;
CC = load_DGD_colors(1); cx = 1;
yrange = [0 50]; % percentage

for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    phase_offset = phase_gabor_eachArea{ax}(:,1) - phase_gabor_eachArea{ax}(:,2) ;
    phase_offset (phase_offset < -90) = phase_offset (phase_offset < -90) + 360;
    
    Count_bins_phase_offset = 100*histcounts(phase_offset, phase_bin_edges, 'Normalization','probability');
    
    figure;
    hsb = superbar( phase_bin_centers, Count_bins_phase_offset,...
        ...'E', CountPrefIOPhase_norm_semAcrossExps_eachArea(ax,:),...
        'BarWidth',0.9*PhaseStep, 'BarFaceColor', CC.(Areas{ax}){cx}, 'BarEdgeColor','none',...
        'ErrorbarRelativeWidth',0.0, 'ErrorbarLineWidth',2);
    set(hsb, 'ShowBaseLine', 'off')

    xlim([phase_bin_edges(1) phase_bin_edges(end)])
    set(gca, 'XTick',phase_bin_centers)
    ylim(yrange)
    set(gca, 'YTick',[yrange(1):10:yrange(end)])
    
    ylabel('Proportion of cells (%)' );
    xlabel('Phase offset (phase deg)');
        
end

%% Histogram amplitude ratio frequency distribution (Fig. 8G)
% Before this section, run the section:
% Get phase offset and amplitude ratio for corr-RDS and anticorr RDS (to run for Figure 8E-G plots)

AmpStep = 1/4;
amp_ratio_bin_edges = [0:AmpStep:2,  Inf];
amp_ratio_bin_ticklabels = {'0', '0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', '2', 'Inf'};
amp_ratio_bin_centers = [1:length(amp_ratio_bin_edges)-1] - 0.5;
amp_ratio_bin_ticks = 0:amp_ratio_bin_centers(end)+0.5;

CC = load_DGD_colors(1); cx = 1;
yrange = [0 50]; % percentage

for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    amp_ratio_thisArea = amp_ratio_eachArea{ax} ;
    
    Count_bins_phase_offset = 100*histcounts(amp_ratio_thisArea, amp_ratio_bin_edges, 'Normalization','probability');
    
    figure;
    hsb = superbar( amp_ratio_bin_centers, Count_bins_phase_offset,...
        ...'E', CountPrefIOPhase_norm_semAcrossExps_eachArea(ax,:),...
        'BarWidth',0.9*1, 'BarFaceColor', CC.(Areas{ax}){cx}, 'BarEdgeColor','none',...
        'ErrorbarRelativeWidth',0.0, 'ErrorbarLineWidth',2);
    set(hsb, 'ShowBaseLine', 'off')

    xlim([amp_ratio_bin_ticks(1) amp_ratio_bin_ticks(end)])
    set(gca, 'XTick',amp_ratio_bin_ticks, 'XTickLabel',amp_ratio_bin_ticklabels)
    ylim(yrange)
    set(gca, 'YTick',[yrange(1):10:yrange(end)])
    
    ylabel('Proportion of cells (%)' );
    xlabel('Amplitude ratio');
    
end
