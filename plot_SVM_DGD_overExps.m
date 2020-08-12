function S = plot_SVM_DGD_overExps( Area, SF, folder )


if nargin < 3 || isempty(folder)
    folder = [ pwd filesep 'svm_' Area '_DGD' filesep ];
end
filename = ['svm_' Area '_SF' SF 'cpd_expid*' '.mat'];
list = dir( [ folder filesep filename ]);
nExps = length(list);
% max_nPredictors = 4;
% nDirections = 4;

% ER_meanAcrossStimPairs_perDir_acrossComboPred = nan(max_nPredictors, nDirections, nExps);

Predictors  = cell(nExps,1);
nPredictors_plot = nan(nExps,1);
nPredictors_count = [];
% ER_MeanPerExp = nan(nExps,1);
% ER_SemPerExp  = nan(nExps,1);
% ER_sh_5prctileAcrossStimPairs_AcrossComboPred = nan(nExps,1);
ER_sh = cell(1);

for e = 1:nExps
    
    % load variable svm:
    warning('off','MATLAB:load:classNotFound');
        load([folder filesep list(e).name], 'svm');
    warning('on','MATLAB:load:classNotFound');
    
    Predictors{e} = [svm.nPredictors];
    nPredictors_plot(e) = length(Predictors{e});
    
    if e == 1
        try
            nDirections = svm.Nr_Groups;
        catch
            nDirections = svm.Nr_Directions;
        end
        ER_MeanPerExp = nan(nExps,1,nDirections);
        ER_SemPerExp  = nan(nExps,1,nDirections);
%         ER_sh_5prctileAcrossStimPairs_AcrossComboPred = nan(nExps,1,nDirections);
    end
    
        Predictors{e} = [svm.nPredictors];
        nPredictors_plot(e) = length(Predictors{e});
        if length(nPredictors_count) < nPredictors_plot(e)
            extendVariable = 1;
            nPredictors_count(end+1:nPredictors_plot(e)) = 0;
            ER_MeanPerExp(:,end+1:nPredictors_plot(e),:) = nan;
            ER_SemPerExp (:,end+1:nPredictors_plot(e),:) = nan;
        else extendVariable = 0;
        end
        nPredictors_count(1:nPredictors_plot(e)) = nPredictors_count(1:nPredictors_plot(e)) + 1;
            
        ER_meanAcrossStimPairs_perDir_acrossComboPred = reshape([svm.ER_meanAcrossStimPairs_perDir_acrossComboPred], nDirections,[])'; %(nP,dir)
        ER_MeanPerExp(e,:,1:nDirections) = ER_meanAcrossStimPairs_perDir_acrossComboPred;
        
        ER_semAcrossStimPairs_perDir_acrossComboPred = reshape([svm.ER_semAcrossStimPairs_perDir_acrossComboPred], nDirections,[])'; %(nP,dir)
        ER_SemPerExp(e,:,1:nDirections) = ER_semAcrossStimPairs_perDir_acrossComboPred;

        for pr = 1 : length(svm)
            if e==1, ER_sh{pr} = []; end;
            ER_sh{pr} = [ER_sh{pr}; svm(pr).ErrRate_Sh(:)];
        end
end
    
    nPredictors_count = repmat(nPredictors_count, nDirections,1);

    ER_meanAcrossExps = shiftdim(nanmean(ER_MeanPerExp,1))';
    ER_stdAcrossExps  = shiftdim( nanstd(ER_MeanPerExp,0,1))';
    ER_semAcrossExps  = ER_stdAcrossExps ./ sqrt(nPredictors_count);
    
    ER_sh_allAcrossExps_1p = nan(1,length(ER_sh));
    ER_sh_allAcrossExps_01p = nan(1,length(ER_sh));
    for pr = 1 : length(ER_sh)
        ER_sh_allAcrossExps_1p(pr)  = prctile(ER_sh{pr}, 1.0/2);
        ER_sh_allAcrossExps_01p(pr) = prctile(ER_sh{pr}, 0.1/2);
    end
%%
CC = load_DGD_colors(1); cx = 1;

    hf  = figure;
        fig_pos = get(gcf, 'position'); fig_pos(3) = fig_pos(3)*nDirections/1.5;
        set(gcf, 'position',fig_pos); % [left, bottom, width, height]
    max_nPredictors = max(nPredictors_plot);
    Predictors_all = Predictors{nPredictors_plot==max_nPredictors};

    ColorMean = CC.(Area){cx};
    ColorSingleExps = [1 1 1]*0.7;
    acircle = 50;
    angles_cartesian = svm.Directions;
%     x_vals = [1:max_nPredictors];
    x_vals = Predictors{e};
    for d = 1 : nDirections
        subaxis(1,nDirections, d, 'ML',0.05,'MR',0.02, 'SH',0.01);
        hold on;
        for e = 1:nExps
%             scatter      ( x_vals, 1-ER_MeanPerExp(e,:,d), acircle/2,ColorSingleExps,'filled' );
            plot         ( x_vals, 1-ER_MeanPerExp(e,:,d), 'Color',ColorSingleExps );
%             AddErrorArea ( x_vals, 1-ER_MeanPerExp(e,:,d),...
%                               ER_SemPerExp(e,:,d), ColorSingleExps, 0.3 );
        end
%             scatter      ( x_vals, 1-ER_meanAcrossExps(d,:), acircle,ColorMean,'filled' );
%             p1 = plot    ( x_vals, 1-ER_meanAcrossExps(d,:), 'Color',ColorMean, 'LineWidth',2.5 );
            AddErrorArea ( x_vals, 1-ER_meanAcrossExps(d,:),...
                              ER_semAcrossExps(d,:), ColorMean, 0.5 );
    %         scatter      ( x_vals, 1-ER_sh_5prctile_meanAcrossExps, acircle/2,[1 1 1]*0.3,'filled' );
    %         p2 = plot    ( x_vals, 1-ER_sh_5prctile_meanAcrossExps, 'Color',[1 1 1]*0.3, 'LineStyle',':','LineWidth',1 );

        plot ( x_vals, 1-ER_sh_allAcrossExps_01p, 'Color',ColorSingleExps, 'LineStyle',':' );
    
        label_nExps={}; for i=1:max_nPredictors, label_nExps{i} = num2str(nPredictors_count(d,i)); end;
        text( x_vals, repmat(0.05,max_nPredictors,1), label_nExps,...
            'HorizontalAlignment','center', 'FontSize',7);
        text( 1, 0.07, 'nExps', 'HorizontalAlignment','center', 'FontSize',7);
        ylim([0 1]);
        xlim( [1 , max(x_vals)+1] )
        box off
        set(gca, 'YTick', [0:0.2:1], 'YTickLabel', [])
        set(gca, 'XTick', x_vals, 'XTickLabel', Predictors_all)
        title([' direction ' num2str(angles_cartesian(d))])
        if d == 1
            set(gca, 'YTick', [0:0.2:1], 'YTickLabel', [0:0.2:1])
            xlabel('Nr. neurons included')
            ylabel('Classification accuracy')
%             legend(p1,'Mean across exps', 'Location','best')
        end
    end
%     mtit(['Area ' Area ' - ' 'DGD' ' - SF ' SF 'cpd'])

    
    
    
    
%%
S.Area               = Area;
S.SF                 = SF;
S.nExps              = nExps;
S.angles_cartesian   = angles_cartesian;
S.Predictors_all     = Predictors_all;
S.nPredictors_count  = nPredictors_count;
S.max_nPredictors    = max_nPredictors;
S.ER_MeanPerExp      = ER_MeanPerExp;
S.ER_SemPerExp       = ER_SemPerExp;
S.ER_meanAcrossExps  = ER_meanAcrossExps;
S.ER_semAcrossExps   = ER_semAcrossExps;

% if strcmp(SF,'0.05'), sx=1; elseif strcmp(SF,'0.10'), sx=2; elseif strcmp(SF,'0.01'), sx=3; end
% S.(Area)(sx).ER_MeanPerExp      = ER_MeanPerExp;
% S.(Area)(sx).ER_SemPerExp       = ER_SemPerExp;
% S.(Area)(sx).ER_meanAcrossExps  = ER_meanAcrossExps;
% S.(Area)(sx).ER_semAcrossExps   = ER_semAcrossExps;
% S.(Area)(sx).angles_cartesian  = angles_cartesian;
% S.(Area)(sx).max_nPredictors   = max_nPredictors;
% S.(Area)(sx).Predictors_all    = Predictors_all;


