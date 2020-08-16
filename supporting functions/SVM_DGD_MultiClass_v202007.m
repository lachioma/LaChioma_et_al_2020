

function [svm, svm_CB] = ...
            SVM_DGD_MultiClass_v202007 ...
                (aDGD, thr_std, thr_type, nPredictors, nCombosPred, Directions, Split,...
                 nShuffles, NoiseCorrelationBlind,...
                 nCVPpermutations, CombosStims, StimNames, show_fig, RoiNrsToUseAll_ALL)

% [scriptPath, scriptName] = fileparts(mfilename('fullpath'));


if isempty(nPredictors)
    nPredictors =  [10] ;
end
if isempty(nCombosPred)
    % nr. combinations per nr. of predictors used:
    nCombosPred = 20;
end
if ~exist('Directions','var') || isempty(Directions)
    Directions = [1];
end
if ~exist('Split','var') || isempty(Split)
    Split = 0;
end
if ~exist('nShuffles','var') || isempty(nShuffles)
    nShuffles = 0;
end
if ~exist('NoiseCorrelationBlind','var') || isempty(NoiseCorrelationBlind)
    % For each cell, shuffle the trials within each stimulus, to remove the
    % correlations across cells.
    NoiseCorrelationBlind = 0;
end
if ~exist('show_fig','var') || isempty(show_fig)
    show_fig = 1;
end
if ~exist('RoiNrsToUseAll_ALL','var')
    RoiNrsToUseAll_ALL = [];
end
if ~exist('nCVPpermutations','var') || isempty(nCVPpermutations)
    nCVPpermutations = 1;
end
if iscell(aDGD.StimSettings)
    StimSettings = aDGD.StimSettings{1};
    Param        = aDGD.Param{1};
else
    StimSettings = aDGD.StimSettings;
    Param        = aDGD.Param;
end
if ~exist('CombosStims','var') || isempty(CombosStims)
    CombosStims = [1 : StimSettings.nPhases];
end
if ~exist('StimNames','var') || isempty(StimNames)
    StimNames = cellstr(num2str(StimSettings.IOPhaseDifferences(CombosStims)'));
end
StimNames = strtrim(StimNames);
svm_CB = [];

CrossVal = 'on';
nFolds   =  10 ;
GetConfusionMatrix = 1;
     
KS = [ 0 ];
BC = [ 0 ];
% KS = log10( linspace(1,1000,7) );
% BC = log10( linspace(1,10000,7) );
% KS = log10( linspace(10^(0),10^(3),4) );
% BC = log10( linspace(10^(0),10^(4),5) );
cKSBC = AllCombosOf2Vectors(KS,BC);
Nr_cKSBC = size(cKSBC,1);

    

if Split > 1
    aDGD = Bin_Fstim_eachTrial(aDGD, Param, Split);
else
    [aDGD.ROIs.mean_Fstim_eachTrial2] = aDGD.ROIs.mean_Fstim_eachTrial;
end

% For reproducibility:
% rng('default');
rng(4);

% nTotNeurons = length(aDGD.ROIs);
% nReps    = size(aDGD.ROIs(1).mean_Fstim_eachTrial,3);
% Combinations of neighboring stims:
% CombosStims = [1:nStims; circshift(1:nStims,[0 -1])]'; 
nStims = size(CombosStims,2);
Nr_CombosStims = size(CombosStims,1);
Nr_Directions          = numel(Directions);
Nr_Directions_perGroup = size(Directions,2);
Nr_Groups              = size(Directions,1);
angles = StimSettings.angles_cartesian(Directions(:)); angles = reshape(angles, Nr_Groups,Nr_Directions_perGroup);
% nTrials = nReps*Nr_Directions_perGroup;

% Total number of runs, excluding the cKSBC:
TotNrRuns = length(nPredictors)*nCombosPred*Nr_CombosStims*Nr_Groups;

if show_fig
    hf  = figure;
        fig_pos = get(gcf, 'position'); fig_pos(3) = fig_pos(3)*Nr_Groups/1.5;
        set(gcf, 'position',fig_pos); % [left, bottom, width, height]
    cmap = colormap('jet'); ix = round( linspace(1,64,Nr_Groups)); Colors = cmap(ix,:);
    acircle = 40;
end


ticstart = tic;
multiWaitbar( 'CloseAll' );
cnt = 0;
multiWaitbar( ['Overall progress = ' num2str(cnt) '/' num2str(TotNrRuns)], 0, 'Color', [0.4 0.1 0.5] );

ER_meanAcrossStimPairs_perDir_acrossComboPred = nan(length(nPredictors), Nr_Groups);
ER_stdAcrossStimPairs_perDir_acrossComboPred  = nan(length(nPredictors), Nr_Groups);
ER_semAcrossStimPairs_perDir_acrossComboPred  = nan(length(nPredictors), Nr_Groups);
if NoiseCorrelationBlind
    ER_meanAcrossStimPairs_perDir_acrossComboPred_CB = nan(length(nPredictors), Nr_Groups);
    ER_stdAcrossStimPairs_perDir_acrossComboPred_CB  = nan(length(nPredictors), Nr_Groups);
    ER_semAcrossStimPairs_perDir_acrossComboPred_CB  = nan(length(nPredictors), Nr_Groups);
end
ER_sh_5prctile_mean_perDir_AcrossComboPred    = nan(length(nPredictors), Nr_Groups);
ER_sh_1prctile_mean_perDir_AcrossComboPred    = nan(length(nPredictors), Nr_Groups);

for d = 1 : Nr_Groups
    RoiNrs{d} = above_stdThreshold_DGD(aDGD,  thr_std, thr_type, angles(d,:));
    nNeurons(d) = length(RoiNrs{d});
end


for np = 1:length(nPredictors)
    
    svm(np).nCombosPred = nCombosPred;
    svm(np).CombosStims = CombosStims;
    svm(np).Nr_CombosStims = Nr_CombosStims;
    svm(np).Directions = angles;
    svm(np).Nr_Directions = Nr_Directions;
    svm(np).Nr_Directions_perGroup = Nr_Directions_perGroup;
    svm(np).Nr_Groups = Nr_Groups;
    svm(np).cKSBC = cKSBC;
    svm(np).Nr_cKSBC = Nr_cKSBC;
    svm(np).nShuffles = nShuffles;
    svm(np).NoiseCorrelationBlind = 0;
    svm(np).Split = Split;
    %
    svm(np).ConfMat_perc_perComboPred = nan(nStims,nStims,Nr_Groups,nCombosPred);
    %
    if NoiseCorrelationBlind
        svm_CB(np).nCombosPred = nCombosPred;
        svm_CB(np).CombosStims = CombosStims;
        svm_CB(np).Nr_CombosStims = Nr_CombosStims;
        svm_CB(np).Directions = angles;
        svm_CB(np).Nr_Directions = Nr_Directions;
        svm_CB(np).Nr_Directions_perGroup = Nr_Directions_perGroup;
        svm_CB(np).Nr_Groups = Nr_Groups;
        svm_CB(np).cKSBC = cKSBC;
        svm_CB(np).Nr_cKSBC = Nr_cKSBC;
        svm_CB(np).nShuffles = nShuffles;
        svm_CB(np).NoiseCorrelationBlind = 1;
        svm_CB(np).Split = Split;
        svm_CB(np).ConfMat_perc_perComboPred = nan(nStims,nStims,Nr_Groups,nCombosPred);
    end
    
    
    ER_meanAcrossStimPairs_perDir_perComboPred = nan(nCombosPred, Nr_Groups);
    if NoiseCorrelationBlind
        ER_meanAcrossStimPairs_perDir_perComboPred_CB = nan(nCombosPred, Nr_Groups);
    end
    if nShuffles > 0
        ER_sh_5prctile = nan(Nr_CombosStims, nCombosPred, Nr_Groups);
        ER_sh_1prctile = nan(Nr_CombosStims, nCombosPred, Nr_Groups);
    end
    
    for d = 1 : Nr_Groups
        
        nPredictors_toUse = nPredictors(np);

        if nPredictors_toUse > nNeurons(d)
            continue
        end
        
        if nPredictors_toUse == nNeurons(d)
            nCombosPred_toUse(d) = 1;
        elseif nPredictors_toUse == nNeurons(d)-1 || nPredictors_toUse == nNeurons(d)-2
            nCombosPred_toUse(d) = min(nCombosPred, nchoosek(nNeurons(d),nPredictors_toUse));
        else
            nCombosPred_toUse(d) = nCombosPred;
        end
        
        
        if isempty(RoiNrsToUseAll_ALL)
            RoiNrsToUseAll = zeros(nCombosPred_toUse(d), nPredictors_toUse);
            % The following is in a separate for-loop because IxCellsAll is not
            % compatible with parfor:
            for nc = 1 : nCombosPred_toUse(d)
                RedoRandPerm = 1;
                while RedoRandPerm == 1
                    IxCells = sort(randperm(nNeurons(d), nPredictors_toUse));
                    if ~ismember(IxCells, RoiNrsToUseAll, 'rows')
                        RedoRandPerm = 0;
                        RoiNrsToUseAll(nc,:) = RoiNrs{d}(IxCells);
                    end
                end
            end
        else
            RoiNrsToUseAll = RoiNrsToUseAll_ALL;
        end
        
        svm(np).nCombosPred_toUse = nCombosPred_toUse;
        svm(np).nPredictors       = nPredictors_toUse;
        svm(np).RoiNrsIncluded{d} = RoiNrsToUseAll;
        if NoiseCorrelationBlind
            svm_CB(np).nCombosPred_toUse = nCombosPred_toUse;
            svm_CB(np).nPredictors       = nPredictors_toUse;
            svm_CB(np).RoiNrsIncluded{d} = RoiNrsToUseAll;
        end
    
        ER           = nan(Nr_CombosStims, nCombosPred, 1);
        if NoiseCorrelationBlind
            ER_CB    = nan(Nr_CombosStims, nCombosPred, 1);
        end
        if nShuffles > 0
            ER_Sh    = nan(Nr_CombosStims, nCombosPred, nShuffles);
        end
        
        for nc = 1 : nCombosPred_toUse(d)
            
            RoiNrsUsed = RoiNrsToUseAll(nc,:);
            
            multiWaitbar( ['Nr_Groups = ' num2str(Nr_Groups)], d/Nr_Groups, 'Color', [0.9 0.8 0.2] );
            multiWaitbar( ['nPredictors = ' num2str(nPredictors)], np/length(nPredictors_toUse), 'Color', [0.8 0.0 0.1] );
            multiWaitbar( ['nCombosPred = ' num2str(nCombosPred_toUse(d))], nc/nCombosPred_toUse(d), 'Color', [1.0 0.4 0.0] );
            
            for cs = 1 : Nr_CombosStims
                
                multiWaitbar( ['Nr_CombosStims = ' num2str(Nr_CombosStims)], cs/Nr_CombosStims, 'Color', [0.2 0.9 0.3] );
                
                cdata = []; cdata_CB = [];
                label = [];
                
                for i = 1:length(RoiNrsUsed)
                        DirIx = Directions(d,:);
                        nr = RoiNrsUsed(i);
                        data_roi = [];  data_roi_CB = [];
                        
                        for s = 1 : nStims
                            
                            six = CombosStims(cs,s);
                            
                            data_tmp = squeeze( aDGD.ROIs(nr).mean_Fstim_eachTrial2(DirIx,six,:) );
                            data_tmp = reshape(data_tmp', [],1);
                            data_roi  = [data_roi; data_tmp];
                            
                            if i==1
                                label_tmp = repmat(StimNames(six), length(data_tmp), 1);
                                label = [label; label_tmp];
                            end
                            
                            if NoiseCorrelationBlind
                                data_tmp_CB = data_tmp( randperm(length(data_tmp)) );
                                data_roi_CB = [data_roi_CB; data_tmp_CB];
                            end
                            
                        end
                        
                        
                        cdata(:,i) = data_roi;
                        
                        if NoiseCorrelationBlind
                            cdata_CB(:,i) = data_roi_CB;
                        end
                        
                end
                
                MisclassRate = nan(Nr_cKSBC, 1);
                
                multiWaitbar( ['Overall progress = ' num2str(cnt) '/' num2str(TotNrRuns)], (cnt+1)/TotNrRuns, 'Relabel', ['Overall progress = ' num2str(cnt+1) '/' num2str(TotNrRuns)] );
                cnt = cnt + 1;
                multiWaitbar( ['Nr_cKSBC = ' num2str(Nr_cKSBC)], 'Busy', 'Color', [0.1 0.5 0.8] );
                % this is alos to prevent parfor warning of braodcast
                % variable:
                KSs = 10.^(cKSBC(:,1));
                BCs = 10.^(cKSBC(:,2));
                % For fitecoc:
                cdata = cdata';
                if NoiseCorrelationBlind
                    cdata_CB = cdata_CB';
                end
                
                ConfMat_perc_eachcKSBC = nan(nStims,nStims, Nr_cKSBC);
                
                for j = 1 : Nr_cKSBC
                    
                    KSj = KSs(j);
                    BCj = BCs(j);
                    
                    [MisclassRate(j),~,~,ConfMat_perc_eachcKSBC(:,:,j)] = ...
                        SVM_Perform_MultiClass( cdata, label, KSj,BCj, CrossVal, nFolds, nCVPpermutations, StimNames, GetConfusionMatrix ); %#ok<*AGROW>
                    
                    disp(['  --comboKS-BC ' num2str(j) '/' num2str(Nr_cKSBC) '  -->  ErrorRate = ' num2str(MisclassRate(j)) ]);

                end
                
                multiWaitbar( ['Nr_cKSBC = ' num2str(Nr_cKSBC)], 1, 'Color', [0.1 0.5 0.8] );
                

                [ER(cs,nc), ~] = min(MisclassRate);
                IxBestClass = find(MisclassRate==min(MisclassRate),1,'last');
                
                KSbest = 10.^[cKSBC(IxBestClass,1)];
                BCbest = 10.^[cKSBC(IxBestClass,2)];
                svm(np).KSbest(nc,cs,d) = KSbest;
                svm(np).BCbest(nc,cs,d) = BCbest;
                
                ConfMat_perc = ConfMat_perc_eachcKSBC(:,:,IxBestClass);
                svm(np).ConfMat_perc_perComboPred(:,:,d,nc) = ConfMat_perc; %(nClasses,nClasses,d,nc)

                
                
                
                if NoiseCorrelationBlind
                    
                    [MisclassRate_CB,~,~,ConfMat_perc_CB] = ...
                        SVM_Perform_MultiClass( cdata_CB, label, KSbest,BCbest, CrossVal, nFolds, nCVPpermutations, StimNames, GetConfusionMatrix ); %#ok<*AGROW>
                    
                    svm_CB(np).KSbest(nc,cs,d) = KSbest;
                    svm_CB(np).BCbest(nc,cs,d) = BCbest;
                    ER_CB(cs,nc) = MisclassRate_CB;
                    svm_CB(np).ConfMat_perc_perComboPred(:,:,d,nc) = ConfMat_perc_CB; %(nClasses,nClasses,d,nc)

                end
                
                
                
                if nShuffles > 0
                    
                    disp(['  ... performing n=' num2str(nShuffles) ' shuffle controls ...'])
                    
                    parfor ns = 1 : nShuffles
                        
                        label_sh = label( randperm(length(label)) );
                        MisclassRate_Sh_tmp = nan(Nr_cKSBC,1);

                        for j = 1 : Nr_cKSBC
                    
                            KSj = KSs(j);
                            BCj = BCs(j);
                            MisclassRate_Sh_tmp(j) = SVM_Perform_MultiClass( cdata, label_sh, KSj,BCj, CrossVal, nFolds, nCVPpermutations, StimNames, 0 );
                            
                        end
                        
                        MisclassRate_Sh = min(MisclassRate_Sh_tmp);
                        ER_Sh(cs,nc,ns) = MisclassRate_Sh;
                        
                    end
                    
                    disp(['  ... 95th percentile shuffle control = ' num2str( prctile(ER_Sh(cs,nc,:),5) )])
                    
                end
                
                
            end

        end
        
        ER_meanAcrossStimPairs_perDir_perComboPred(:,d) = nanmean(ER,1); % (nc,d)
        svm(np).ErrRate(:,:,d) = ER;
        
        if NoiseCorrelationBlind
            ER_meanAcrossStimPairs_perDir_perComboPred_CB(:,d) = nanmean(ER_CB,1); % (nc,d)
            svm_CB(np).ErrRate(:,:,d) = ER_CB;
        end
                
        if nShuffles > 0
            svm(np).ErrRate_Sh(:,:,:,d) = ER_Sh; %(cs,nc,ns)
            if NoiseCorrelationBlind
                svm_CB(np).ErrRate_Sh(:,:,:,d) = ER_Sh;
            end
            svm(np).ER_sh_5prctile(:,:,d) = prctile(ER_Sh, 5, 3); %(cs,nc,d)
            svm(np).ER_sh_1prctile(:,:,d) = prctile(ER_Sh, 1, 3); %(cs,nc,d)
        end


    end
    
    
    ER_meanAcrossStimPairs_perDir_acrossComboPred(np,:) = nanmean(ER_meanAcrossStimPairs_perDir_perComboPred,1); % (np,d)
    ER_stdAcrossStimPairs_perDir_acrossComboPred (np,:) = nanstd (ER_meanAcrossStimPairs_perDir_perComboPred,0,1); % (np,d)
    ER_semAcrossStimPairs_perDir_acrossComboPred (np,:) = ER_stdAcrossStimPairs_perDir_acrossComboPred (np,:) ./ sqrt(nCombosPred_toUse); % (np,d)
    svm(np).ER_meanAcrossStimPairs_perDir_perComboPred    = ER_meanAcrossStimPairs_perDir_perComboPred;
    svm(np).ER_meanAcrossStimPairs_perDir_acrossComboPred = ER_meanAcrossStimPairs_perDir_acrossComboPred(np,:);
    svm(np).ER_stdAcrossStimPairs_perDir_acrossComboPred  = ER_stdAcrossStimPairs_perDir_acrossComboPred (np,:);
    svm(np).ER_semAcrossStimPairs_perDir_acrossComboPred  = ER_semAcrossStimPairs_perDir_acrossComboPred (np,:);

    svm(np).ConfMat_perc_acrossComboPred = nanmean(svm(np).ConfMat_perc_perComboPred, 4); %(:,:,d)
    
    
    if NoiseCorrelationBlind
        ER_meanAcrossStimPairs_perDir_acrossComboPred_CB(np,:) = nanmean(ER_meanAcrossStimPairs_perDir_perComboPred_CB,1); % (np,d)
        ER_stdAcrossStimPairs_perDir_acrossComboPred_CB (np,:) = nanstd (ER_meanAcrossStimPairs_perDir_perComboPred_CB,0,1); % (np,d)
        ER_semAcrossStimPairs_perDir_acrossComboPred_CB (np,:) = ER_stdAcrossStimPairs_perDir_acrossComboPred (np,:) ./ sqrt(nCombosPred_toUse); % (np,d)
        svm_CB(np).ER_meanAcrossStimPairs_perDir_perComboPred    = ER_meanAcrossStimPairs_perDir_perComboPred_CB;
        svm_CB(np).ER_meanAcrossStimPairs_perDir_acrossComboPred = ER_meanAcrossStimPairs_perDir_acrossComboPred_CB(np,:);
        svm_CB(np).ER_stdAcrossStimPairs_perDir_acrossComboPred  = ER_stdAcrossStimPairs_perDir_acrossComboPred_CB (np,:);
        svm_CB(np).ER_semAcrossStimPairs_perDir_acrossComboPred  = ER_semAcrossStimPairs_perDir_acrossComboPred_CB (np,:);
        svm_CB(np).ConfMat_perc_acrossComboPred = nanmean(svm_CB(np).ConfMat_perc_perComboPred, 4); %(:,:,d)
    end
    
    if nShuffles > 0
        ER_sh_5prctile_mean_perDir_AcrossComboPred(np,:) = nanmean(nanmean(svm(np).ER_sh_5prctile,1),2); %(np,d)
        ER_sh_1prctile_mean_perDir_AcrossComboPred(np,:) = nanmean(nanmean(svm(np).ER_sh_1prctile,1),2); %(np,d)
        svm(np).ER_sh_5prctile_mean_perDir_AcrossComboPred = ER_sh_5prctile_mean_perDir_AcrossComboPred(np,:); %(1,d)
        svm(np).ER_sh_1prctile_mean_perDir_AcrossComboPred = ER_sh_1prctile_mean_perDir_AcrossComboPred(np,:); %(1,d)
        if NoiseCorrelationBlind
            svm_CB(np).ER_sh_5prctile  = ER_sh_5prctile; %(cs,nc,d)
            svm_CB(np).ER_sh_1prctile  = ER_sh_1prctile; %(cs,nc,d)
            svm_CB(np).ER_sh_5prctile_mean_perDir_AcrossComboPred = ER_sh_5prctile_mean_perDir_AcrossComboPred(np,:); %(1,d)
            svm_CB(np).ER_sh_1prctile_mean_perDir_AcrossComboPred = ER_sh_1prctile_mean_perDir_AcrossComboPred(np,:); %(1,d)
        end
    end
    
    
    
    if show_fig
        
        figure(hf)
        for d = 1:Nr_Groups
            subaxis(1,Nr_Groups,d, 'ML',0.05,'MR',0.03, 'SH',0.03);
            cla; hold on;
            if np==1
                scatter( [1:np], 1-ER_meanAcrossStimPairs_perDir_acrossComboPred([1:np],d), acircle,Colors(d,:),'filled' );
                errorbar([1:np], 1-ER_meanAcrossStimPairs_perDir_acrossComboPred([1:np],d),...
                                  ER_semAcrossStimPairs_perDir_acrossComboPred([1:np],d), 'Color',Colors(d,:) );
                if nShuffles > 0
                scatter( [1:np], 1-ER_sh_5prctile_mean_perDir_AcrossComboPred([1:np],d), acircle/3, 'k', 'filled');
                end
            else
                scatter      ( [1:np], 1-ER_meanAcrossStimPairs_perDir_acrossComboPred([1:np],d), acircle,Colors(d,:),'filled' );
                plot         ( [1:np], 1-ER_meanAcrossStimPairs_perDir_acrossComboPred([1:np],d), 'Color',Colors(d,:) );
                AddErrorArea ( [1:np], 1-ER_meanAcrossStimPairs_perDir_acrossComboPred([1:np],d),...
                                  ER_semAcrossStimPairs_perDir_acrossComboPred([1:np],d), Colors(d,:) );
                if nShuffles > 0
                scatter( [1:np], 1-ER_sh_5prctile_mean_perDir_AcrossComboPred([1:np],d), acircle/3, 'k', 'filled');
                plot   ( [1:np], 1-ER_sh_5prctile_mean_perDir_AcrossComboPred([1:np],d), 'Color','k' );
                end
            end
            ylim([0 1]);
            xlim([0.5 np+0.5])
            set(gca, 'YTick', [0:0.2:1], 'YTickLabel', [])
            set(gca, 'XTick', [1:np], 'XTickLabel', nPredictors(1:np))
            title([' direction ' num2str(angles(d,:))])
            if d==1
                set(gca, 'YTick', [0:0.2:1], 'YTickLabel', [0:0.2:1])
                xlabel('Nr. neurons included')
                ylabel('Classification accuracy')
            end
        end
    end
    
    
end

multiWaitbar( 'CloseAll' );

disp( '   ==== done! ====')
disp(['   ==== total time elapsed ' num2str(toc(ticstart)/60) ' min'])     
    


% if exist('scriptName','var')
%     BackupRunningScript( [], scriptPath, scriptName )
% end
