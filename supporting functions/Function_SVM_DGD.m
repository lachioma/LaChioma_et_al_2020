function Function_SVM_DGD_appended_v202007(aDGD, SF, thr_type, thr_std, ExpIDs,...
            nPredictors, nCombosPred, Directions, Split,...
            nShuffles, NoiseCorrelationBlind, nCVPpermutations, CombosStims, Method, SaveFolder)

% This function is basically like Script_SVM_DGD_appended but as a
% function, so it can be run by multiple Matlab sessions with different
% parameters.

if ~exist('SF','var') || isempty(SF)
    SF = '0.05';
end
Area = aDGD.Area;
if ~exist('thr_type','var') || isempty(thr_type)
    thr_type = 1;
end
if ~exist('thr_std','var') || isempty(thr_std)
    thr_std = 4;
end
if ~exist('nPredictors','var') || isempty(nPredictors)
    nPredictors = [2, 5, 10, 20, 40, 80 ] ;
end
if ~exist('nCombosPred','var') || isempty(nCombosPred)
    nCombosPred =  10;
end
if ~exist('Directions','var') || isempty(Directions)
    Directions = [1:4];
end
if ~exist('Split','var') || isempty(Split)
    Split = 0;
end
if ~exist('ExpIDs','var') || isempty(ExpIDs)
    ExpIDsall = [aDGD.ROIs.ExpID];
    ExpIDs = unique(ExpIDsall);
    % ExpIDs = [8,9];
end
if ~exist('nShuffles','var') || isempty(nShuffles)
    nShuffles = 0;
end
if ~exist('NoiseCorrelationBlind','var') || isempty(NoiseCorrelationBlind)
    % For each cell, shuffle the trials within each stimulus, to remove the
    % noise correlations across cells.
    NoiseCorrelationBlind = 0;
end
plot_bestKS_BC = 0;

if ~exist('Method','var') || isempty(Method)
    Method = 2;
end
if ~exist('SaveFolder','var')
    SaveFolder = pwd;
end

ticstart = tic;
varfilenameAllExps = [SaveFolder filesep 'svm' '_' Area '_' 'SF' num2str(SF,'%3.2f') 'cpd' '_AllExps_' datestr(now,'yyyy-mm-dd_HH-MM-SS') '.mat'];
ee = 0;

for e = 1:length(ExpIDs)
    
    expid = ExpIDs(e);
    Info = aDGD.Info{e};
%     exptag = Info.Exp(1:5);
    ExpIDRoiNrs = find([aDGD.ROIs.ExpID] == expid);
    
    SpatFreq = str2double(SF);
    SFRoiNrs = find(ismember([aDGD.ROIs(ExpIDRoiNrs).SpatFreq], SpatFreq));
    if isempty(SFRoiNrs)
        disp(['  - ExpID = ' num2str(expid) ' has SF = ' num2str(aDGD.ROIs(ExpIDRoiNrs(1)).SpatFreq) ', instead of the requested ' SF ', so this exp will be excluded.'])
        expExcluded_skipforloop = 1;
        continue
    else
        expExcluded_skipforloop = 0;
        ee = ee + 1;
    end
    
    % Prepare variable aDGDtemp:
    aDGDtmp.ROIs  = aDGD.ROIs(ExpIDRoiNrs);
    aDGDtmp.Param = aDGD.Param(e);
    aDGDtmp.StimSettings = aDGD.StimSettings(e);

    switch Method 
        case 1
            
            show_fig = 0;
            
            [svm] = ...
                SVM_DGD_MultiClass_v202007 ...
                    (aDGDtmp, thr_std, thr_type, nPredictors, nCombosPred, Directions, Split, nShuffles, NoiseCorrelationBlind, nCVPpermutations, CombosStims, [], show_fig);

    end
    
    varfilename = [SaveFolder filesep 'svm' '_' Area '_' 'SF' num2str(SF,'%3.2f') 'cpd' '_expid' num2str(expid) '_' datestr(now,'yyyy-mm-dd_HH-MM-SS') '.mat'];
    save(varfilename, 'svm',...
        'SF','thr_type','thr_std','thr_type',...
        'nPredictors', 'nCombosPred', 'Directions', 'ExpIDs', 'nShuffles', 'NoiseCorrelationBlind');

    switch plot_bestKS_BC case 1
        best=[];
        best(:,1) = log10(reshape([svm.KSbest],[],1));
        best(:,2) = log10(reshape([svm.BCbest],[],1));
        [xyUnique, ignore, ixs] = unique(best,'rows');
        [nRows, nCols] = size(xyUnique);
        xyCount = hist(ixs,nRows);
        figure;
        scatter(...
            xyUnique(:,1), ...  %X values
            xyUnique(:,2), ...  %Y values
            xyCount*30, ...     %Individual marker sizes, note scale factor to make this visible
            'b'); ...           %Marker colors
        xlabel('KernelScale'); ylabel('BoxConstraint');    
    end

    svmAll{ee,1}  = svm;
%     svm2All{ee,1} = svm2;
    save(varfilenameAllExps, 'svmAll',...
        'SF','thr_type','thr_std','thr_type',...
        'nPredictors', 'nCombosPred', 'Directions', 'ExpIDs', 'nShuffles', 'NoiseCorrelationBlind');
end

disp(' ')
disp('   ------- All exps done ! --------');
disp(' ')
disp(['  Total time elapsed = ' num2str(toc(ticstart)/60/60) ' hours.']);
disp(' ')
