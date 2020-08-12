function [MisclassRate, SVMModel, MisclassRate_eachPerm, ConfMat_perc] = ...
    SVM_Perform_MultiClass( cdata, label, KS,BC, CrossVal, nFolds, nCVPpermutations, ClassNames, GetConfusionMatrix )


% cdata must be in the matrix form (features, observations).
% Having this format gives faster computation for linear models.


nObservations = size(cdata,2);
nNeurons = size(cdata,1);

if ~exist('KS','var') || isempty(KS)
    KS = 1; % default=1
end
if ~exist('BC','var') || isempty(BC)
    BC = 1; % default=1
end
if ~exist('CrossVal','var') || isempty(CrossVal)
    CrossVal = 'on';
end
if ~exist('nFolds','var') || isempty(nFolds)
    nFolds = 10;
elseif strcmpi(nFolds, 'leaveout')
    nFolds = nObservations;
end
if ~exist('nCVPpermutations','var') || isempty(nCVPpermutations)
    nCVPpermutations = 1;
end
if ~exist('ClassNames','var') || isempty(ClassNames)
    ClassNames = [1:length(unique(label))];
end
if ~exist('GetConfusionMatrix','var') || isempty(GetConfusionMatrix)
    GetConfusionMatrix = 0;
end
label = strtrim(label);
ClassNames = strtrim(ClassNames);
nClasses = numel(ClassNames);

Kernel = 'linear';
Solver = 'SMO';
if license('test','Distrib_Computing_Toolbox')
     UseParallel = true;
else UseParallel = false;
end


    
    t = templateSVM('Standardize',true,...
        'KernelScale',KS,'BoxConstraint',BC,...
        'KernelFunction',Kernel, 'Solver', Solver,...
        'CrossVal', 'off',...
        'CacheSize','maximal'...
        );

    SVMModel = fitcecoc(cdata,label,...
        'Coding','onevsone',...
        'Learners',t,...
        'ObservationsIn','columns',...
        'ClassNames',ClassNames,...
        'Options',statset('UseParallel',UseParallel)...
        );


    

if strcmp(CrossVal, 'on')

    cvp = cvpartition(label, 'KFold',nFolds);
    
    MisclassRate_eachPerm = nan(nCVPpermutations,1);
    ConfMat_eachPerm = nan(nClasses, nClasses, nCVPpermutations);

    for p = 1 : nCVPpermutations
        
        if p>1
            cvp = repartition(cvp);
        end
        CVMdl = crossval(SVMModel, 'cvpartition', cvp,...
            'Options',statset('UseParallel',UseParallel));

        MisclassRate_eachPerm(p) = kfoldLoss(CVMdl);
        
        if GetConfusionMatrix
            oofLabel = kfoldPredict(CVMdl,'Options',statset('UseParallel',UseParallel));
            ConfMat_eachPerm(:,:,p) = confusionmat(label, oofLabel, 'order',ClassNames);
        end
        
    end
    
    MisclassRate = nanmean(MisclassRate_eachPerm);
    %
    ConfMat_sumAcrossPerms = sum(ConfMat_eachPerm,3);
    ConfMat_perc = ConfMat_sumAcrossPerms ./ repmat(sum(ConfMat_sumAcrossPerms,2),1,nClasses) * 100;
    
else
    MisclassRate = nan;
    ConfMat_perc = nan;
end

