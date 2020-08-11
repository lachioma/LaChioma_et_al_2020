function abovestd = above_stdThreshold_DGD(aDGx, stdThreshold, type, WhichDirs, thr_rep_fraction, ApplySmooth)

% WhichDirs: include only these direction angles (deafult: all).
if nargin < 5 || isempty(thr_rep_fraction)
    thr_rep_fraction = 2/3; % fraction of trials per stim that needs to cross thr_std
end
if nargin < 3 || isempty(type)
        type = 1;
    %     type = 2;
%     type = 3;
end

if nargin < 6 || isempty(ApplySmooth)
    ApplySmooth = 0;
end
if ApplySmooth
    SamplingFreq = aDGx.Info.SamplingFreq;
    smoothingFactor = round(SamplingFreq) +1;
end
                            


if iscell(aDGx.Param)
    frames_prestim = aDGx.Param{1}.frames_prestim ;
    frames_stim    = aDGx.Param{1}.frames_stim ;
else
    frames_prestim = aDGx.Param.frames_prestim ;
    frames_stim    = aDGx.Param.frames_stim ;
end

range1 = frames_prestim+1 : frames_prestim+frames_stim ;

if nargin < 4 || isempty(WhichDirs)
%     WhichDirs_ix = 1:length(angles); % nr.directions
    WhichDirs_ix = 1 : size(aDGx.ROIs(1).dFoF, 2);
else
    if iscell(aDGx.StimSettings)
        angles = aDGx.StimSettings{1}.angles;
    else
        angles = aDGx.StimSettings.angles;
    end
    WhichDirs_ix = find(ismember(angles,WhichDirs));
end

    abovestd = [];
    
    for ix = 1:length(aDGx.ROIs)
        if   isfield(aDGx.ROIs, 'reps' );
            reps = aDGx.ROIs(ix).reps;
        else
            reps = size(aDGx.ROIs(ix).dFoF, 4);
        end
        
        % go through directions and interocular phases to find responsive
        % ROIs:
        ForLoops;
        
    end
    
    disp(['  Responsive ROIs ' inputname(1) ' ' num2str(length(abovestd)) '/' num2str(length(aDGx.ROIs))...
        ' - thr_std=' num2str(stdThreshold) ' thr_type=' num2str(type) ]);

    
    function ForLoops
        for d = 1:size(aDGx.ROIs(ix).dFoF, 2) % directions
            % Skip this direction if is not included in WhichDirs
            if ~ismember(d,WhichDirs_ix)
                continue
            end
            for p = 1:size(aDGx.ROIs(ix).dFoF, 3) % interocular phases
                if isfield(aDGx.ROIs, 'nTrials' )
                    reps = aDGx.ROIs(ix).nTrials(d,p);
                end
                switch type
                    case 1
                        % get 1 max value for each repetition (n_reps values in total):
                        dFoF = aDGx.ROIs(ix).dFoF(range1,d,p,:);
                        if ApplySmooth
                            dFoF = smooth(dFoF, smoothingFactor, 'sgolay',2);
                        end
                        maxima = squeeze( nanmax( dFoF ,[],1));
                    case 2
                        maxima = squeeze(  aDGx.ROIs(ix).max_Fstim_eachTrial(d,p,:) );
                    case 3
                        maxima = squeeze(  aDGx.ROIs(ix).mean_Fstim_eachTrial(d,p,:) );
                end
                if isempty( maxima(~isnan(maxima)) ) % if tmp is all NaN (meaning that there is no trace for that ROI)
                    continue; end;
                maxima = maxima(~isnan(maxima));
                nRepsThreshold = round(thr_rep_fraction*reps);
                if reps == 1, nRepsThreshold = 1; elseif reps==2, nRepsThreshold = 2; end;
                if nnz( maxima > stdThreshold*aDGx.ROIs(ix).std_prestim ) >= nRepsThreshold %max(2,ceil(reps/2))
                    abovestd = [abovestd; ix]; %#ok<AGROW>
                    return
                end
            end
        end
    end
        
end

