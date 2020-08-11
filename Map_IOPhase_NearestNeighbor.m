function [NNexp, NNall] = Map_IOPhase_NearestNeighbor( aDGDap, Angles, thr_type, thr_std, thr_di, SF, Edges, thr_dist, PhasesToExclude, plot_figs)

    save_figures = 0 ;
    
    % Shuffling = 1 for significance of tuning difference vs. cortical distance - exp-wise
    % Shuffling = 2 for significance of tuning difference vs. cortical distance - cell-wise
    Shuffling = 0;  RSiters = 1000;
    
    
ShiftingMethod = 'BasedOnAverageTuningCurve';

if ~exist('SF','var') || isempty(SF)
    error(' input missing: specify SF. Aborting...');
end
if ~exist('Edges','var') || isempty(Edges)
    Edges = [0:20:200, Inf];
end
if ~exist('thr_dist','var') || isempty(thr_dist)
    thr_dist = 0;
end
nBins = length(Edges)-2; % last bin is excluded
if ~exist('PhasesToExclude','var') || isempty(PhasesToExclude)
    PhasesToExclude = NaN;
end
% if isnan(PhasesToExclude)
%     PhasesToInclude = 'all';
% else
%     PhasesToInclude = 0:45:315;
%     PhasesToInclude( ismember(PhasesToInclude, PhasesToExclude) ) = [];
% end
if ~exist('plot_figs','var') || isempty(plot_figs)
    plot_figs = [1:10];
end

IOPhases  = aDGDap.StimSettings{1}.IOPhaseDifferences;
PhaseStep = aDGDap.StimSettings{1}.PhaseStep;
nIOPhases = aDGDap.StimSettings{1}.nPhases;
deltaIOPhases = [0:PhaseStep:180];
type = 'mean';

% NNall.BinEdges = Edges;
ExpIDsall = [aDGDap.ROIs.ExpID];
SpatFreq = str2double(SF);
SFRoiNrs = find([aDGDap.ROIs.SpatFreq] == SpatFreq);
ExpIDs = unique(ExpIDsall(SFRoiNrs));
nExps  = numel(ExpIDs);


% deltaIOPhases = [0:45:180];
% Shift IOPhase distributions to have peak at 180 (5th bin):
% Peak = 5; 
Peak = find(IOPhases == 180); % 7

% mm = 0; % cumulative indexing for pooling all ROIs
% deltaIOPhaseNN_all = [];
% Dist_um_all = [];
% deltaIOPhaseNN_eachExp = cell(nExps, 1);
% Dist_um_eachExp = cell(nExps, 1);
deltaIOPhaseNN_bootstr_binmean_prctile05_eachExp = nan(nExps, nBins);
deltaIOPhaseNN_bootstr_binmean_prctile95_eachExp = nan(nExps, nBins);
deltaIOPhaseNN_binmean_eachExp = [];
Dist_um_byDeltaIOPhase_meanEachExp_mat = nan(nExps, length(deltaIOPhases));
D_1stNearest_perdIOP_meanEachExp = nan(nExps, length(deltaIOPhases));
D_2ndNearest_perdIOP_meanEachExp = nan(nExps, length(deltaIOPhases));
D_1stNearest_perdIOP_meanEachExp_shuffle = nan(nExps, length(deltaIOPhases));
D_2ndNearest_perdIOP_meanEachExp_shuffle = nan(nExps, length(deltaIOPhases));
D_1stNearest_perdIOP_95thEachExp_shuffle = nan(nExps, length(deltaIOPhases));
D_1stNearest_perdIOP_05thEachExp_shuffle = nan(nExps, length(deltaIOPhases));
D_1stNearest_perdIOP_AllCells = []; D_2ndNearest_perdIOP_AllCells = [];
D_1stNearest_perdIOP_AllCells_shuffle = []; D_2ndNearest_perdIOP_AllCells_shuffle = [];

    bin_hist   = [deltaIOPhases, 180+PhaseStep] - PhaseStep/2;
    nBins_dIOP = length(deltaIOPhases);
%     dIOP_allCells_count_eachdIOP  =  cell(nExps, nBins_dIOP);
%     dIOP_1stNearest_eachdIOP      =  cell(nExps, nBins_dIOP);
    D_1stNearest_perdIOP_eachExp          =  cell(nExps, 1);
    dIOP_1stNearest_eachExp       =  cell(nExps, 1);
    D_1stNearest_perPrefIOPhase             =  cell(nExps, nIOPhases, 1);
    dIOP_1stNearest_perPrefIOPhase          =  cell(nExps, nIOPhases);
    dIOP_allOtherCells_allCells             =  cell(nIOPhases, 1);
    dIOP_allOtherCells_meanAllCells         =  nan(nIOPhases, 1);
    dIOP_allOtherCells_stdAllCells          =  nan(nIOPhases, 1);
    dIOP_allOtherCells_semAllCells          =  nan(nIOPhases, 1);
    dIOP_allOtherCells_count_perPrefIOPhase =  nan(nExps, nIOPhases, nBins_dIOP);
    dIOP_allOtherCells_prob_perPrefIOPhase  =  nan(nExps, nIOPhases, nBins_dIOP);
    dIOP_1stNearest_perPrefIOPhase_prob     =  nan(nExps, nIOPhases, nBins_dIOP);
    dIOP_1stNearest_perPrefIOPhase_allCells     = cell(nIOPhases, 1);
    dIOP_1stNearest_perPrefIOPhase_meanAllCells =  nan(nIOPhases, 1);
    dIOP_1stNearest_perPrefIOPhase_stdAllCells  =  nan(nIOPhases, 1);
    dIOP_1stNearest_perPrefIOPhase_semAllCells  =  nan(nIOPhases, 1);
    

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
    Info = aDGDap.Info{e};
    exptag = Info.Exp(1:5);
    % Prepare variable aDGDtemp:
    ExpIDRoiNrs = find([aDGDap.ROIs.ExpID] == expid);
    SFExpIDRoiNrs = intersect(SFRoiNrs, ExpIDRoiNrs);
    aDGDtmp.ROIs = aDGDap.ROIs(SFExpIDRoiNrs);
    % I need this because GetDI_vs_Ori needs it:
    aDGDtmp.StimSettings = aDGDap.StimSettings(e); %it's a cell array
    aDGDtmp.Param = aDGDap.Param(e);
    
    Ori = [];
    i = 1;
    Ori = GetDI_vs_Ori( Ori, Angles, SF, aDGDtmp, thr_type, thr_std, i);

    ts = 1;
    OSFrespRoiNrs = Ori(i).OSFrespRoiNrs{ts};
    DI = Ori(i).StatsDI{ts}.DI;
    DIOSFrespRoiNrs = OSFrespRoiNrs(DI > thr_di);
    nCells = length(DIOSFrespRoiNrs);
%     PhaseStep = aDGDap(e).StimSettings.PhaseStep;
%     nIOPhases = aDGDap(e).StimSettings.nPhases;
    
    angle_num = Angles(1); if size(Angles,1)==2, angle_num = Angles(1,1)+Angles(2,1)/100; end;
    fprintf(' - Nr.responsive ROIs with thr_di=%3.2f,Ori=%3.2f -> %3.0f\n', thr_di,angle_num,nCells);

    for m = 1 : nCells
        nr = DIOSFrespRoiNrs(m);
        DirMaxIx = aDGDtmp.ROIs(nr).DirMaxIx;
        if strcmp(type,'mean')
            PhaseMax  = aDGDtmp.ROIs(nr).PhaseMax_mean;
        elseif strcmp(type,'max')
            PhaseMax  = aDGDtmp.ROIs(nr).PhaseMax_max;
        end
        PrefIOPhase(m) = PhaseMax(DirMaxIx); %#ok<*AGROW>
        PrefIOPhaseIx(m) = PrefIOPhase(m)/PhaseStep +1;
    end
    
% % %     CountPrefIOPhase = histcounts(PrefIOPhase,[0-PhaseStep/2:PhaseStep:360-PhaseStep+PhaseStep/2]);
% % %     % If there is no single peak:
% % %     if sum(CountPrefIOPhase==max(CountPrefIOPhase)) > 1
% % % %         CountPrefIOPhase_sm = smooth(CountPrefIOPhase,3);
% % %         CountPrefIOPhase_interp = interp(CountPrefIOPhase,360/8,3);
% % %         CountPrefIOPhase_interpext = [CountPrefIOPhase_interp(end-99:end),CountPrefIOPhase_interp,CountPrefIOPhase_interp(1:100)];
% % %         CountPrefIOPhase_sm = smooth(CountPrefIOPhase_interpext, 90);
% % %         CountPrefIOPhase_sm = CountPrefIOPhase_sm(100+[1:360]);
% % %         [~,peak_distrib_deg] = max(CountPrefIOPhase_sm);
% % %         [~,~,peak_distrib_ix] = histcounts(peak_distrib_deg,[0-PhaseStep/2:PhaseStep:360-PhaseStep+PhaseStep/2]);
% % %     else % otherwise take directly the peak:
% % %         [~,peak_distrib_ix] = max(CountPrefIOPhase);
% % %     end
% % %     shift_peak = Peak - peak_distrib_ix; % positive value -> shift to right
% % % %     CountPrefIOPhase_shift(i,:) = circshift(CountPrefIOPhase, [0,shift_peak]);
    
    if ~isnan(PhasesToExclude)
        % in this case there is no need to shift the tuning curves
        switch ShiftingMethod

            case 'BasedOnAverageTuningCurve'

                TuningCurve_AllCells = nan(nCells(i),nIOPhases);
                for m = 1 : nCells(i)
                    nr = DIOSFrespRoiNrs(m);
                        if strcmp(type,'mean')
                            TuningCurve = aDGDtmp.ROIs(nr).DisparityCurve_mean;
                        elseif strcmp(type,'max')
                            TuningCurve = aDGDtmp.ROIs(nr).DisparityCurve_max;
                        end
                        % Normalize curves:
                        TuningCurve(TuningCurve<0) = 0;
                        TuningCurve = TuningCurve / max(TuningCurve);
                        TuningCurve_AllCells(m,:) = TuningCurve;
                end
                TuningCurve_AllCells_mean = nanmean(TuningCurve_AllCells, 1);
                [~, peak_tuning_ix] = max(TuningCurve_AllCells_mean);
                shift_peak(i) = Peak - peak_tuning_ix; % positive value -> shift to right


            case 'BasedOnIOPhaseDistribution'

                CountPrefIOPhase = histcounts(PrefIOPhase,[0-PhaseStep/2:PhaseStep:360-PhaseStep+PhaseStep/2]);
                % If there is no single peak:
                if sum(CountPrefIOPhase==max(CountPrefIOPhase)) > 1
        %         CountPrefIOPhase_sm = smooth(CountPrefIOPhase,3);
                    CountPrefIOPhase_interp = interp(CountPrefIOPhase,360/8,3);
                    CountPrefIOPhase_interpext = [CountPrefIOPhase_interp(end-99:end),CountPrefIOPhase_interp,CountPrefIOPhase_interp(1:100)];
                    CountPrefIOPhase_sm = smooth(CountPrefIOPhase_interpext, 90);
                    CountPrefIOPhase_sm = CountPrefIOPhase_sm(100+[1:360]);
                    [~,peak_distrib_deg] = max(CountPrefIOPhase_sm);
                    [~,~,peak_distrib_ix] = histcounts(peak_distrib_deg,[0-PhaseStep/2:PhaseStep:360-PhaseStep+PhaseStep/2]);
                else % otherwise take directly the peak:
                    [~,peak_distrib_ix] = max(CountPrefIOPhase);
                end
                shift_peak(i) = Peak - peak_distrib_ix; % positive value -> shift to right
        end
    else
        shift_peak = 0;
    end
    
    PrefIOPhase_shift = PrefIOPhase + PhaseStep*shift_peak;
    PrefIOPhase_shift( PrefIOPhase_shift >= 360) = PrefIOPhase_shift( PrefIOPhase_shift >= 360) - 360;
    PrefIOPhase_shift( PrefIOPhase_shift <    0) = PrefIOPhase_shift( PrefIOPhase_shift <    0) + 360;
    
    mi = 0;
    for m = 1 : nCells
        if ~ismember(PrefIOPhase_shift(m), PhasesToExclude)
            nr = DIOSFrespRoiNrs(m);
            mi = mi + 1;
            DIOSFrespRoiNrs_excl(mi) = nr;
            PrefIOPhase_excl(mi) = PrefIOPhase_shift(m);
        end
    end
    
    nCells_excl = length(PrefIOPhase_excl);
    
    if nCells_excl < 2
        disp('   -- exp skipped, less than 2 cells are avaialble')
        continue
    end
    
    % Get pixel size:
%     save_dir  = ['I:\Alessandro La Chioma\AnalyzedData\Alessandro\' Info.Mouse '\' Info.Spot '\' Info.Date '\'];
%     FP = load([save_dir 'FieldParameters-' exptag '.mat']);
    FP = aDGDap.FP{e};
    umppx_x = FP.px_size_x;
    umppx_y = FP.px_size_y;
    % Load ROI variable:
%     load([save_dir 'ROI-' exptag '.mat']);
    ROI = aDGDap.ROI{e};
    X = [ROI(DIOSFrespRoiNrs_excl).x]'; %you need square brackets here and a column vector.
    Y = [ROI(DIOSFrespRoiNrs_excl).y]';
    X_um = X * umppx_x;
    Y_um = Y * umppx_y;
    Dall_um = pdist([X_um,Y_um]); % distances are arranged in the order (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2), ..., (m,m–1))
    % Make D square matrix where i,j is distance between i-th and j-th of XY:
    Dall_um = squareform(Dall_um);
    % Dall_px = pdist([X   ,Y]); 
    % Dall_px = squareform(Dall_px);
    
    
    D = nan; D_thr = nan;
    deltaIOPhaseNN = nan;
    deltaIOPhaseNN_binmean = nan(nBins, 1);
    deltaIOPhaseNN_binsem  = nan(nBins, 1);
    
%     nNeighbors = nCells_excl-1;
    if nCells_excl > 1
        
        CellPairs_ix = combnk(1:nCells_excl, 2);
        nCellPairs = size(CellPairs_ix, 1);
        
        D                = nan(nCellPairs,1);
        deltaIOPhaseNN   = nan(nCellPairs,1);
        
        for nc = 1 : nCellPairs
            
            ix1 = CellPairs_ix(nc,1);
            ix2 = CellPairs_ix(nc,2);
%             nr1 = DIOSFrespRoiNrs_excl( ix1 );
%             nr2 = DIOSFrespRoiNrs_excl( ix2 );
            
            D(nc) = Dall_um(ix1,ix2);
            
            deltaIOPhaseNN(nc) = round( circ_dist(PrefIOPhase_excl(ix1)*pi/180, PrefIOPhase_excl(ix2)*pi/180) *180/pi );

        end
%         deltaIOPhaseNN = int32(deltaIOPhaseNN);


        D_ix_toInclude = find( D > thr_dist );
        D_thr = D( D_ix_toInclude );
        deltaIOPhaseNN_thr = deltaIOPhaseNN( D_ix_toInclude );
%         Dist_um_eachExp{e} = D;
%         deltaIOPhaseNN_eachExp{e} = deltaIOPhaseNN;
        

        
        [~,~,bins] = histcounts(D_thr, Edges);
        for b = 1:nBins
            binmean = circ_mean((abs(deltaIOPhaseNN_thr(bins==b))*pi/180))*180/pi ;
            [~, s0] =  circ_std( abs(deltaIOPhaseNN_thr(bins==b))*pi/180 ) ;
            binsem = s0*180/pi / sqrt(sum(bins==b));
            if isempty(binmean), binmean = NaN; binsem = NaN;end;
            deltaIOPhaseNN_binmean(b) = binmean;
            deltaIOPhaseNN_binsem(b)  = binsem;
        end
        
        
        if ismember(1, Shuffling )
            deltaIOPhaseNN_bootstr_binmean = nan(RSiters, nBins);
            for rsi = 1 : RSiters
                D_bootstr = datasample(D_thr, length(D_thr));

                [~,~,bins] = histcounts(D_bootstr, Edges);
                for b = 1 : nBins
                    binmean = circ_mean((abs(deltaIOPhaseNN_thr(bins==b))*pi/180))*180/pi ;
%                     [~, s0] =  circ_std( abs(deltaIOPhaseNN(bins==b))*pi/180 ) ;
%                     binsem = s0*180/pi / sqrt(sum(bins==b));
                    if isempty(binmean), binmean = NaN; binsem = NaN;end;
                    deltaIOPhaseNN_bootstr_binmean(rsi,b) = binmean;
                end
            end

            deltaIOPhaseNN_bootstr_binmean_prctile05_eachExp(e,:) = prctile(deltaIOPhaseNN_bootstr_binmean,  2.5, 1);
            deltaIOPhaseNN_bootstr_binmean_prctile95_eachExp(e,:) = prctile(deltaIOPhaseNN_bootstr_binmean, 97.5, 1);
        else
            deltaIOPhaseNN_bootstr_binmean_prctile05_eachExp(e,:) = nan(1, nBins);
            deltaIOPhaseNN_bootstr_binmean_prctile95_eachExp(e,:) = nan(1, nBins);
        end
        
        
        
        D_1stNearest_perdIOP = nan(nCells_excl, numel(deltaIOPhases));
        D_2ndNearest_perdIOP = nan(nCells_excl, numel(deltaIOPhases));
        
        for c = 1 : nCells_excl
            for dp = 1 : numel(deltaIOPhases)
                dIOP = deltaIOPhases(dp);
%                 Ix = find(PrefIOPhase_excl == mod(PrefIOPhase_excl(c) + dIOP, 360));
%                 Ix = [ find( PrefIOPhase_excl == mod(PrefIOPhase_excl(c) + dIOP, 360)) , ...
%                        find( PrefIOPhase_excl == mod(PrefIOPhase_excl(c) - dIOP, 360)) ];
%                 Ix = unique(Ix); % remove double cells, you need it for when dIOP==0 or dIOP==180
                Ix = [ find( PrefIOPhase_excl == mod(PrefIOPhase_excl(c) + dIOP, 360)  | ...) , ...
                             PrefIOPhase_excl == mod(PrefIOPhase_excl(c) - dIOP, 360)) ];
                Ix(Ix==c) = [];  % remove the same cell, you need it for when dIOP==0
                if isempty(Ix)
                    D_1stNearest_perdIOP(c,dp) = nan;
                    D_2ndNearest_perdIOP(c,dp) = nan;
                else
                    D_forEachIx_sorted = sort(Dall_um(c,Ix),'ascend');
                    D_1stNearest_perdIOP(c,dp) = D_forEachIx_sorted(1);
                    if length(Ix) > 1
                        D_2ndNearest_perdIOP(c,dp) = D_forEachIx_sorted(2);
                    else
                        D_2ndNearest_perdIOP(c,dp) = nan;
                    end
                end
            end
        end

    
    
        deltaIOP_mat = pdist([PrefIOPhase_excl]'); % distances are arranged in the order (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2), ..., (m,m–1))
        % Make deltaIOP_mat square matrix where i,j is distance between i-th and j-th of XY:
        deltaIOP_mat = squareform(deltaIOP_mat);
        deltaIOP_mat = abs(angwrap(deltaIOP_mat));

        D_1stNearest_perdIOP = nan(nCells_excl, numel(deltaIOPhases));
        D_2ndNearest_perdIOP = nan(nCells_excl, numel(deltaIOPhases));

        for dp = 1 : numel(deltaIOPhases)
            dIOP = deltaIOPhases(dp);

            CellPairs_ix_withThisdIOP = find( abs(deltaIOPhaseNN) == dIOP );

            CellPairs_ix_list_Nearest1st_alreadyUsed = [];
            CellPairs_ix_list_Nearest2nd_alreadyUsed = [];


            for c = 1 : nCells_excl

                CellPairs_ix_withThisCell = union( find( CellPairs_ix(:,1) == c ), find( CellPairs_ix(:,2) == c ) );
                CellPairs_ix_withThisdIOP_withThisCell = intersect(CellPairs_ix_withThisdIOP, CellPairs_ix_withThisCell);
%                 CellPairs_ix_withThisdIOP_withThisCell_notUsedYet = intersect(CellPairs_ix_withThisdIOP_withThisCell, setxor(CellPairs_ix_withThisdIOP_withThisCell, CellPairs_ix_list_Nearest1st_alreadyUsed) );

                if isempty(CellPairs_ix_withThisdIOP_withThisCell)
                    D_1stNearest_perdIOP(c,dp) = nan;
                    D_2ndNearest_perdIOP(c,dp) = nan;
                else

                    Distances_fromThisCell = D(CellPairs_ix_withThisdIOP_withThisCell);
                    Distances_fromThisCell_sorted = sort(Distances_fromThisCell, 'ascend');

                    Dist_Nearest1st_fromThisCell = Distances_fromThisCell_sorted(1);
                    Ix1 = find( Distances_fromThisCell == Dist_Nearest1st_fromThisCell );
                    CellPairs_ix_Nearest1st = CellPairs_ix_withThisdIOP_withThisCell(Ix1);

                    if ismember(CellPairs_ix_Nearest1st, CellPairs_ix_list_Nearest1st_alreadyUsed);
                        continue
                    else
                        CellPairs_ix_list_Nearest1st_alreadyUsed = [CellPairs_ix_list_Nearest1st_alreadyUsed; CellPairs_ix_Nearest1st];
                        D_1stNearest_perdIOP(c,dp) = Dist_Nearest1st_fromThisCell;
                    end

                    if length(CellPairs_ix_withThisdIOP_withThisCell) > 1
                        Dist_Nearest2nd_fromThisCell = Distances_fromThisCell_sorted(2);
                        Ix2 = find( Distances_fromThisCell == Dist_Nearest2nd_fromThisCell );
                        CellPairs_ix_Nearest2nd = CellPairs_ix_withThisdIOP_withThisCell(Ix2);
                        if ismember(CellPairs_ix_Nearest2nd, CellPairs_ix_list_Nearest2nd_alreadyUsed);
                            continue
                        else
                            CellPairs_ix_list_Nearest2nd_alreadyUsed = [CellPairs_ix_list_Nearest2nd_alreadyUsed; CellPairs_ix_withThisdIOP_withThisCell(Ix2)];
                            D_2ndNearest_perdIOP(c,dp) = Dist_Nearest2nd_fromThisCell;
                        end
                    else
                        D_2ndNearest_perdIOP(c,dp) = nan;
                    end
                end
            end
        end

        
        
        
        
        
        
        
        
        
        %% Shuffle test
        nShuffles = 1000;
        
        D_1stNearest_perdIOP_shuffle = nan(nCells_excl, numel(deltaIOPhases), nShuffles);
        D_2ndNearest_perdIOP_shuffle = nan(nCells_excl, numel(deltaIOPhases), nShuffles);
        D_1stNearest_perdIOP_meanShuffle = nan(nCells_excl, numel(deltaIOPhases));
        D_2ndNearest_perdIOP_meanShuffle = nan(nCells_excl, numel(deltaIOPhases));
        
        
        for sh = 1 : nShuffles
        
        for dp = 1 : numel(deltaIOPhases)
            dIOP = deltaIOPhases(dp);

            CellPairs_ix_withThisdIOP = find( abs(deltaIOPhaseNN) == dIOP );

            CellPairs_ix_list_Nearest1st_alreadyUsed = [];
            CellPairs_ix_list_Nearest2nd_alreadyUsed = [];


            for c = 1 : nCells_excl

                CellPairs_ix_withThisCell = union( find( CellPairs_ix(:,1) == c ), find( CellPairs_ix(:,2) == c ) );
                CellPairs_ix_withThisdIOP_withThisCell = intersect(CellPairs_ix_withThisdIOP, CellPairs_ix_withThisCell);
%                 CellPairs_ix_withThisdIOP_withThisCell_notUsedYet = intersect(CellPairs_ix_withThisdIOP_withThisCell, setxor(CellPairs_ix_withThisdIOP_withThisCell, CellPairs_ix_list_Nearest1st_alreadyUsed) );

                if isempty(CellPairs_ix_withThisdIOP_withThisCell)
                    D_1stNearest_perdIOP_meanShuffle(c,dp) = nan;
                    D_2ndNearest_perdIOP_meanShuffle(c,dp) = nan;
                else
                    
%                     for sh = 1 : nShuffles
                        
                        D_shuffle = D( randperm(length(D)) );
                        
                        Distances_fromThisCell = D_shuffle(CellPairs_ix_withThisdIOP_withThisCell);
                        Distances_fromThisCell_sorted = sort(Distances_fromThisCell, 'ascend');

                        Dist_Nearest1st_fromThisCell = Distances_fromThisCell_sorted(1);
                        Ix1 = find( Distances_fromThisCell == Dist_Nearest1st_fromThisCell );
                        CellPairs_ix_Nearest1st = CellPairs_ix_withThisdIOP_withThisCell(Ix1);

                        if ismember(CellPairs_ix_Nearest1st, CellPairs_ix_list_Nearest1st_alreadyUsed);
                            continue
                        else
                            CellPairs_ix_list_Nearest1st_alreadyUsed = [CellPairs_ix_list_Nearest1st_alreadyUsed; CellPairs_ix_Nearest1st];
                            D_1stNearest_perdIOP_shuffle(c,dp,sh) = Dist_Nearest1st_fromThisCell;
                        end

                        if length(CellPairs_ix_withThisdIOP_withThisCell) > 1
                            Dist_Nearest2nd_fromThisCell = Distances_fromThisCell_sorted(2);
                            Ix2 = find( Distances_fromThisCell == Dist_Nearest2nd_fromThisCell );
                            CellPairs_ix_Nearest2nd = CellPairs_ix_withThisdIOP_withThisCell(Ix2);
                            if ismember(CellPairs_ix_Nearest2nd, CellPairs_ix_list_Nearest2nd_alreadyUsed);
                                continue
                            else
                                CellPairs_ix_list_Nearest2nd_alreadyUsed = [CellPairs_ix_list_Nearest2nd_alreadyUsed; CellPairs_ix_withThisdIOP_withThisCell(Ix2)];
                                D_2ndNearest_perdIOP_shuffle(c,dp,sh) = Dist_Nearest2nd_fromThisCell;
                            end
                        else
                            D_2ndNearest_perdIOP_shuffle(c,dp,sh) = nan;
                        end
%                     end
                end
            end
        end
        
        end
    
        D_1stNearest_perdIOP_meanShuffle = nanmean(D_1stNearest_perdIOP_shuffle, 3);
        D_2ndNearest_perdIOP_meanShuffle = nanmean(D_2ndNearest_perdIOP_shuffle, 3);
        [s1, s2, s3] = size( D_1stNearest_perdIOP_shuffle );
        D_1stNearest_perdIOP_shuffle_resh = reshape( permute( D_1stNearest_perdIOP_shuffle, [1 3 2] ), s1*s3, s2, 1 );
        D_2ndNearest_perdIOP_shuffle_resh = reshape( permute( D_2ndNearest_perdIOP_shuffle, [1 3 2] ), s1*s3, s2, 1 );

%%
        
        
  
    
%     else
%         D_1stNearest_perdIOP = nan(1, numel(deltaIOPhases));
%         D_2ndNearest_perdIOP = nan(1, numel(deltaIOPhases));
%     end
    
    
%     cc = 0;
%     Ix_list = [];
%     Ix_Nearest1st_list = [];
%     PrefIOPhase_diff_Nearest1st = [];
%     
%     for c = 1 : nCells_excl
%         Distances_fromThisCell        = Dall_um(:,c);
%         Distances_fromThisCell_sorted = sort(Distances_fromThisCell, 'ascend');
%         Dist_Nearest1st_fromThisCell  = Distances_fromThisCell_sorted(2);
%         Ix_Nearest1st_fromThisCell    = find( Distances_fromThisCell == Dist_Nearest1st_fromThisCell );
%         if ismember(c, Ix_Nearest1st_list);
%             continue
%         else
%             cc = cc + 1;
%             Ix_list(cc,1) = c; % useful for debugging
%             Ix_Nearest1st_list = [Ix_Nearest1st_list; Ix_Nearest1st_fromThisCell];
%             
%             PrefIOPhase_diff_Nearest1st(cc,1) = round( abs( circ_dist( PrefIOPhase_excl(c)*pi/180, PrefIOPhase_excl(Ix_Nearest1st_fromThisCell)*pi/180 ) *180/pi ) );
%         end
%     end
%     
%     bin_hist = [deltaIOPhases, 180+PhaseStep] - PhaseStep/2;
% %     PrefIOPhase_diff_Nearest1st_count = histcounts(PrefIOPhase_diff_Nearest1st, bin_hist);
%     PrefIOPhase_diff_Nearest1st_count = histcounts(PrefIOPhase_diff_Nearest1st, bin_hist, 'Normalization','probability');
%     
%     NNexp(e).PrefIOPhase_diff_Nearest1st_eachExp = PrefIOPhase_diff_Nearest1st;
%     NNexp(e).PrefIOPhase_diff_Nearest1st_count_eachExp = PrefIOPhase_diff_Nearest1st_count;
%     
%     figure;
%     histogram(PrefIOPhase_diff_Nearest1st, bin_hist);
%     
%     figure;
%     histogram(abs(deltaIOPhaseNN_thr), bin_hist);
    
    


    
        
    D_1stNearest                   =   nan(nCells_excl, 1);
    dIOP_1stNearest                =   nan(nCells_excl, 1);
%     dp_ix                          = zeros(nBins_dIOP, 1);
    dIOP_allOtherCells_count_perPrefIOPhase_thisExp = cell( nIOPhases, 1);
    dIOP_allOtherCells_prob_perPrefIOPhase_thisExp  = cell( nIOPhases, 1);
    dIOP_allOtherCells_count       =   nan(nCells_excl, nBins_dIOP);
    dIOP_allOtherCells_prob        =   nan(nCells_excl, nBins_dIOP);
    
    for c = 1 : nCells_excl
        % get index of IOPhase (from 1 to 12 (or nPhases)):
        px = 1 + PrefIOPhase_excl(c) / PhaseStep ;

        [D_forEachIx_sorted, D_Ix] = sort(Dall_um(c,:),'ascend');

        D_1stNearest(c) = D_forEachIx_sorted(2);
        D_1stNearest_perPrefIOPhase{e,px} = [D_1stNearest_perPrefIOPhase{e,px}; D_1stNearest(c)];
        
        dIOP_1stNearest(c) = PrefIOPhase_excl(c) - PrefIOPhase_excl(D_Ix(2));
        dIOP_1stNearest(c) = abs( angwrap( dIOP_1stNearest(c) ) );
        dIOP_1stNearest_perPrefIOPhase{e,px} = [dIOP_1stNearest_perPrefIOPhase{e,px}; dIOP_1stNearest(c)];
        
        PrefIOPhase_excl_minusThisCell = PrefIOPhase_excl;
        PrefIOPhase_excl_minusThisCell(c) = [];

        dIOP_allOtherCells = PrefIOPhase_excl(c) - PrefIOPhase_excl_minusThisCell;
        dIOP_allOtherCells = abs( angwrap(dIOP_allOtherCells) );
        dIOP_allOtherCells = reshape(dIOP_allOtherCells, [],1);
        dIOP_allOtherCells_allCells{px} = [dIOP_allOtherCells_allCells{px}; dIOP_allOtherCells];
        
        dIOP_allOtherCells_count(c,:) = histcounts(dIOP_allOtherCells, bin_hist,...
            'Normalization','count');
        dIOP_allOtherCells_prob(c,:)  = histcounts(dIOP_allOtherCells, bin_hist,...
            'Normalization','probability');
        dIOP_allOtherCells_count_perPrefIOPhase_thisExp{px} = [ dIOP_allOtherCells_count_perPrefIOPhase_thisExp{px} ; dIOP_allOtherCells_count(c,:) ];
        dIOP_allOtherCells_prob_perPrefIOPhase_thisExp{px}  = [ dIOP_allOtherCells_prob_perPrefIOPhase_thisExp{px}  ; dIOP_allOtherCells_prob(c,:) ];

%         dp = 1+ dIOP_1stNearest(c) / PhaseStep ;
%         dp_ix(dp) = dp_ix(dp) + 1;
%         dIOP_allCells_count_eachdIOP{e,dp}(dp_ix(dp),:) = dIOP_allOtherCells_count;
%         dIOP_1stNearest_eachdIOP{e,dp}(dp_ix(dp),1) = dIOP_1stNearest(c);
    end
    
    
    
    else
        D_1stNearest_perdIOP = nan(1, numel(deltaIOPhases));
        D_2ndNearest_perdIOP = nan(1, numel(deltaIOPhases));
    end
    
    
    
    
    D_1stNearest_eachExp{e}    = D_1stNearest;
    dIOP_1stNearest_eachExp{e} = dIOP_1stNearest;

    for px = 1 : nIOPhases
        if isempty( dIOP_1stNearest_perPrefIOPhase{e,px} )
            dIOP_1stNearest_perPrefIOPhase{e,px} = nan;
            dIOP_allOtherCells_prob_perPrefIOPhase_thisExp{px} = nan(1,nBins_dIOP);
        end
        
        dIOP_1stNearest_perPrefIOPhase_prob(e,px,:) = histcounts( dIOP_1stNearest_perPrefIOPhase{e,px}, bin_hist,...
            'Normalization','probability');
        
        dIOP_allOtherCells_prob_perPrefIOPhase(e,px,:) = mean(dIOP_allOtherCells_prob_perPrefIOPhase_thisExp{px}, 1);
    end










    
    NNexp(e).deltaIOPhaseNN_binMean_eachExp = deltaIOPhaseNN_binmean;
    NNexp(e).deltaIOPhaseNN_binSEM_eachExp  = deltaIOPhaseNN_binsem;
    NNexp(e).Dist_um_eachExp            = D;
    NNexp(e).deltaIOPhaseNN_eachExp     = deltaIOPhaseNN;
    NNexp(e).Dist_um_thr_eachExp        = D_thr;
    NNexp(e).deltaIOPhaseNN_thr_eachExp = deltaIOPhaseNN_thr;
    
    
    deltaIOPhaseNN_binmean_eachExp = [deltaIOPhaseNN_binmean_eachExp; deltaIOPhaseNN_binmean'];
    
    
    Dist_um_byDeltaIOPhase = cell(numel(deltaIOPhases), 1);
    for dp = 1 : numel(deltaIOPhases)
        dIOP = deltaIOPhases(dp);
        Dist_um_byDeltaIOPhase{dp} = D_thr( abs(int32(deltaIOPhaseNN_thr))==dIOP );
        Dist_um_byDeltaIOPhase_meanEachExp(dp) = mean(Dist_um_byDeltaIOPhase{dp});
        Dist_um_byDeltaIOPhase_stdEachExp(dp)  =  std(Dist_um_byDeltaIOPhase{dp});
        Dist_um_byDeltaIOPhase_semEachExp(dp)  = Dist_um_byDeltaIOPhase_stdEachExp(dp)/sqrt(numel(Dist_um_byDeltaIOPhase{dp}));
    end

    NNexp(e).Dist_um_byDeltaIOPhase_meanEachExp = Dist_um_byDeltaIOPhase_meanEachExp;
    NNexp(e).Dist_um_byDeltaIOPhase_stdEachExp  = Dist_um_byDeltaIOPhase_stdEachExp;
    NNexp(e).Dist_um_byDeltaIOPhase_semEachExp  = Dist_um_byDeltaIOPhase_semEachExp;
    Dist_um_byDeltaIOPhase_meanEachExp_mat(e,:) = Dist_um_byDeltaIOPhase_meanEachExp;

    
    D_1stNearest_perdIOP_meanEachExp(e,:) = nanmean(D_1stNearest_perdIOP, 1);
    D_2ndNearest_perdIOP_meanEachExp(e,:) = nanmean(D_2ndNearest_perdIOP, 1);
    D_1stNearest_perdIOP_AllCells = [D_1stNearest_perdIOP_AllCells; D_1stNearest_perdIOP];
    D_2ndNearest_perdIOP_AllCells = [D_2ndNearest_perdIOP_AllCells; D_2ndNearest_perdIOP];
    % shuffles
    D_1stNearest_perdIOP_meanEachExp_shuffle(e,:) = nanmean(D_1stNearest_perdIOP_meanShuffle, 1);
    D_2ndNearest_perdIOP_meanEachExp_shuffle(e,:) = nanmean(D_2ndNearest_perdIOP_meanShuffle, 1);
    D_1stNearest_perdIOP_AllCells_shuffle = [D_1stNearest_perdIOP_AllCells_shuffle; D_1stNearest_perdIOP_shuffle_resh];
    D_2ndNearest_perdIOP_AllCells_shuffle = [D_2ndNearest_perdIOP_AllCells_shuffle; D_2ndNearest_perdIOP_shuffle_resh];
%     D_1stNearest_perdIOP_95thEachExp_shuffle(e,:) = prctile(D_1stNearest_perdIOP_meanShuffle, 95, 1);
%     D_1stNearest_perdIOP_05thEachExp_shuffle(e,:) = prctile(D_1stNearest_perdIOP_meanShuffle,  5, 1);
    D_1stNearest_perdIOP_95thEachExp_shuffle(e,:) = nanmean( prctile(D_1stNearest_perdIOP_shuffle, 97.5, 3), 1 );
    D_1stNearest_perdIOP_05thEachExp_shuffle(e,:) = nanmean( prctile(D_2ndNearest_perdIOP_shuffle,  2.5, 3), 1 );
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for px = 1 : nIOPhases
    dIOP_1stNearest_perPrefIOPhase_prob_mean(px,:) = nanmean(dIOP_1stNearest_perPrefIOPhase_prob(:,px,:), 1);
    dIOP_1stNearest_perPrefIOPhase_prob_std(px,:)  =  nanstd(dIOP_1stNearest_perPrefIOPhase_prob(:,px,:), 0,1);
    samplesize_perbin = reshape( sum( ~isnan(dIOP_1stNearest_perPrefIOPhase_prob(:,px,:)), 1), 1,[],1 );
    dIOP_1stNearest_perPrefIOPhase_prob_sem(px,:)  =  dIOP_1stNearest_perPrefIOPhase_prob_std(px,:) ./ sqrt(samplesize_perbin);
    
    dIOP_1stNearest_perPrefIOPhase_allCells{px} = cell2mat(dIOP_1stNearest_perPrefIOPhase(:,px));
    dIOP_1stNearest_perPrefIOPhase_meanAllCells(px) = nanmean(dIOP_1stNearest_perPrefIOPhase_allCells{px}, 1);
    dIOP_1stNearest_perPrefIOPhase_stdAllCells(px) =   nanstd(dIOP_1stNearest_perPrefIOPhase_allCells{px}, 0,1);
    samplesize = sum(~isnan(dIOP_1stNearest_perPrefIOPhase_allCells{px}));
    dIOP_1stNearest_perPrefIOPhase_semAllCells(px) = dIOP_1stNearest_perPrefIOPhase_stdAllCells(px) / sqrt(samplesize);
    
    dIOP_allOtherCells_prob_perPrefIOPhase_mean(px,:) = nanmean(dIOP_allOtherCells_prob_perPrefIOPhase(:,px,:), 1);
    dIOP_allOtherCells_prob_perPrefIOPhase_std(px,:)  =  nanstd(dIOP_allOtherCells_prob_perPrefIOPhase(:,px,:), 0,1);
    samplesize_perbin = reshape( sum( ~isnan(dIOP_allOtherCells_prob_perPrefIOPhase(:,px,:)), 1), 1,[],1 );
    dIOP_allOtherCells_prob_perPrefIOPhase_sem(px,:)  =  dIOP_allOtherCells_prob_perPrefIOPhase_std(px,:) ./ sqrt(samplesize_perbin);
    
    dIOP_allOtherCells_meanAllCells(px) = nanmean(dIOP_allOtherCells_allCells{px});
    dIOP_allOtherCells_stdAllCells(px)  =  nanstd(dIOP_allOtherCells_allCells{px});
    samplesize = sum(~isnan(dIOP_allOtherCells_allCells{px}));
    dIOP_allOtherCells_semAllCells(px)  = dIOP_allOtherCells_stdAllCells(px)/sqrt(samplesize);
end
% figure;
% for px = 1 : nIOPhases
%     subaxis(2, nIOPhases/2, px);
%     hold on
%     plot(dIOP_1stNearest_perPrefIOPhase_prob_mean(px,:));
%     AddErrorArea([],dIOP_1stNearest_perPrefIOPhase_prob_mean(px,:), dIOP_1stNearest_perPrefIOPhase_prob_sem(px,:));
%     %
%     plot(dIOP_allOtherCells_prob_perPrefIOPhase_mean(px,:));
%     AddErrorArea([],dIOP_allOtherCells_prob_perPrefIOPhase_mean(px,:), dIOP_allOtherCells_prob_perPrefIOPhase_sem(px,:));
% end

%% Average distance of 1st neighbor
% for px = 1 : nIOPhases
%     D_1stNearest_perPrefIOPhase_mat{px} = cell2mat( D_1stNearest_perPrefIOPhase(:,px) );
%     D_1stNearest_perPrefIOPhase_mat_mean(px,1) = nanmean(D_1stNearest_perPrefIOPhase_mat{px});
%     D_1stNearest_perPrefIOPhase_mat_std(px,1)  =  nanstd(D_1stNearest_perPrefIOPhase_mat{px});
%     samplesize = sum(~isnan(D_1stNearest_perPrefIOPhase_mat{px}));
%     D_1stNearest_perPrefIOPhase_mat_sem(px,1)  =  D_1stNearest_perPrefIOPhase_mat_std(px,1) / sqrt(samplesize);
% end
%%

%% Mean tuning difference vs. cortical distance - cell-wise
Dist_um_allCells        = vertcat(NNexp.Dist_um_thr_eachExp);
deltaIOPhaseNN_allCells = vertcat(NNexp.deltaIOPhaseNN_thr_eachExp);

deltaIOPhaseNN_allCells_binmean =  nan(nBins, 1);
deltaIOPhaseNN_allCells_binsem  =  nan(nBins, 1);
deltaIOPhaseNN_allCells_eachBin = cell(nBins, 1);
deltaIOPhaseNN_allCells_eachBin_Group = [];

[~,~,bins] = histcounts(Dist_um_allCells, Edges);
for b = 1 : nBins
    binmean = circ_mean((abs(deltaIOPhaseNN_allCells(bins==b))*pi/180))*180/pi ;
    [~, s0] =  circ_std( abs(deltaIOPhaseNN_allCells(bins==b))*pi/180 ) ;
    binsem = s0*180/pi / sqrt(sum(bins==b));
    if isempty(binmean), binmean = NaN; binsem = NaN;end;
    deltaIOPhaseNN_allCells_binmean(b) = binmean;
    deltaIOPhaseNN_allCells_binsem(b)  = binsem;
    % for the following anova test:
    deltaIOPhaseNN_allCells_eachBin{b}      = abs(deltaIOPhaseNN_allCells(bins==b));
    deltaIOPhaseNN_allCells_eachBin_Group   = [deltaIOPhaseNN_allCells_eachBin_Group; repmat( (Edges(1+b)), length(deltaIOPhaseNN_allCells_eachBin{b}),1) ];
end
NNall.deltaIOPhaseNN_allCells_binmean = deltaIOPhaseNN_allCells_binmean;
NNall.deltaIOPhaseNN_allCells_binsem  = deltaIOPhaseNN_allCells_binsem;

%% Mean tuning difference vs. cortical distance - exp-wise
deltaIOPhaseNN_binmean_meanAcrossExps = nanmean(deltaIOPhaseNN_binmean_eachExp, 1);
deltaIOPhaseNN_binmean_stdAcrossExps =  nanstd(deltaIOPhaseNN_binmean_eachExp, 0,1);
n_vals = sum(~isnan(deltaIOPhaseNN_binmean_eachExp), 1);
deltaIOPhaseNN_binmean_semAcrossExps = deltaIOPhaseNN_binmean_stdAcrossExps./sqrt(n_vals);


%% Mean cortical distance vs. tuning difference - exp-wise
Dist_um_byDeltaIOPhase_meanAcrossExps = nanmean(Dist_um_byDeltaIOPhase_meanEachExp_mat, 1);
Dist_um_byDeltaIOPhase_stdAcrossExps  =  nanstd(Dist_um_byDeltaIOPhase_meanEachExp_mat, 0,1);
samplesizeperbin4 = sum(~isnan(Dist_um_byDeltaIOPhase_meanEachExp_mat), 1);
Dist_um_byDeltaIOPhase_semAcrossExps  = Dist_um_byDeltaIOPhase_stdAcrossExps./sqrt(samplesizeperbin4);
NNall.Dist_um_byDeltaIOPhase_meanAcrossExps = Dist_um_byDeltaIOPhase_meanAcrossExps;
NNall.Dist_um_byDeltaIOPhase_stdAcrossExps  = Dist_um_byDeltaIOPhase_stdAcrossExps;
NNall.Dist_um_byDeltaIOPhase_semAcrossExps  = Dist_um_byDeltaIOPhase_semAcrossExps;


%% Mean cortical distance vs. tuning difference for 1st and 2nd neighbor - exp-wise
D_1stNearest_perdIOP_meanAcrossExps = nanmean(D_1stNearest_perdIOP_meanEachExp, 1);
D_1stNearest_perdIOP_stdAcrossExps  =  nanstd(D_1stNearest_perdIOP_meanEachExp, 0,1);
samplesizeperbin = sum(~isnan(D_1stNearest_perdIOP_meanEachExp), 1);
D_1stNearest_perdIOP_semAcrossExps  = D_1stNearest_perdIOP_stdAcrossExps ./ sqrt(samplesizeperbin);
D_2ndNearest_perdIOP_meanAcrossExps = nanmean(D_2ndNearest_perdIOP_meanEachExp, 1);
D_2ndNearest_perdIOP_stdAcrossExps  =  nanstd(D_2ndNearest_perdIOP_meanEachExp, 0,1);
samplesizeperbin = sum(~isnan(D_2ndNearest_perdIOP_meanEachExp), 1);
D_2ndNearest_perdIOP_semAcrossExps  = D_2ndNearest_perdIOP_stdAcrossExps ./ sqrt(samplesizeperbin);
% shuffles:
% D_1stNearest_perdIOP_meanAcrossExps_shuffle = nanmean(D_1stNearest_perdIOP_meanEachExp_shuffle, 1);
% D_1stNearest_perdIOP_stdAcrossExps_shuffle  =  nanstd(D_1stNearest_perdIOP_meanEachExp_shuffle, 0,1);
% samplesizeperbin = sum(~isnan(D_1stNearest_perdIOP_meanEachExp_shuffle), 1);
% D_1stNearest_perdIOP_semAcrossExps_shuffle  = D_1stNearest_perdIOP_stdAcrossExps_shuffle ./ sqrt(samplesizeperbin);
% D_2ndNearest_perdIOP_meanAcrossExps_shuffle = nanmean(D_2ndNearest_perdIOP_meanEachExp_shuffle, 1);
% D_2ndNearest_perdIOP_stdAcrossExps_shuffle  =  nanstd(D_2ndNearest_perdIOP_meanEachExp_shuffle, 0,1);
% samplesizeperbin = sum(~isnan(D_2ndNearest_perdIOP_meanEachExp_shuffle), 1);
% D_2ndNearest_perdIOP_semAcrossExps_shuffle  = D_2ndNearest_perdIOP_stdAcrossExps_shuffle ./ sqrt(samplesizeperbin);
D_1stNearest_perdIOP_95thAcrossExps_shuffle = nanmean(D_1stNearest_perdIOP_95thEachExp_shuffle, 1);
D_1stNearest_perdIOP_05thAcrossExps_shuffle = nanmean(D_1stNearest_perdIOP_05thEachExp_shuffle, 1);
%% Mean cortical distance vs. tuning difference for 1st and 2nd neighbor - cell-wise
D_1stNearest_perdIOP_meanAllCells = nanmean(D_1stNearest_perdIOP_AllCells, 1);
D_1stNearest_perdIOP_stdAllCells  =  nanstd(D_1stNearest_perdIOP_AllCells, 0,1);
samplesizeperbin = sum(~isnan(D_1stNearest_perdIOP_AllCells), 1);
D_1stNearest_perdIOP_semAllCells  = D_1stNearest_perdIOP_stdAllCells ./ sqrt(samplesizeperbin);
D_2ndNearest_perdIOP_meanAllCells = nanmean(D_2ndNearest_perdIOP_AllCells, 1);
D_2ndNearest_perdIOP_stdAllCells  =  nanstd(D_2ndNearest_perdIOP_AllCells, 0,1);
samplesizeperbin = sum(~isnan(D_2ndNearest_perdIOP_AllCells), 1);
D_2ndNearest_perdIOP_semAllCells  = D_2ndNearest_perdIOP_stdAllCells ./ sqrt(samplesizeperbin);
% shuffles
D_1stNearest_perdIOP_meanAllCells_shuffle = nanmean(D_1stNearest_perdIOP_AllCells_shuffle, 1);
D_1stNearest_perdIOP_stdAllCells_shuffle  =  nanstd(D_1stNearest_perdIOP_AllCells_shuffle, 0,1);
samplesizeperbin = sum(~isnan(D_1stNearest_perdIOP_AllCells_shuffle), 1);
D_1stNearest_perdIOP_semAllCells_shuffle  = D_1stNearest_perdIOP_stdAllCells_shuffle ./ sqrt(samplesizeperbin);
D_2ndNearest_perdIOP_meanAllCells_shuffle = nanmean(D_2ndNearest_perdIOP_AllCells_shuffle, 1);
D_2ndNearest_perdIOP_stdAllCells_shuffle  =  nanstd(D_2ndNearest_perdIOP_AllCells_shuffle, 0,1);
samplesizeperbin = sum(~isnan(D_2ndNearest_perdIOP_AllCells_shuffle), 1);
D_2ndNearest_perdIOP_semAllCells_shuffle  = D_2ndNearest_perdIOP_stdAllCells_shuffle ./ sqrt(samplesizeperbin);
D_1stNearest_perdIOP_95thAcrossCells_shuffle = prctile(D_1stNearest_perdIOP_AllCells_shuffle, 97.5, 1);
D_1stNearest_perdIOP_05thAcrossCells_shuffle = prctile(D_1stNearest_perdIOP_AllCells_shuffle,  2.5, 1);

%%
% PrefIOPhase_diff_Nearest1st_allCells = vertcat(NNexp.PrefIOPhase_diff_Nearest1st_eachExp);
% 
%     NNexp(e).PrefIOPhase_diff_Nearest1st_count_eachExp;
%     
%     figure;
%     histogram(PrefIOPhase_diff_Nearest1st_allCells, bin_hist, 'Normalization','probability');
%     
%     figure;
%     histogram(abs(deltaIOPhaseNN_allCells), bin_hist, 'Normalization','probability');


%% Statistic tests

disp('ANOVA test across all cells: deltaIOPhaseNN_allCells_eachBin')
deltaIOPhaseNN_allCells_eachBin_mat = cell2mat(deltaIOPhaseNN_allCells_eachBin);
try
    [p, stats, tbl] = perform_anova1(deltaIOPhaseNN_allCells_eachBin_mat, deltaIOPhaseNN_allCells_eachBin_Group );
catch
    disp(['  /!\ something went wrong: probably not enough data for anova']);
    p = nan; stats = nan; tbl = nan;
end
deltaIOPhaseNN_allCells_eachBin_anova1_p     = p;


if ismember(2, Shuffling)
    deltaIOPhaseNN_allCells_bootstr_binmean = nan(RSiters, nBins);

    for rsi = 1 : RSiters
        Dist_um_allCells_bootstr = datasample(Dist_um_allCells, length(Dist_um_allCells));

        [~,~,bins] = histcounts(Dist_um_allCells_bootstr, Edges);
        for b = 1 : nBins
            binmean = circ_mean((abs(deltaIOPhaseNN_allCells(bins==b))*pi/180))*180/pi ;
            [~, s0] =  circ_std( abs(deltaIOPhaseNN_allCells(bins==b))*pi/180 ) ;
            binsem = s0*180/pi / sqrt(sum(bins==b));
            if isempty(binmean), binmean = NaN; binsem = NaN;end;
            deltaIOPhaseNN_allCells_bootstr_binmean(rsi,b) = binmean;
    %         deltaIOPhaseNN_allCells_binsem(b)  = binsem;
        end
    end

    deltaIOPhaseNN_allCells_bootstr_binmean_prctile05 = prctile(deltaIOPhaseNN_allCells_bootstr_binmean,  2.5, 1);
    deltaIOPhaseNN_allCells_bootstr_binmean_prctile95 = prctile(deltaIOPhaseNN_allCells_bootstr_binmean, 97.5, 1);
else
    deltaIOPhaseNN_allCells_bootstr_binmean_prctile05 = nan(1, nBins);
    deltaIOPhaseNN_allCells_bootstr_binmean_prctile95 = nan(1, nBins);
end

NNall.deltaIOPhaseNN_allCells_bootstr_binmean_prctile05 = deltaIOPhaseNN_allCells_bootstr_binmean_prctile05;
NNall.deltaIOPhaseNN_allCells_bootstr_binmean_prctile95 = deltaIOPhaseNN_allCells_bootstr_binmean_prctile95;
    
% % Linear regression with all ROIs pooled together
%  % to check whether regression slope is significantly different from 0 (check if 95% CI contains 0)
% [R,P,RL,RU] = corrcoef(Dist_um_allCells, abs(deltaIOPhaseNN_allCells), 'alpha',0.05);
% deltaIOPhaseNN_allCells_eachBin_Corrcoef_R         = R(2);
% deltaIOPhaseNN_allCells_eachBin_Corrcoef_p         = P(2);
% deltaIOPhaseNN_allCells_eachBin_Corrcoef_LowerCI95 = RL(2);
% deltaIOPhaseNN_allCells_eachBin_Corrcoef_UpperCI95 = RU(2);


disp('ANOVA test exp by exp: deltaIOPhaseNN_binmean_eachExp')
Groups = cellstr(num2str(Edges(:))); Groups = Groups(2:end-1);
try
    [p, stats, tbl] = perform_anova1(deltaIOPhaseNN_binmean_eachExp, Groups );
catch
    disp(['  /!\ something went wrong: probably not enough data for anova']);
    p = nan; stats = nan; tbl = nan;
end
deltaIOPhaseNN_binmean_eachExp_anova1_p     = p;
% NNall.deltaIOPhaseNN_binmean_acrossCells_eachExp_anova1_p     = p;
% NNall.deltaIOPhaseNN_binmean_acrossCells_eachExp_anova1_stats = stats;
% NNall.deltaIOPhaseNN_binmean_acrossCells_eachExp_anova1_tbl   = tbl;


disp('ANOVA test exp by exp: Dist_um_byDeltaIOPhase')
Groups = deltaIOPhases;
try
    [p, stats, tbl] = perform_anova1(Dist_um_byDeltaIOPhase_meanEachExp_mat, Groups );
catch
    disp(['  /!\ something went wrong: probably not enough data for anova']);
    p = nan; stats = nan; tbl = nan;
end
Dist_um_byDeltaIOPhase_meanEachExp_anova1_p     = p;
% NNall.Dist_um_byDeltaIOPhase_meanEachExp_anova1_p     = p;
% NNall.Dist_um_byDeltaIOPhase_meanEachExp_anova1_stats = stats;
% NNall.Dist_um_byDeltaIOPhase_meanEachExp_anova1_tbl   = tbl;


disp('ANOVA test exp by exp: Dist_um_1stNeighbor_byDeltaIOPhase')
Groups = deltaIOPhases;
try
    [p, stats, tbl] = perform_anova1(D_1stNearest_perdIOP_meanEachExp, Groups );
catch
    disp(['  /!\ something went wrong: probably not enough data for anova']);
    p = nan; stats = nan; tbl = nan;
end
D_1stNearest_perdIOP_meanEachExp_anova1_p     = p;


disp('ANOVA test exp by exp: Dist_um_2ndNeighbor_byDeltaIOPhase')
try
    [p, stats, tbl] = perform_anova1(D_2ndNearest_perdIOP_meanEachExp, Groups );
catch
    disp(['  /!\ something went wrong: probably not enough data for anova']);
    p = nan; stats = nan; tbl = nan;
end
D_2ndNearest_perdIOP_meanEachExp_anova1_p     = p;


disp('ANOVA test across all cells: Dist_um_1stNeighbor_byDeltaIOPhase')
Groups = deltaIOPhases;
try
    [p, stats, tbl] = perform_anova1(D_1stNearest_perdIOP_AllCells, Groups );
catch
    disp(['  /!\ something went wrong: probably not enough data for anova']);
    p = nan; stats = nan; tbl = nan;
end
D_1stNearest_perdIOP_AllCells_anova1_p     = p;


disp('ANOVA test across all cells: Dist_um_2ndNeighbor_byDeltaIOPhase')
try
    [p, stats, tbl] = perform_anova1(D_2ndNearest_perdIOP_AllCells, Groups );
catch
    disp(['  /!\ something went wrong: probably not enough data for anova']);
    p = nan; stats = nan; tbl = nan;
end
D_2ndNearest_perdIOP_AllCells_anova1_p     = p;

%% Figures

Area = aDGDap(1).Area;
Edges_plot = Edges(2:end-1) - diff(Edges(1:end-1))/2;
NNall.BinEdges_plot = Edges_plot; NNall.BinEdges = Edges;
CC = load_DGD_colors(2);
color_sign = [1 1 1]*0.5;

if any(plot_figs == 1)
    hf(1) = figure;
    figname{1} = ['deltaIOPhase_vs_distance' '_Area-' Area '_SF-' SF '_Ori-' num2str(Angles(1)) ];
    hold on;
    title({['Area ' Area ' - SF ' SF ' - Ori ' num2str(Angles(1)) ];...
           ['anova(exp-by-exp) p=' num2str(deltaIOPhaseNN_binmean_eachExp_anova1_p)];...
           ['anova(all-cells) p='  num2str(deltaIOPhaseNN_allCells_eachBin_anova1_p)]})
    % exp-by-exp:
%     hp1 = plot( Edges_plot, deltaIOPhaseNN_binmean_meanAcrossExps, '-o', 'Color','k', 'LineWidth',2);
%     supererr(   Edges_plot, deltaIOPhaseNN_binmean_meanAcrossExps, [], deltaIOPhaseNN_binmean_semAcrossExps')
    % all-cells:
    hp2 = plot( Edges_plot, deltaIOPhaseNN_allCells_binmean, '-o', 'Color',CC.(Area){1}, 'LineWidth',2);
    supererr(   Edges_plot, deltaIOPhaseNN_allCells_binmean, [], deltaIOPhaseNN_allCells_binsem)
    % CI 95% of shuffles:
%     hc1 = plot(    Edges_plot, mean(deltaIOPhaseNN_bootstr_binmean_prctile05_eachExp,1), '--', 'Color','k');
%           plot(    Edges_plot, mean(deltaIOPhaseNN_bootstr_binmean_prctile95_eachExp,1), '--', 'Color','k');
    hc2 = plot(    Edges_plot, deltaIOPhaseNN_allCells_bootstr_binmean_prctile05, '--', 'Color',color_sign);
          plot(    Edges_plot, deltaIOPhaseNN_allCells_bootstr_binmean_prctile95, '--', 'Color',color_sign);

    set(gca, 'XTick',Edges(1:end-1), 'XTickLabel',Edges(1:end-1));
%     legend([hp1,hp2,hc2],'exp-by-exp','all-cells','shuffle 95th CI','Location','SouthEast')
    legend([hp2,hc2],'all-cells','shuffle 95th CI','Location','SouthEast')
end



%% Save figures
if save_figures == 1  &&  plot_figs > 0
    fig_dir = [pwd filesep 'figs' filesep];
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
    for f = 1 : length(figname)
        if ~isempty(figname{f}) 
            saveas(hf(plot_figs(f)), [fig_dir figname{plot_figs(f)} '.png'], 'png');
            saveas(hf(plot_figs(f)), [fig_dir figname{plot_figs(f)} '.fig'], 'fig');
            disp([ ' * ' figname{f} ' saved']);
        end
    end
end
