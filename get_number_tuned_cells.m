function N_resp = get_number_tuned_cells(Fit, Areas, RDStype1, RDStype2, RDStype12, ...
                    thr_Rsquared, Rois_to_intersect)

% Rois_to_intersect

nAreas = numel(Areas);

if nargin < 7 || isempty(Rois_to_intersect)
     Rois_to_intersect = {};
else
    % verify that Rois_to_intersect is in the form {nAreas, nExps}:
    [m,n] = size(Rois_to_intersect);
    assert(m == nAreas, ' ''Rois_to_intersect'' must be cell array {nAreas, nExps}');
    for ax = 1 : nAreas
        nExps(ax) = length(Fit(ax).(RDStype12)); %#ok<*AGROW>
    end
    assert(n == max(nExps), ' ''Rois_to_intersect'' must be cell array {nAreas, nExps}');
%     for ax = 1 : nAreas
%         Area = Areas{ax};
%         ai = find(strcmp(Area,{Fit.Area}));
%         nExps = length(Fit(ai).(RDStype12));
%         isnotempty = @(c) ~isempty(c);
%         nExps_rois_to_intersect = sum(cellfun(isnotempty, Rois_to_intersect(ax,:)));
%         assert(nExps_rois_to_intersect == nExps, ' ''Rois_to_intersect'' must be cell array {nAreas, nExps}');
%     end
end

for ax = 1 : nAreas
    Area = Areas{ax};
    ai = find(strcmp(Area,{Fit.Area}));
    
    nExps = length(Fit(ai).(RDStype12));

    for e = 1 : nExps
        
%         respRDS_1  = above_stdThreshold_DGD(aD_RDS(ai).(RDStype1)(e) , thr_std, thr_type, [], thr_rep_fraction_rds);
%         respRDS_2  = above_stdThreshold_DGD(aD_RDS(ai).(RDStype2)(e) , thr_std, thr_type, [], thr_rep_fraction_rds);
        respRDS_1  = above_RsquaredThreshold(Fit(ai).(RDStype1)(e).TC, -Inf, 'Rsquared_meanTC');
        respRDS_2  = above_RsquaredThreshold(Fit(ai).(RDStype2)(e).TC, -Inf, 'Rsquared_meanTC');
        respRDS_1_only = respRDS_1( ~ismember(respRDS_1, respRDS_2) );
        respRDS_2_only = respRDS_2( ~ismember(respRDS_2, respRDS_1) );
        respRDS_any =     union(respRDS_1, respRDS_2);
        respRDS     = intersect(respRDS_1, respRDS_2);
        rsqRDS_1      = above_RsquaredThreshold(Fit(ai).(RDStype1)(e).TC, thr_Rsquared, 'Rsquared_meanTC');
        rsqRDS_2      = above_RsquaredThreshold(Fit(ai).(RDStype2)(e).TC, thr_Rsquared, 'Rsquared_meanTC');
        rsqRDS_12_1   = above_RsquaredThreshold(Fit(ai).(RDStype12  )(e).TC, thr_Rsquared, 'Rsquared_meanTC_1');
        rsqRDS_12_2   = above_RsquaredThreshold(Fit(ai).(RDStype12  )(e).TC, thr_Rsquared, 'Rsquared_meanTC_2');
        
        respRDS_12_1_tmp  = above_RsquaredThreshold(Fit(ai).(RDStype12 )(e).TC, -Inf, 'Rsquared_meanTC_1');
        respRDS_12_2_tmp  = above_RsquaredThreshold(Fit(ai).(RDStype12 )(e).TC, -Inf, 'Rsquared_meanTC_2');
        respRDS_12        =  intersect(respRDS_12_1_tmp, respRDS_12_2_tmp);

        if ~isempty(Rois_to_intersect) ...&& ...
           ...~isempty(Rois_to_intersect{ai,e}) % this line is necessary because when there are no RPS cells for a given exp, Rois_to_intersect{ai,e} is defined to be [].
           % errata corrige: that line needed to be removed because when
           % there are no RPS cells for a given exp, the following
           % intersections were not done!
            Rois_to_intersect_thisExp = Rois_to_intersect{ai,e};
            respRDS_1       = intersect(respRDS_1       , Rois_to_intersect_thisExp);
            respRDS_2       = intersect(respRDS_2       , Rois_to_intersect_thisExp);
            respRDS_1_only  = intersect(respRDS_1_only  , Rois_to_intersect_thisExp);
            respRDS_2_only  = intersect(respRDS_2_only  , Rois_to_intersect_thisExp);
            respRDS_any     = intersect(respRDS_any     , Rois_to_intersect_thisExp);
            respRDS         = intersect(respRDS         , Rois_to_intersect_thisExp);
            rsqRDS_1        = intersect(rsqRDS_1        , Rois_to_intersect_thisExp);
            rsqRDS_2        = intersect(rsqRDS_2        , Rois_to_intersect_thisExp);
            rsqRDS_12_1     = intersect(rsqRDS_12_1     , Rois_to_intersect_thisExp);
            rsqRDS_12_2     = intersect(rsqRDS_12_2     , Rois_to_intersect_thisExp);
            respRDS_12      = intersect(respRDS_12      , Rois_to_intersect_thisExp);
        end
        
        respRDS_1_only_untuned = respRDS_1_only( ~ismember(respRDS_1_only, rsqRDS_1) );
        respRDS_1_only_tuned   = intersect( respRDS_1_only, rsqRDS_1 );
        respRDS_2_only_untuned = respRDS_2_only( ~ismember(respRDS_2_only, rsqRDS_2) );
        respRDS_2_only_tuned   = intersect( respRDS_2_only, rsqRDS_2 );
        
        n_respRDS_1{ax}(e,1) = length(respRDS_1);
        n_respRDS_2{ax}(e,1) = length(respRDS_2);
        n_respRDS_1_only_eachExp{ax}(e,1)   = length(respRDS_1_only);
        n_respRDS_2_only_eachExp{ax}(e,1)   = length(respRDS_2_only);
        n_respRDS_1_only_untuned_eachExp{ax}(e,1) = length(respRDS_1_only_untuned);
        n_respRDS_1_only_tuned_eachExp{ax}(e,1)   = length(respRDS_1_only_tuned);
        n_respRDS_2_only_untuned_eachExp{ax}(e,1) = length(respRDS_2_only_untuned);
        n_respRDS_2_only_tuned_eachExp{ax}(e,1)   = length(respRDS_2_only_tuned);
                
        respRDS_12_tuned_12 = mintersect(respRDS_12, rsqRDS_12_1, rsqRDS_12_2);
        respRDS_12_tuned_to_notboth = respRDS_12( ~ismember(respRDS_12, respRDS_12_tuned_12) );
        respRDS_12_tuned_1_only = intersect(respRDS_12_tuned_to_notboth, rsqRDS_12_1);
        respRDS_12_tuned_2_only = intersect(respRDS_12_tuned_to_notboth, rsqRDS_12_2);
        respRDS_12_tuned    = unionm(respRDS_12_tuned_12, respRDS_12_tuned_1_only, respRDS_12_tuned_2_only);
        respRDS_12_untuned  = respRDS_12( ~ismember(respRDS_12, unionm(respRDS_12_tuned_12, respRDS_12_tuned_1_only, respRDS_12_tuned_2_only)) );
        
        respRDS_any_untuned = unionm(respRDS_1_only_untuned, respRDS_2_only_untuned, respRDS_12_untuned);
        respRDS_any_tuned   = unionm(respRDS_1_only_tuned, respRDS_2_only_tuned, respRDS_12_tuned_12, respRDS_12_tuned_1_only, respRDS_12_tuned_2_only);
        
        n_respRDS_any_eachExp{ax}(e,1)            = length(respRDS_any);
        n_respRDS_any_untuned_eachExp{ax}(e,1)    = length(respRDS_any_untuned);
        n_respRDS_any_tuned_eachExp{ax}(e,1)      = length(respRDS_any_tuned);
        n_respRDS_1_only_eachExp{ax}(e,1)         = length(respRDS_1_only);
        n_respRDS_2_only_eachExp{ax}(e,1)         = length(respRDS_2_only);

        n_respRDS_12_eachExp{ax}(e,1)              = length(respRDS_12);
        n_respRDS_12_tuned_any_eachExp{ax}(e,1)    = length(respRDS_12_tuned);
        n_respRDS_12_untuned_eachExp{ax}(e,1)      = length(respRDS_12_untuned);
        n_respRDS_12_tuned_12_eachExp{ax}(e,1)     = length(respRDS_12_tuned_12);
        n_respRDS_12_tuned_1_only_eachExp{ax}(e,1) = length(respRDS_12_tuned_1_only);
        n_respRDS_12_tuned_2_only_eachExp{ax}(e,1) = length(respRDS_12_tuned_2_only);
        
        n_respRDS_any_tuned_1_only_eachExp{ax}(e,1)  = n_respRDS_1_only_tuned_eachExp{ax}(e,1) + n_respRDS_12_tuned_1_only_eachExp{ax}(e,1);
        n_respRDS_any_tuned_2_only_eachExp{ax}(e,1)  = n_respRDS_2_only_tuned_eachExp{ax}(e,1) + n_respRDS_12_tuned_2_only_eachExp{ax}(e,1);
        
        n_totCells_eachExp{ax}(e,1)            = length(Fit(ai).(RDStype1)(e).TC);
    end    
end

% Get all variable containing 'n_resp':
namesWorkspace = who; outStr = regexpi(namesWorkspace, 'n_resp'); ind = ~cellfun('isempty',outStr); varlist = namesWorkspace(ind);
fieldNames = [{'fieldNames'}; varlist]; % v2struct wants 'fieldNames' in this way
N_resp = v2struct(fieldNames);
% Get n_totCells_eachExp:
N_resp.n_totCells_eachExp = n_totCells_eachExp;
