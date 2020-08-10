function selROIs = above_RsquaredThreshold(TC, thr_Rsquared, GoF_fieldname)

if nargin < 3
    GoF_fieldname = 'Rsquared_meanTC';
end

selROIs = [];

cnt = 0;
for nr = 1 : length(TC)
    if ~isempty(fieldnames(TC(nr).GoF)) && ...
            TC(nr).GoF.(GoF_fieldname) > thr_Rsquared
        cnt = cnt + 1;
        selROIs(cnt,1) = nr; %#ok<AGROW>
    end
end