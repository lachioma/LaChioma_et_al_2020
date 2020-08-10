function [p, p_pairwise, stats, tbl_disp] = perform_anova1(Y, Group, TestType, CType, display_anova, display_multcompare)

% Input Y: each column is a group.

if nargin < 2
    Group = [];
end
if nargin < 3 || isempty(TestType)
    TestType = 'ANOVA';
    TestType = 'KW'; % Kruskal-Wallis
end
if nargin < 4 || isempty(CType)
    CType = 'tukey-kramer';
end
if nargin < 5
    display_anova = 'off';
end
if nargin < 6
    display_multcompare = 'off';
end
TestType = upper(TestType);

if iscell(Y)
    if length(Y) ~= length(Group)
        disp(' Input Y (cell array) and Group must have the same length, aborting');
        return
    end
    for i = 1 : length(Y)
        Y{i} = reshape(Y{i},[],1);
        if iscell(Group)
            str_group{i} = Group{i};
            Group_cell{i,1} = repmat(Group(i), length(Y{i}), 1); %#ok<*AGROW>
        else
            str_group{i} = num2str(Group(i));
            Group_cell{i,1} = repmat({Group(i)}, length(Y{i}), 1);
        end
    end
    Y = reshape(Y,[],1);
    Y_vect = cell2mat(Y);
    Group_vect = vertcat(Group_cell{:});
else
    Y_vect     = reshape(Y,[],1);
    Group_vect = reshape(Group,[],1);
end

if iscell(Y)
    for i = 1 : length(Y)
        Y_mean(i)  = nanmean(Y{i});
        Y_sem(i)   =  nanstd(Y{i},0)/sqrt(sum(~isnan(Y{i})));
        Y_med(i)   = nanmedian(Y{i});
        Y_25th(i)  = prctile(Y{i}, 25);
        Y_75th(i)  = prctile(Y{i}, 75);
        % make strings:
        str_mean{i}  = num2str(Y_mean(i));
        str_sem{i}   = num2str(Y_sem(i));
        str_med{i}   = num2str(Y_med(i));
        str_25th{i}  = num2str(Y_25th(i));
        str_75th{i}  = num2str(Y_75th(i));
    end
    tbl_means =  [{'Group'},  {'Mean'}, {'SEM'}, {'Median'},{'25th perc'},{'75th perc'};...
                  str_group', str_mean', str_sem', str_med',  str_25th',   str_75th', ];
    disp(tbl_means);
end
switch TestType
    case 'ANOVA'
        [p,tbl_anova,stats]  = anova1(Y_vect, Group_vect, display_anova);
        disp([' Test ANOVA one-way']);
        disp([' ' tbl_anova{1,5} '(' num2str(tbl_anova{2,3}) ',' num2str(tbl_anova{3,3}) ') = ' num2str(tbl_anova{2,5}) '  -  p = ' num2str(p)])
    case 'KW'
        [p,tbl_anova,stats]  = kruskalwallis(Y_vect, Group_vect, display_anova);
        disp([' Test Kruskal-Wallis']);
        disp([' ' tbl_anova{1,5} '(' num2str(tbl_anova{2,3}) ') = ' num2str(tbl_anova{2,5}) '  -  p = ' num2str(p)])
end

if p==0
    disp([' Test ' TestType ' not possible to perform !'])
    p = nan;
    p_pairwise = nan;
    stats = [];
    tbl_disp = '';
    return
end
    



if p < 0.05
    
    [c,~,~,gnames] = multcompare(stats, 'CType',CType, 'Display',display_multcompare);
    if ischar(gnames)
        gnames = cellstr(gnames);
    end
%     tbl = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,6)), num2cell(c(:,3:5))];
%     tbl_disp = [{'Group1','Group2','p-value','CI%95 lower','delta means','CI%95 upper'}; tbl];
    tbl = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,6))];
    tbl_disp = [{'Group1','Group2','p-value'}; tbl];

    % p_pairwise = c(:,6);
    nGroups = length(unique(Group));
    p_pairwise = ones(nGroups);
    p_pairwise(sub2ind([nGroups,nGroups],c(:,1),c(:,2))) = c(:,6);
    p_pairwise = min(p_pairwise, p_pairwise');
    % P( P>0.05 ) = nan;

    
    disp([' Post-hoc multiple pairwise comparison tests (' CType '):'])
    disp(tbl_disp)
    
else
    nGroups = length(unique(Group));
    p_pairwise = ones(nGroups);
    tbl_disp = '';
end