




%% DI distribution (Fig. 2B)
% You just need to load variable Disp:
% load([pwd '\Vars\Disp_DGD_V1_LM_RL.mat']);

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


%% ODI (Fig. 3A)
% You just need to load variable OD
% load([pwd '\Vars\OD_DG_V1_LM_RL.mat']);


Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
CC = load_DGD_colors(1);
edges = -1:2/7:1; dx = (edges(2)-edges(1)); xticks = edges(2:end)-dx/2;
xrange = [-1.1 1.1]; yrange = [-0.8 35]; dy = [10]; 
c = 2; % color shade

for ax = 1 : nAreas
    Area = Areas{ax};
    ODI = OD.(Area);
    CountODI_EachExp_perc = 100*cat(1, ODI.CountODI_EachExp_perc);
    CountODI_meanAcrossExps_perc = mean(CountODI_EachExp_perc, 1);
    CountODI_stdAcrossExps_perc  = std (CountODI_EachExp_perc, 0, 1);
    nExps_tmp = size(CountODI_EachExp_perc, 1);
    CountODI_semAcrossExps_perc  = CountODI_stdAcrossExps_perc / sqrt(nExps_tmp);


    % Figure: ODI distribution across exps
    figure;
        fig_pos = get(gcf, 'position'); fig_pos(3) = fig_pos(3)*1.3; %fig_pos(2) = fig_pos(2)-fig_pos(4)*0.5;
        set(gcf, 'position',fig_pos); % [left, bottom, width, height]

    hb = superbar( xticks, CountODI_meanAcrossExps_perc,...
        'E', CountODI_semAcrossExps_perc,...
        'BarWidth',0.9*dx, 'BarFaceColor', CC.(Areas{ax}){c}, 'BarEdgeColor','none',...
        'ErrorbarRelativeWidth',0.0, 'ErrorbarLineWidth',3);

    set(gca, 'XTick',[-1 0 1],'XLim',xrange,'XTickLabel',{'-1','0','1'}, 'TickLength',[0.02 0.025])
    set(gca, 'YTick',[0:dy:200], 'YLim',yrange);
    xlabel('Ocular Dominance Index (ODI)');
    ylabel('Proportion of cells (%)')
    hh = get( hb(1), 'BaseLine');
    hh.LineStyle = 'none';

end


%% DI vs ODI (Fig. 3B)
% You just need to load variables OD and DIvar
% load([pwd '\Vars\OD_DG_V1_LM_RL.mat']);
% load([pwd '\Vars\DIvar_DGD_V1_LM_RL.mat']);


Areas = {'V1','LM','RL'};
SFs   = {'0.01'};
nSFs = length(SFs);
nAreas = numel(Areas);
CC  = load_DGD_colors(1); CC2  = load_DGD_colors(1);
acircle = 42;
rangex = [-1.05 1.05]; dx = 0.5; ticklabels = {'-1','','0','','1'};
rangey = [-0.05 1.05];
errorbar_width = 0.1; markesize = 10;

Edges = linspace(-1, 1, 7+1);
nBins = length(Edges) - 1;
xValBins = mean([Edges(1:end-1);Edges(2:end)]);
p_anova = nan(nAreas,1);

for ax = 1 : nAreas
    Area = Areas{ax};
    ODI = cat(1, OD.(Area).ODI);
    MI  = abs(ODI);
    resp_ODI = find(~isnan(ODI));
    
    ODIselall = []; DIselall = [];
    
    for i = 1 : nSFs
        SF = SFs{i};
        s = find( strcmp(DIvar.SF, SF) );
        
        DI = DIvar.(Area)(s).DI;
        
        resp_DGDsf = DIvar.(Area)(s).resp_DGDthisSF;
        resp = intersect(resp_DGDsf, resp_ODI);
        
        ODIsel = ODI(resp);
        MIsel  = MI(resp);
        BIsel  = 1-MIsel;
        DIsel  = DI(resp);
    
        ODIselall = [ODIselall; ODIsel];
        DIselall  = [DIselall ; DIsel ];
        
    end
        
    ODIsel_binmean = nan(nBins,1);
    ODIsel_binsem  = nan(nBins,1);

    [~,~,bins] = histcounts(ODIselall, Edges);
    for b = 1:nBins
        binmean = mean( DIselall(bins==b) ) ;
        binsem  =  std( DIselall(bins==b) ) / sqrt(sum(bins==b));
        if isempty(binmean), binmean = NaN; binsem = NaN;end;
        ODIsel_binmean(b) = binmean;
        ODIsel_binsem(b)  = binsem;
    end

    %%%% ANOVA test:
    Group_all = [];
    DIsel_allbins = [];
    [~,~,bins] = histcounts(ODIselall, Edges);
    for b = 1:nBins
        DIsel_thisbin =  DIselall(bins==b) ;
        DIsel_allbins = [DIsel_allbins; DIsel_thisbin];
        Group_thisbin = repmat(b, length(DIsel_thisbin),1);
        Group_all = [Group_all; Group_thisbin];
    end
    TestType = 'KW';
    p_anova(ax) = perform_anova1(DIsel_allbins, Group_all, TestType);
    %%%%
    
    
    [rhoMIDI,pMIDI] = corr( DIselall, (ODIselall), 'type', 'Pearson' );
    


    figure; hold on;
        fig_pos = get(gcf, 'position'); fig_pos(3) = fig_pos(3)*1.0; %fig_pos(2) = fig_pos(2)-fig_pos(4)*0.5;
        set(gcf, 'position',fig_pos); % [left, bottom, width, height]

    scatter(ODIselall, DIselall, acircle, CC.(Area){2}, 'filled')
        
    supererr(xValBins, ODIsel_binmean, [], ODIsel_binsem,...
        errorbar_width, 'LineWidth',2, 'Color',CC2.(Area){1});
    plot( xValBins, ODIsel_binmean,...
        'LineWidth',2,'Color',CC2.(Area){1},'Marker','o','MarkerSize',markesize,'MarkerFaceColor',CC.(Area){1})


    xlabel('ODI'); 
    ylabel('DI');
    xlim(rangex); ylim(rangey);
    set(gca, 'XTick',[-1:dx:1],'YTick',[0:dx:1],'XTickLabel',ticklabels);

    str_title = {...
           ['Area ' Area ' - SF:' ' all pooled' ];...
           ['CorrCoeff = ' num2str(rhoMIDI,'%3.2f') ' (p=' num2str(pMIDI,'%4.3f') ')' ]};
    title(str_title, 'Interpreter','none');
end


%% Facilitation index (Fig. 3C)
% You just need to load variable R:
% load([pwd '\Vars\R__DGD_V1_LM_RL.mat']);

Which_R = 'Ratio_F_DGD_F_DG_sumLR'; yrange = [0 3]; dy = 0.5;
% Which_R = 'Ratio_F_DGD_F_DG_maxLR'; yrange = [0 4];  dy = 0.5;

% Which_Mean = 'exps_pooled'; %not implemented
Which_Mean = 'exp_wise';

% Bar_Orientation = 'v' ;
Bar_Orientation = 'h' ;

Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
SFs   = {'0.01'};
nSFs = numel(SFs);
CC = load_DGD_colors(2);

for ax = 1 : nAreas
    Area = Areas{ax};
    
    R_meanAcrossExps = nan(nSFs,1);
    R_SEMAcrossExps  = nan(nSFs,1);
    R_AllSFs=[]; Group=[];
            
    figure; hold on;
    
    for i = 1 : nSFs
        SF = SFs{i};
        s  = find(strcmp({R.(Area).SF}, SF));

        if     strcmp(Which_R, 'Ratio_F_DGD_F_DG_sumLR')
            if     strcmp(Which_Mean, 'exp_wise')
                R_meanAcrossExps(i) = R.(Area)(s).Ratio_F_DGD_F_DG_sumLR_meanAcrossExps;
                R_SEMAcrossExps(i)  = R.(Area)(s).Ratio_F_DGD_F_DG_sumLR_semAcrossExps;
                R_singlevals        = R.(Area)(s).Ratio_F_DGD_F_DG_sumLR_meanEachExp;
            end
        elseif strcmp(Which_R, 'Ratio_F_DGD_F_DG_maxLR')
            if     strcmp(Which_Mean, 'exp_wise')
                R_meanAcrossExps(i) = R.(Area)(s).Ratio_F_DGD_F_DG_maxLR_meanAcrossExps;
                R_SEMAcrossExps(i)  = R.(Area)(s).Ratio_F_DGD_F_DG_maxLR_semAcrossExps;
                R_singlevals        = R.(Area)(s).Ratio_F_DGD_F_DG_maxLR_meanEachExp;
            end
        end
        
        
        % for following Significance test
        R_AllSFs = [R_AllSFs; R_singlevals];
        Group_tmp = repmat({SF}, numel(R_singlevals),1);
        Group = [Group; Group_tmp];
        
    end
    
    [~,P] = perform_anova1(R_AllSFs, Group);
    P( P>0.05 ) = nan;
%     P = nan;
    
%         hb = bar( i, R_meanAcrossExps, R_SEMAcrossExps, 'FaceColor',CC.(Area){i});
        superbar( [1:nSFs], R_meanAcrossExps,...
            'E', R_SEMAcrossExps,...
            'Orientation', Bar_Orientation,...
            'BarWidth',[], 'BarFaceColor', CC.(Areas{ax}), 'BarEdgeColor','none',...
            'ErrorbarRelativeWidth',0.0, 'ErrorbarLineWidth',4,...
            'P', P,...
            'PStarShowNS',false,'PLineWidth',1,'PStarBackgroundColor','none','PLineBackingColor','none',...
            'PStarFontSize',10, 'PLineOffset',0.2 ...
        );

    
    if strcmp(Bar_Orientation,'v')
        set(gca, 'YLim',yrange,'YTick',[0:dy:5]);
        set(gca, 'XTick',[],'XTickLabel',[])
    else
        set(gca, 'XLim',yrange,'XTick',[0:dy:5]);
        set(gca, 'YTick',[],'YTickLabel',[])
    end
    title('Facilitation index')

end

%% Suppression index (Fig. 3C)
% You just need to load variable R:
% load([pwd '\Vars\R__DGD_V1_LM_RL.mat']);

% Which_R = 'Ratio_F_DGD_min_F_DG_sumLR'; yrange = [0 0.6];
Which_R = 'Ratio_F_DGD_min_F_DG_maxLR'; yrange = [0 1]; dy = 0.25;

% Bar_Orientation = 'v' ;
Bar_Orientation = 'h' ;

Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
SFs   = {'0.01'};
nSFs = numel(SFs);
CC = load_DGD_colors(2);

for ax = 1 : nAreas
    Area = Areas{ax};
    
    R_meanAcrossExps = nan(nSFs,1);
    R_SEMAcrossExps  = nan(nSFs,1);
    R_AllSFs=[]; Group=[];
    
    figure; hold on;
    
    for i = 1 : nSFs
        SF = SFs{i};
        s  = find(strcmp({R.(Area).SF}, SF));

        if     strcmp(Which_R, 'Ratio_F_DGD_min_F_DG_sumLR')
            R_meanAcrossExps(i) = R.(Area)(s).Ratio_F_DGD_min_F_DG_sumLR_meanAcrossExps;
            R_SEMAcrossExps(i)  = R.(Area)(s).Ratio_F_DGD_min_F_DG_sumLR_semAcrossExps;
            R_singlevals        = R.(Area)(s).Ratio_F_DGD_min_F_DG_sumLR_meanEachExp';
        elseif strcmp(Which_R, 'Ratio_F_DGD_min_F_DG_maxLR')
            R_meanAcrossExps(i) = R.(Area)(s).Ratio_F_DGD_min_F_DG_maxLR_meanAcrossExps;
            R_SEMAcrossExps(i)  = R.(Area)(s).Ratio_F_DGD_min_F_DG_maxLR_semAcrossExps;
            R_singlevals        = R.(Area)(s).Ratio_F_DGD_min_F_DG_maxLR_meanEachExp';
        end

        
        % for following Significance test
        R_AllSFs = [R_AllSFs; R_singlevals];
        Group_tmp = repmat({SF}, numel(R_singlevals),1);
        Group = [Group; Group_tmp];
        
    end

    [~,P] = perform_anova1(R_AllSFs, Group);
    P( P>0.05 ) = nan;
%     P = nan;
    
%         hb = bar( i, R_meanAcrossExps, R_SEMAcrossExps, 'FaceColor',CC.(Area){i});
        hb = superbar( [1:nSFs], R_meanAcrossExps,...
            'E', R_SEMAcrossExps,...
            'Orientation', Bar_Orientation,...
            'BarWidth',[], 'BarFaceColor', 'none', 'BarEdgeColor',CC.(Areas{ax}),'BarLineWidth',10,...
            'ErrorbarRelativeWidth',0.0, 'ErrorbarLineWidth',4,...
            'P', P,...
            'PStarShowNS',false,'PLineWidth',1,'PStarBackgroundColor','none','PLineBackingColor','none',...
            'PStarFontSize',10, 'PLineOffset',0.1 ...
        );
    
    
    
    hh = get( hb(1), 'BaseLine'); hh.LineStyle = 'none';
    if strcmp(Bar_Orientation,'v')
        set(gca, 'YLim',yrange,'YTick',[-2:0.25:5]);
        set(gca, 'XTick',[],'XTickLabel',[])
    else
        set(gca, 'XLim',yrange,'XTick',[0:dy:5]);
        set(gca, 'YTick',[],'YTickLabel',[])
    end 
    title('Suppression index')

end



%% Facilitation and suppression: F_DGD vs. F_DG (Fig. 3C)
% You just need to load variable R:
% load([pwd '\Vars\R__DGD_V1_LM_RL.mat']);

Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
SFs   = {'0.01'};
nSFs = numel(SFs);
thr_std = 4;
CC = load_DGD_colors(2);
acircle = 72; color_diagonal = [1 1 1]*0.5;

for ax = 1 : nAreas
    Area = Areas{ax};
    for i = 1 : nSFs
        SF = SFs{i};
        s  = find(strcmp({Rall.(Area).SF}, SF));
        
        F_DGD_max_bothDG_DGD = Rall.(Area)(s).F_DGD_max_bothDG_DGD;
        F_DGD_min_bothDG_DGD = Rall.(Area)(s).F_DGD_min_bothDG_DGD;
        F_DG_maxLR_bothDG_DGD = Rall.(Area)(s).F_DG_maxLR_bothDG_DGD;
        F_DG_sumLR_bothDG_DGD = Rall.(Area)(s).F_DG_sumLR_bothDG_DGD;
        
        figure; hs=[];
            cbar_range = [0 1];
            axes_range_1 = [-2.5/30 2.5]; axes_range_2 = [-1/30 1];

        hold on
        scatter(F_DG_sumLR_bothDG_DGD, F_DGD_max_bothDG_DGD, acircle, CC.(Area){i}, 'filled');
        h = plot([0 10],[0 10],'LineStyle','--','Color',color_diagonal);
        uistack(h, 'down');
        xlim(axes_range_1); ylim(axes_range_1);
        set(gca, 'XTick',[0:10], 'YTick',[0:10])
        xlabel('\DeltaF/F_{Monoc Left+Right}'); %xlabel('F DG sumLR (\DeltaF/F)')
        ylabel('\DeltaF/F_{Binoc Peak}'); %ylabel('F DGD peak (\DeltaF/F)')
        axis square
            if SaveFigures
%                 title(''); xlabel(''); ylabel(''); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
                title(''); axis_pos = get(gca, 'Position');
                filename = ['I:\Alessandro La Chioma\AnalyzedData\Alessandro\DGD\2017-01-06\Figures Poster\'...
                                'F_DGD_peak vs F_DG ' 'Area ' Area ' SF ' SF '.eps'];
                set(gca, 'Color', 'none'); set(gcf, 'Color', 'none');
                set(gca, 'Position',axis_pos);
                export_fig(gcf, filename, '-eps', '-q101', '-transparent')
            end
        %
        figure;
%         hs(2) = subaxis(1,2,2);
        hold on
        scatter(F_DG_maxLR_bothDG_DGD, F_DGD_min_bothDG_DGD, acircle-10, CC.(Area){i},'LineWidth',1);
%         scatter(F_DG_sumLR_bothDG_DGD, F_DGD_min_bothDG_DGD, acircle-10, CC.(Area){i},'LineWidth',1);
        h = plot([0 10],[0 10],'LineStyle','--','Color',color_diagonal);
        uistack(h, 'down');
        xlim(axes_range_2); ylim(axes_range_2);
        set(gca, 'XTick',[0:0.5:10], 'YTick',[0:0.5:10])
        xlabel('\DeltaF/F_{Monoc Max}'); %xlabel('F DG maxLR (\DeltaF/F)')
        ylabel('\DeltaF/F_{Binoc Null}'); %ylabel('F DGD min (\DeltaF/F)')
        axis square

    end
end




%% Neighbor analysis (Fig. 4B)
Areas = {'V1'};%, 'LM', 'RL'};
SFs = {'0.01'};
thr_type = 1; thr_std  = 4;
thr_di   = 0.3;
plot_figs = [];
Edges = [0,10,20,40:40:200,300,400,Inf];
PhasesToExclude = [ ];
nAreas = numel(Areas);
nSFs = numel(SFs);
NNexp=[]; NNall=[]; 

for i = 1 : nSFs
    SF = SFs{i};
    NNexp(i).SF = SF; NNall(i).SF = SF;
    for ax = 1 : nAreas
        Area = Areas{ax};

        filename = [pwd '\Vars\' 'var_DGD_' Area '.mat'];
        aDGDx = load( filename );
        aDGDx = aDGDx.('aDGD');
        
        [NNexp(i).(Area)] = Map_IOPhase_NearestNeighbor( aDGDx, [90,270], thr_type, thr_std, thr_di, SF, Edges, PhasesToExclude, plot_figs);
    end
end


%% Occurrences of disparity preference (Fig. 4C)
% You just need to load variable TC_A
% load([pwd '\Vars\TC_A__DGD_V1_LM_RL.mat']);

 
Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
SFs   = {'0.01'};
nSFs = numel(SFs);
CC = load_DGD_colors(1);
IOPhases = [0:45:360]; xrange = [-30 390];
yrange = [0 50; 0 50; 0 120]; dy = [10, 10, 20]; 

for ax = 1 : nAreas
    Area = Areas{ax};
    for i = 1 : nSFs
        SF = SFs{i};
        s = find( strcmp(TC_A.SF, SF) );
        
        CountPrefIOPhase_norm_meanAcrossExps = TC_A.(Area)(s).DD_Centered.CountPrefIOPhase_norm_meanAcrossExps;
        CountPrefIOPhase_norm_semAcrossExps  = TC_A.(Area)(s).DD_Centered.CountPrefIOPhase_norm_semAcrossExps;
        % make distribution symmetric, adding x=360=0 at the end:
        CountPrefIOPhase_norm_meanAcrossExps = 100*[CountPrefIOPhase_norm_meanAcrossExps, CountPrefIOPhase_norm_meanAcrossExps(1)];
        CountPrefIOPhase_norm_semAcrossExps  = 100*[CountPrefIOPhase_norm_semAcrossExps , CountPrefIOPhase_norm_semAcrossExps(1) ];

        hdi = figure;
        superbar( IOPhases, CountPrefIOPhase_norm_meanAcrossExps,...
            'E', CountPrefIOPhase_norm_semAcrossExps,...
            'BarWidth',0.85*45, 'BarFaceColor', CC.(Areas{ax}){s}, 'BarEdgeColor','none',...
            'ErrorbarRelativeWidth',0.0, 'ErrorbarLineWidth',2);
        
        set(gca, 'XTick',IOPhases,'XTickLabel',IOPhases,'XLim',xrange)
        set(gca, 'YTick',[0:dy(s):100], 'YLim',yrange(s,:));
        ylabel('Proportion of cells (%)' );
        xlabel('Disparity preference (phase deg)');
        title({ ['Area ' Area ' - SF ' SF ' cpd'] })

    end
end



%% Noise correlations (distribution) (Fig. 5A)
% You just need to load variables NC
% load([pwd '\Vars\NC_DGD_V1_LM_RL.mat']);

        NC = NC00;
        
Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
CC = load_DGD_colors(1); c = 2;
xrange = [-0.5 0.5]; dx = 0.25;

for ax = 1:nAreas
    Area = Areas{ax};
        
    NCv_allExpsPooled = NC.(Area).NCv_allExpsPooled;

    [NC_bincounts, edges] = histcounts(NCv_allExpsPooled, [-0.5:0.01:0.5]);

    figure; % #1
        fig_pos = get(gcf, 'position'); % fig_pos(3) = fig_pos(3)*1.5;
        set(gcf, 'position',fig_pos); % [left, bottom, width, height]
        ddx = (edges(2)-edges(1));
        xticks = edges(2:end)-ddx/2;
    superbar( xticks, NC_bincounts,...
        'BarWidth',1*ddx, 'BarFaceColor', CC.(Areas{ax}){i}, 'BarEdgeColor','none');

    title({['nPairs=' num2str(numel(NCv_allExpsPooled))];...
           ['MeanAll=' num2str(mean(NCv_allExpsPooled)) ' - ' 'Mean=' num2str(NC.(Area).NC_meanAcrossExps,'%4.3f') ]...
           });
    xlabel('Pairwise noise correlation')
    ylabel('# pairs')
    set(gca,'XLim',xrange,'XTick',[xrange(1):dx:xrange(2)],'XTickLabel',num2cell([xrange(1):dx:xrange(2)]))

end


%% Noise correlations vs. Distance vs. disparity preference difference (Fig. 5B)
% You just need to load variables NC
% load([pwd '\Vars\NC_DGD_V1_LM_RL.mat']);

        NC = NC03;

Areas = {'V1','LM','RL'};
nAreas = numel(Areas);

nBinsDist = NC.(Areas{1})(1).nBinsDist; nBinsIOP = NC.(Areas{1})(1).nBinsIOP;
EdgesDist = NC.(Areas{1})(1).EdgesDist; xvalsEdges = EdgesDist(1:end-1) ;
dist_limit = 250; dix = find(EdgesDist>dist_limit, 1) -1;

Palettes = {'Blues','Greens','Oranges'};
for ax = 1 : nAreas
    colors = cbrewer2(Palettes{ax}, 64);
    colors = flipud(colors); % puts red on top, blue at the bottom
    dc = 0.02;
    colors = colors - repmat([dc dc dc],64,1); colors(colors<0)=0;
    CCmap{ax} = colors;
end

for ax = 1:nAreas
    Area = Areas{ax};

    NC_vs_Dist_vs_IOP_allPairsPooled = NC.(Area).NC_vs_Dist_vs_IOP_allPairsPooled;

    figure; % #6
    hold on;
    imagesc( NC_vs_Dist_vs_IOP_allPairsPooled((1:dix),:) );
    cbar_h = colorbar;
    cbar_h.TickDirection = 'in'; cbar_h.Box = 'off';
    ylabel(cbar_h,'Noise correlation','rotation',270, 'VerticalAlignment','bottom', 'FontSize',14  )
%         colormap('hot');
    colormap(CCmap{ax});
    set(gca,'YDir','normal')
    set(gca,'XTick',[1:nBinsIOP], 'XTickLabel',[0:45:180])
    axis tight
    axis square;
    set(gca, 'YTick',[0:dix]+0.5, 'YTickLabel',xvalsEdges(1:dix+1))
    xlabel('Preferred disparity difference (deg)');
    ylabel(['Intra-pair distance (' char(181) 'm)'])
    set(gca,'LineWidth',2)
    box on 

end

%% Noise Correlations, perform shuffling (in preparation for Fig. 5C,D)

% load([pwd '\Vars\NC_DGD_V1_LM_RL.mat']);

        NC = NC03;

nShuffles = 1000;
Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
for ax = 1 : nAreas
    Area = Areas{ax};
    NCsh.(Area) = NoiseCorr_DGD_Shuffling(NC.(Area), nShuffles);
end

%% Noise Correlations vs. cortical distance (w/ shuffles) (Fig. 5C)
% Before this section, run the section:
% %% Noise Correlations, perform shuffling

Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
or = 2;
%
EdgesDist = NC.(Areas{1})(1).EdgesDist; dx = diff(EdgesDist(1:2)); xvalsEdges = EdgesDist(1:end-1) +dx/2 ;
dist_limit = 200; dix = find(EdgesDist>dist_limit, 1);
linestyle = '-'; markersize = 12; markerstyle = 'o';
errorbar_width = 0;
yrange = [-0.01 0.05]; dy = 0.01;

for ax = 1:nAreas
    Area = Areas{ax};

    MarkerFaceColor = CC.(Area){i};

    NCmean_Dist_perBin_acrossExps = NCsh.(Area).NCmean_Dist_perBin_acrossExps(:,or);
    NCsem_Dist_perBin_acrossExps  = NCsh.(Area).NCsem_Dist_perBin_acrossExps(:,or);

    figure; hold on;

    h = plot( xvalsEdges(1:dix), NCmean_Dist_perBin_acrossExps(1:dix),...
        'LineWidth',3, 'LineStyle',linestyle, 'Color',CC.(Area){i},...
        'Marker',markerstyle, 'MarkerFaceColor',MarkerFaceColor,'MarkerSize',markersize);
    supererr( xvalsEdges(1:dix), NCmean_Dist_perBin_acrossExps(1:dix), [], NCsem_Dist_perBin_acrossExps(1:dix),...
        errorbar_width, 'LineWidth',2, 'Color',CC.(Area){i});
    plot( xvalsEdges(1:dix), NCsh.(Area).NCmean_Dist_perBin_eachExp_shuf.prctile975(1:dix,or),...
        'LineWidth',1, 'LineStyle',':', 'Color',CC.(Area){i},...
        'Marker','none', 'MarkerFaceColor',MarkerFaceColor,'MarkerSize',markersize/2);
    plot( xvalsEdges(1:dix), NCsh.(Area).NCmean_Dist_perBin_eachExp_shuf.prctile025(1:dix,or),...
        'LineWidth',1, 'LineStyle',':', 'Color',CC.(Area){i},...
        'Marker','none', 'MarkerFaceColor',MarkerFaceColor,'MarkerSize',markersize/2);
    uistack(h, 'top');

    set(gca,'XTick',EdgesDist, 'XLim',[0 EdgesDist(dix)+dx]);
    set(gca,'YLim',yrange);
%         set(gca,'YTick',[-1:dy:1]), 'YLim',yrange);
    xlabel(['Intra pair distance (' char(181) 'm)'])
    ylabel({['Pairwise noise correlation']})

end


%% Noise Correlations vs. disparity preference difference (w/ shuffles) (Fig. 5D)
% Before this section, run the section:
% %% Noise Correlations, perform shuffling

        NC = NC03;
        
Areas = {'V1','LM','RL'};
nAreas = numel(Areas);
CC = load_DGD_colors(1);
or = 2;
linestyle = {'-'}; markersize = 12; markerstyle = {'o'};
errorbar_width = 0;
yrange = [-0.0101 0.05]; dy = 0.01;

for ax = 1 : nAreas
    Area = Areas{ax};
        
    NCmean_IOPdiff_perBin_acrossExps = NC.(Area).NCmean_IOPdiff_perBin_acrossExps;
    NCsem_IOPdiff_perBin_acrossExps  = NC.(Area).NCsem_IOPdiff_perBin_acrossExps;

    figure; hold on;
    MarkerFaceColor = CC.(Area){i};

        h= plot( [0:45:180], NCmean_IOPdiff_perBin_acrossExps(:,or),...
            'LineWidth',3, 'LineStyle',linestyle, 'Color',CC.(Area){i},...
            'Marker',markerstyle, 'MarkerFaceColor',MarkerFaceColor,'MarkerSize',markersize);
        supererr( [0:45:180], NCmean_IOPdiff_perBin_acrossExps(:,or), [], NCsem_IOPdiff_perBin_acrossExps(:,or),...
            errorbar_width, 'LineWidth',1.5, 'Color',CC.(Area){i});

        plot( [0:45:180], NCsh.(Area).NCmean_IOPdiff_perBin_perOri_eachExp_shuf.prctile975(:,or),...
            'LineWidth',1, 'LineStyle',':', 'Color',CC.(Area){i},...
            'Marker','none', 'MarkerFaceColor',MarkerFaceColor,'MarkerSize',markersize/2);
        plot( [0:45:180], NCsh.(Area).NCmean_IOPdiff_perBin_perOri_eachExp_shuf.prctile025(:,or),...
            'LineWidth',1, 'LineStyle',':', 'Color',CC.(Area){i},...
            'Marker','none', 'MarkerFaceColor',MarkerFaceColor,'MarkerSize',markersize/2);

        uistack(h, 'top');

    set(gca,'XTick',[0:45:180], 'XLim',[-20 200]);
    set(gca,'YTick',[-1:dy:1], 'YLim',yrange);
    xlabel('Preferred disparity difference')
    ylabel({['Pairwise noise correlation']})

end


%% Population decoding (SVM)
% Run this section to perform SVM decoding from scratch (this will take
% many hours, with the default settings). 
% Alternativerly, run directly the next section, which loads and plots 
% SVM classifications already computed.

% Areas = {'V1', 'LM', 'RL'};
Areas = {'V1'};
nAreas = numel(Areas);
SF = '0.01';

thr_type = 1;
thr_std  = 4;
nPredictors = [2, 5, 10, 15, 20, 40];
nCombosPred = 20;
% Directions = [2,4]'; % treat each direction (90 and 270) separately
Directions = [2,4]; % treat directions (90 and 270) together
ExpIDs = [];
Split = 6;
nShuffles = 100;
NoiseCorrelationBlind = 0;
nCVPpermutations = 1;
CombosStims = [1:8];
Method = 1;

newdir_name = [pwd '\Vars\' 'SVM_discrimination_new\' 'mc_' datestr(now,'yyyy-mm-dd_HH-MM-SS')];
mkdir(newdir_name);


for ax = 1 : nAreas
    Area = Areas{ax};

    filename = [pwd '\Vars\' 'var_DGD_' Area '.mat'];
    aDGDx = load( filename );
    aDGDx = aDGDx.('aDGD');
    Function_SVM_DGD(aDGDx, SF, thr_type, thr_std, ExpIDs,...
            nPredictors, nCombosPred, Directions, Split,...
            nShuffles, NoiseCorrelationBlind, nCVPpermutations, CombosStims, Method, newdir_name)
end


%% Population decoding (SVM), plot (Fig. 6)

Areas = {'V1', 'LM', 'RL'};

folder = [pwd '\Vars\SVM_discrimination'];
nAreas = numel(Areas);
SF = '0.01';
for ax = 1 : nAreas
    Area = Areas{ax};
    S = plot_SVM_DGD_overExps( Area, SF, folder );
end



%% Plot RDC tuning curves (as in Fig. 7A)

% load([pwd '\Vars\Fit__AsymGauss__aRDS__V1_LM_RL.mat')



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
% load([pwd '\Vars\Fit__AsymGauss__aRDS__V1_LM_RL.mat')

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

% load([pwd '\Vars\Fit__Gabor__aRDS__aRDS2__aRDS_aRDS2__V1_LM_RL.mat')

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
