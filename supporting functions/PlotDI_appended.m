function [StatsDI, CD, StatsDI_perExp] = PlotDI_appended(aDGD, respRois, SpatFreq, varargin)

% [scriptPath, scriptName] = fileparts(mfilename('fullpath'));

show_figures = 1 ;
% save_figures = 0 ;

%     switch save_figures case 1
%         Info = aDGD.Info;
%         fig_dir=['I:\Alessandro La Chioma\AnalyzedData\Alessandro\' Info.Mouse '\' Info.Spot '\' Info.Date...
%                   '\Figures\DisparityTuning\DI-distribution_' Info.Exp filesep];
%         if ~exist(fig_dir,'dir'), mkdir(fig_dir), end
%     end

if ~isempty(varargin)
    thr_type = varargin{1};
    thr_std = varargin{2};
    str_std = [' (thr_std=' num2str(thr_std) ',thr_type=' num2str(thr_type) ')'];
    if length(varargin) == 3
        show_figures = varargin{3};
    end
else
    str_std = '';
end
if iscell(aDGD.StimSettings)
    StimSettings = aDGD.StimSettings{1};
else
    StimSettings = aDGD.StimSettings;
end
if isfield(aDGD.ROIs,'SpatFreq')
    SFRois = find([aDGD.ROIs.SpatFreq] == SpatFreq);
elseif aDGD.StimSettings.spacfreq == SpatFreq
    SFRois = [1:length(aDGD.ROIs)];
end
% SFRois = find([aDGD.ROIs.SpatFreq] == SpatFreq);
respSFRois = intersect(respRois,SFRois);

%% Calculation pooling together ROIs from all exps

DI = [aDGD.ROIs.DI_mean]';
DIsel = DI(respSFRois);
DImean   = mean(DIsel);
DImedian = median(DIsel);

StatsDI.DI     = DIsel;
StatsDI.min    = min(DIsel);
StatsDI.max    = max(DIsel);
StatsDI.mean   = DImean;
StatsDI.median = DImedian;
StatsDI.std    = std(DIsel);

%% Calculation exp wise (only for appended exps)

if isfield(aDGD.ROIs, 'ExpID')
    
    ExpIDsall = [aDGD.ROIs.ExpID];
    ExpIDs = unique(ExpIDsall(SFRois));
    nExps  = length(ExpIDs);

    SFrespRoiNrs_eachExp = cell(nExps,1);
    DI_eachExp           = cell(nExps,1);
    DI_meanEachExp = nan(nExps,1);
    DI_stdEachExp  = nan(nExps,1);
    DI_SEMEachExp  = nan(nExps,1);

    for e = 1 : nExps

        expid = ExpIDs(e);
        ExpIDRoiNrs = find(ExpIDsall == expid);
        respSFExpIDRoiNrs_thisExp = intersect(respSFRois, ExpIDRoiNrs);
        SFrespRoiNrs_eachExp{e} = respSFExpIDRoiNrs_thisExp;

        DIsel_thisExp = DI(respSFExpIDRoiNrs_thisExp);
        DI_eachExp{e} = DIsel_thisExp;
        if ~isempty(DIsel_thisExp)
            DI_meanEachExp(e) = mean(DIsel_thisExp);
            DI_stdEachExp(e)  =  std(DIsel_thisExp);
            DI_SEMEachExp(e)  =  DI_stdEachExp(e) / sqrt(length(DIsel_thisExp));
        else continue
        end

    end

    DI_meanAcrossExps = nanmean(DI_meanEachExp);
    DI_stdAcrossExps  =  nanstd(DI_meanEachExp);
    nExps_noNaN = sum(~isnan(DI_meanEachExp));
    DI_SEMAcrossExps  = DI_stdAcrossExps / sqrt(nExps_noNaN);

    StatsDI_perExp.SFrespRoiNrs_eachExp  = SFrespRoiNrs_eachExp;
    StatsDI_perExp.DI_eachExp     = DI_eachExp;
    StatsDI_perExp.DI_meanEachExp = DI_meanEachExp;
    StatsDI_perExp.DI_stdEachExp  = DI_stdEachExp;
    StatsDI_perExp.DI_SEMEachExp  = DI_SEMEachExp;
    StatsDI_perExp.DI_meanAcrossExps = DI_meanAcrossExps;
    StatsDI_perExp.DI_stdAcrossExps  = DI_stdAcrossExps;
    StatsDI_perExp.DI_SEMAcrossExps  = DI_SEMAcrossExps;
    
else
    StatsDI_perExp = [];
end

%% Cumulative distribution

if isempty(DIsel)
    disp(['       No cells found for SF=' num2str(SpatFreq) 'cpd' str_std]);
    disp([' --->  No cumulative distribution will be calculated nor plotted for this SF']);
    CD = [];
    return
else

    [CD.f,CD.x,CD.LCB,CD.UCB] = ecdf(DIsel);


    switch show_figures case 1

        str_title = { [ 'SF = ' num2str(SpatFreq) ' cpd - TF = ' num2str(StimSettings.cyclespersecond) ' Hz'];...
            [ 'mean DI: ' num2str(DImean,'%3.2f') ' ; median DI: ' num2str(DImedian,'%3.2f') ];...
            [ 'Responsive cells: ' num2str(length(respSFRois)) '/' num2str(length(SFRois))  str_std ]};
        str_title1 = [ {['Disparity selectivity index (DI) distribution']};...
                        str_title ];
        str_title2 = [ {['Disparity selectivity index (DI) cumulative distribution']};...
                        str_title ];

        hdi=figure;
        ax1 = axes;
        histogram(ax1,DIsel, 0:1/7:1, 'Normalization','Count');
        ylabel(ax1,'# cells' );
        title( str_title1, 'Interpreter','none' );
        ax1_lim = get(ax1,'ylim');
        ax1_pos = get(ax1,'Position');
        hold on; box off;
        ax2 = axes('YAxisLocation', 'Right');
        set(ax2, 'color', 'none')
        set(ax2, 'XTick', [])
        set(ax2,'YLim',[0 ax1_lim(2)/length(respSFRois)]);
        set(ax2,'Position',ax1_pos)
        ylabel(ax2,'% cells','Rotation',270, 'VerticalAlignment','bottom')

        hdi2=figure;
        % [h,stats] = cdfplot(DIsel);
        plot(CD.x,CD.f);
        set(gca, 'XLim',[0 1]);
        xlabel('DI');
        ylabel('% cells');
        title( str_title2, 'Interpreter','none' );

    end

% switch save_figures case 1
%     figname  = [ fig_dir 'DI-Distrib_' inputname(1) '_' Info.Exp ];
%     saveas(hdi , figname, 'png'); saveas(hdi , figname ,'fig');
%     figname2 = [ fig_dir 'DI-CumulDistrib_' inputname(1) '_'  Info.Exp ];
%     saveas(hdi2, figname2,'png'); saveas(hdi2, figname2,'fig');
% %     export_fig( figname, '-m1','-png', hf );
%     disp(['Figures DI distributions for ' inputname(1) ' saved in ' fig_dir ]);
% %         close(hf)
% end

end


if exist('scriptName','var')
%     BackupRunningScript( aDGD.Info, scriptPath, scriptName )
end
