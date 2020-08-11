
function Ori = GetDI_vs_Ori( Ori, Angles, SF, aDGDtmp, thr_type, thr_std, i)

for ts = 1:length(thr_std)
    Ori(i).SF = SF;
    SpatFreq = str2double(SF);
    Ori(i).Angles = Angles;
    Ori(i).thr_type = thr_type;
    Ori(i).thr_std(ts) = thr_std(ts);
    resp_std = above_stdThreshold_DGD(aDGDtmp, thr_std(ts), thr_type);
%     Ori(i).ORoiNrs = find([aDGDtmp.ROIs.DirMaxAng] == Ori(i).Angles(1) | [aDGDtmp.ROIs.DirMaxAng] == Ori(i).Angles(2));
    Ori(i).ORoiNrs = find(ismember([aDGDtmp.ROIs.DirMaxAng], Angles));
    if isfield(aDGDtmp.ROIs,'SpatFreq')
        Ori(i).SFRoiNrs = find([aDGDtmp.ROIs.SpatFreq] == SpatFreq);
    elseif aDGDtmp.StimSettings.spacfreq == SpatFreq
        Ori(i).SFRoiNrs = [1:length(aDGDtmp.ROIs)];
    end
    Ori(i).OSFRoiNrs = intersect(Ori(i).ORoiNrs, Ori(i).SFRoiNrs);
    Ori(i).OSFrespRoiNrs{ts} = intersect(resp_std, Ori(i).OSFRoiNrs);
    if isfield(aDGDtmp.ROIs, 'ExpID') % single non-appended exps do not have ExpID
    Ori(i).nExps  = numel(unique([aDGDtmp.ROIs(Ori(i).SFRoiNrs).ExpID]));
    end
    show_figures = 0 ;
    [Ori(i).StatsDI{ts}, Ori(i).CD{ts}, Ori(i).StatsDI_perExp{ts}] =...
        PlotDI_appended(aDGDtmp , Ori(i).OSFrespRoiNrs{ts}, SpatFreq, thr_type,thr_std(ts), show_figures);
end