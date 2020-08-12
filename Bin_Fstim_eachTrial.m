function aDGD = Bin_Fstim_eachTrial(aDGD, Param, Split)

% Split = 6;
if Split > 1
    nRois_init = length(aDGD.ROIs);
    frames_prestim    = Param.frames_prestim;
    frames_stim       = Param.frames_stim;
    addpoststimframes = Param.addpoststimframes;
    skipframes        = Param.skipframes;
    fr_start = frames_prestim+skipframes+1;
    fr_end   = frames_prestim+frames_stim+addpoststimframes;
end

switch Split
    case [0,1]
        [aDGD.ROIs.mean_Fstim_eachTrial2] = aDGD.ROIs.mean_Fstim_eachTrial;
    case 2
        range1 = [fr_start : fr_start + (fr_end-fr_start)/2];
        range2 = [range1(end)+1 : fr_end];
        for i = 1 : length(aDGD.ROIs)
            for d = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 1)
                for p = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 2)
                    rr = 0;
                    for r = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 3)
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range1 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range2 ,d,p,r )  )  ;
                    end
                end
            end
        end
%         for i = 1 : length(aDGD.ROIs)
%             for d = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 1)
%                 for p = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 2)
%                     rr = 0;
%                     for r = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 3)
%                         rr = rr + 1;
%                         aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range1 ,d,p,r )  )  ;
%                         aDGD.ROIs(nRois_init+i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range2 ,d,p,r )  )  ;
%                     end
%                 end
%             end
%         end
    case 3
        fr_int = (fr_end - fr_start)/3;
        range1 = [fr_start : fr_start + fr_int];
        range2 = [range1(end)+1 : range1(end) + fr_int];
        range3 = [range2(end)+1 : fr_end];
        for i = 1 : length(aDGD.ROIs)
            for d = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 1)
                for p = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 2)
                    rr = 0;
                    for r = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 3)
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range1 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range2 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range3 ,d,p,r )  )  ;
                    end
                end
            end
        end
    case 4
        fr_int = (fr_end - fr_start)/4;
        range1 = [fr_start : fr_start + fr_int];
        range2 = [range1(end)+1 : range1(end)+1 + fr_int];
        range3 = [range2(end)+1 : range2(end)+1 + fr_int];
        range4 = [range3(end)+1 : fr_end];
        for i = 1 : length(aDGD.ROIs)
            for d = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 1)
                for p = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 2)
                    rr = 0;
                    for r = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 3)
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range1 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range2 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range3 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range4 ,d,p,r )  )  ;
                    end
                end
            end
        end
    case 5
        fr_int = (fr_end - fr_start)/4;
        range1 = [fr_start : fr_start + fr_int];
        range2 = [range1(end)+1 : range1(end) + fr_int];
        range3 = [range2(end)+1 : range2(end) + fr_int];
        range4 = [range3(end)+1 : fr_end];
        range5 = [fr_end + 1    : fr_end + fr_int];
        for i = 1 : length(aDGD.ROIs)
            for d = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 1)
                for p = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 2)
                    rr = 0;
                    for r = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 3)
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range1 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range2 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range3 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range4 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range5 ,d,p,r )  )  ;
                    end
                end
            end
        end
    case 6
        fr_int = (fr_end - fr_start)/4;
        range1 = [fr_start : fr_start + fr_int];
        range2 = [range1(end)+1 : range1(end) + fr_int];
        range3 = [range2(end)+1 : range2(end) + fr_int];
        range4 = [range3(end)+1 : fr_end];
        range5 = [fr_end + 1    : fr_end + fr_int];
        range6 = [range5(end)+1 : range5(end)+1 + fr_int];
        for i = 1 : length(aDGD.ROIs)
            for d = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 1)
                for p = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 2)
                    rr = 0;
                    for r = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 3)
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range1 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range2 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range3 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range4 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range5 ,d,p,r )  )  ;
                        rr = rr + 1;
                        aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range6 ,d,p,r )  )  ;
                    end
                end
            end
        end
    case 9
        fr_int = 3;
        for ix = 1:Split
            range{ix} = [fr_start + (ix-1)*fr_int : fr_start + ix*fr_int-1];
        end
        for i = 1 : length(aDGD.ROIs)
            for d = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 1)
                for p = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 2)
                    rr = 0;
                    for r = 1 : size(aDGD.ROIs(1).mean_Fstim_eachTrial, 3)
                        for ra = 1 : length(range)
                            rr = rr + 1;
                            aDGD.ROIs(i).mean_Fstim_eachTrial2(d,p,rr) = nanmean(  aDGD.ROIs(i).dFoF( range{ra} ,d,p,r )  )  ;
                        end
                    end
                end
            end
        end
end