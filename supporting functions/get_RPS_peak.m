
function [StripePositions_int_deg_L, StripePositions_int_deg_R,...
          StripePositions_int_deg_LorR, StripePos_delta_deg, StripePositions_deg,...
          StripePositions_int_deg_meanLR] = ...
                get_RPS_peak( aRPS, aRPS2, Param, ListRois, Dir, Plot_RF,...
                interp_factor, smooth_factor, RFposPeak_RPS_fromFit, thr_Rsquared_rps,...
                respRPS, respRPS2)
% Get left and right RF center in degrees.
% Get also the difference StripePos_delta_deg.
% Get the center of either left or right (the one with higher response)
% StripePositions_int_deg_LorR.
% You can omit aRPS2 by leaving empty [].
% RF center is given by the peak of responses to RPS.
% In this code, 0 deg corresponds to the upper margin of the monitor for
% "Elevation', and left margin of the monitor for 'Azimuth'.

% mouse's eyes position on the screen in [azimuth elevation], from the
% [left top] border of the screen:
% mouse_center = [0 15]; % (cm) 

RPS = aRPS.StimSettings;
StripePositions_int_deg_L      = nan(numel(ListRois),1);
StripePositions_int_deg_R      = nan(numel(ListRois),1);
StripePositions_int_deg_LorR   = nan(numel(ListRois),1);
StripePositions_int_deg_meanLR = nan(numel(ListRois),1);
StripePos_delta_deg            = nan(numel(ListRois),1);

if nargin < 6 || isempty(Plot_RF)
    Plot_RF = 0;
end
if nargin < 7 || isempty(interp_factor)
    interp_factor = 1; %3
end
if nargin < 8 || isempty(smooth_factor)
    smooth_factor = 1; %5
end
if nargin < 9 || isempty(RFposPeak_RPS_fromFit)
    RFposPeak_RPS_fromFit = 0;  
end
% if nargin < 10 || isempty(thr_Rsquared_rps)
% end
if nargin < 11% || isempty(respRPS )
    respRPS  = [];
end
if nargin < 12% || isempty(respRPS2)
    respRPS2 = [];
end

[StripePositions_deg, StripeStep_deg, StripeWidth_deg] =...
    GetPositionsRPS(RPS, Dir, aRPS.Info, Param); %#ok<ASGLU>



if strcmp(Dir, 'Elevation')
    dd = 2; %index used for RPS.positions(dd) to get figure size and in Plot_Traces
    if RPS.positions(1)==0
        d = 1;
    else
        d = 2;
    end
elseif strcmp(Dir, 'Azimuth')
    dd = 1; %index used for RPS.positions(dd) to get figure size and in Plot_Traces
    d = 1;
end


if RFposPeak_RPS_fromFit == 0
    X  = [1 : 1 : RPS.positions(dd)];
    % Xi = [1 : (RPS.positions(dd)-1)/(StripePositions_deg(end)-1) : RPS.positions(dd)];
    Xi = linspace(1, RPS.positions(dd), RPS.positions(dd)*interp_factor );
    StripePositions_int_deg = linspace(StripePositions_deg(1),StripePositions_deg(end), RPS.positions(dd)*interp_factor)';
end

for i = 1 : numel(ListRois)
    
    nr = ListRois(i);
    

    % Left eye:
    if RFposPeak_RPS_fromFit
        if ~isempty(aRPS) && ismember(nr,respRPS )  &&...
            isfield(aRPS.TC(nr).GoF,'Rsquared_meanTC') && aRPS.TC(nr).GoF.Rsquared_meanTC > thr_Rsquared_rps
                StripePositions_int_deg_L(i) = aRPS.TC(nr).RFposPeak;
                F_L_int = reshape(aRPS .TC(nr).FittedData,[],1);
                flag(1) = 1;
        else
                StripePositions_int_deg_L(i) = nan;
                flag(1) = 0;
        end
        
        if ~isempty(aRPS2) && ismember(nr,respRPS2) &&...
            isfield(aRPS2.TC(nr).GoF,'Rsquared_meanTC') && aRPS2.TC(nr).GoF.Rsquared_meanTC > thr_Rsquared_rps
                StripePositions_int_deg_R(i) = aRPS2.TC(nr).RFposPeak;
                F_R_int = reshape(aRPS2.TC(nr).FittedData,[],1);
                flag(2) = 1;
            else
                StripePositions_int_deg_R(i) = nan;
                flag(2) = 0;
        end
        
            
        if sum(flag) > 0 % fit RF exists for either L and R
            if sum(flag)==2 % for both L and R
                F_LorR_int   = max([F_L_int, F_R_int],[],2);
                [~, Ix_LorR] = max(F_LorR_int);
                [~, Ix_L   ] = max(F_L_int);
                [~, Ix_R   ] = max(F_R_int);
                Ix_meanLR = round(mean([Ix_L, Ix_R]));
            elseif sum(flag)==1 % for only L or R
                if flag(1)
                    [~, Ix_LorR] = max(F_L_int);
                    Ix_meanLR = Ix_LorR;
                elseif flag(2)
                    [~, Ix_LorR] = max(F_R_int);
                    Ix_meanLR = Ix_LorR;
                end
            end
            StripePositions_int_deg = aRPS.TC(nr).FittedData_x;
            StripePositions_int_deg_LorR(i)   = StripePositions_int_deg(Ix_LorR);
            StripePositions_int_deg_meanLR(i) = StripePositions_int_deg(Ix_meanLR);
            StripePos_delta_deg(i) = StripePositions_int_deg_L(i) - StripePositions_int_deg_R(i);
        else
            StripePositions_int_deg_LorR(i)   = nan;
            StripePositions_int_deg_meanLR(i) = nan;
            StripePos_delta_deg(i)            = nan;
        end
        
        
        
    else % if RFposPeak_RPS_fromFit==0 % take from dFoF
        
        
        
        if ~isempty(aRPS)
            F_L = nan(RPS.positions(dd),1);
        end
        if ~isempty(aRPS2)
            F_R = nan(RPS.positions(dd),1);
        end

        for p = 1 : RPS.positions(dd)
            if ~isempty(aRPS) && ismember(nr,respRPS )
                F_L(p) = aRPS .ROIs(nr).max_Fstim(d,p);
            end
            if ~isempty(aRPS2) && ismember(nr,respRPS2)
                F_R(p) = aRPS2.ROIs(nr).max_Fstim(d,p);
            end
        end
        
        if ~isempty(aRPS ) && ismember(nr,respRPS )
            F_L_int = interp1(X, F_L, Xi, 'pchip');
            F_L_int = smooth(F_L_int, smooth_factor);
        %     figure; plot(F_L);
        %     figure; plot(F_L_int);
            [~, Ix_L] = max(F_L_int);
            StripePositions_int_deg_L(i) = StripePositions_int_deg(Ix_L);
        end
        if ~isempty(aRPS2) && ismember(nr,respRPS2)
            F_R_int = interp1(X, F_R, Xi, 'pchip');
            F_R_int = smooth(F_R_int, smooth_factor);
            [~, Ix_R] = max(F_R_int);
            StripePositions_int_deg_R(i) = StripePositions_int_deg(Ix_R);
        end
        if ~isempty(aRPS ) && ismember(nr,respRPS ) &&...
           ~isempty(aRPS2) && ismember(nr,respRPS2)  
            F_LorR_int = max([F_L_int, F_R_int],[],2);
            [~, Ix_LorR] = max(F_LorR_int);
            StripePositions_int_deg_LorR(i) = StripePositions_int_deg(Ix_LorR);
            %
            StripePos_delta_deg(i) = StripePositions_int_deg_L(i) - StripePositions_int_deg_R(i);
        end
%         F_zint = zscore(F_int);
%         % Apply threshold:
%         F_zint_thr = F_zint;
%         F_zint_thr( F_zint_thr < thr_zscore ) = 0;
        switch Plot_RF, case 1
            figure; hold on;
            plot( StripePositions_int_deg, F_L_int );
            plot( StripePositions_int_deg, F_R_int );
            title([ 'ROI # ' num2str(nr) ]);
        end
    end   
end

% %     if RFposPeak_RPS_fromFit
% %         if isfield(aRPS.TC(nr).GoF,'Rsquared_meanTC') && aRPS.TC(nr).GoF.Rsquared_meanTC > thr_Rsquared_rps
% %             StripePositions_int_deg_L(i) = aRPS.TC(nr).RFposPeak;
% %             flag(1) = 1;
% %         else
% %             StripePositions_int_deg_L(i) = nan;
% %             flag(1) = 0;
% %         end
% %         if ~isempty(aRPS2) && ismember(nr,respRPS2)
% %             if isfield(aRPS2.TC(nr).GoF,'Rsquared_meanTC') && aRPS2.TC(nr).GoF.Rsquared_meanTC > thr_Rsquared_rps
% %                 StripePositions_int_deg_R(i) = aRPS2.TC(nr).RFposPeak;
% %                 flag(2) = 1;
% %             else
% %                 StripePositions_int_deg_R(i) = nan;
% %                 flag(2) = 0;
% %             end
% %             %
% %             F_L_int = reshape(aRPS .TC(nr).FittedData,[],1);
% %             F_R_int = reshape(aRPS2.TC(nr).FittedData,[],1);
% %             if sum(flag) > 0
% %                 if sum(flag)==2
% %                     F_LorR_int = max([F_L_int, F_R_int],[],2);
% %                     [~, Ix_LorR] = max(F_LorR_int);
% %                     [~, Ix_L] = max(F_L_int);
% %                     [~, Ix_R] = max(F_R_int);
% %                     Ix_meanLR = round(mean([Ix_L, Ix_R]), 1); %round to 0.1 to match x of fitted data
% %                 elseif sum(flag)==1
% %                     if flag(1)
% %                         [~, Ix_LorR] = max(F_L_int);
% %                     elseif flag(2)
% %                         [~, Ix_LorR] = max(F_R_int);
% %                     end
% %                 end
% %                 StripePositions_int_deg = aRPS.TC(nr).FittedData_x;
% %                 StripePositions_int_deg_LorR(i) = StripePositions_int_deg(Ix_LorR);
% %                 StripePositions_int_deg_meanLR(i) = StripePositions_int_deg(Ix_meanLR);
% %                 StripePos_delta_deg(i) = StripePositions_int_deg_L(i) - StripePositions_int_deg_R(i);
% %             else
% %                 StripePositions_int_deg_LorR(i)   = nan;
% %                 StripePositions_int_deg_meanLR(i) = nan;
% %                 StripePos_delta_deg(i)            = nan;
% %             end
% %         end
% %         
% %     else % take from dFoF
% %         
% %         F_L = nan(RPS.positions(dd),1);
% %         if ~isempty(aRPS2)
% %             F_R = nan(RPS.positions(dd),1);
% %         end
% % 
% %         for p = 1 : RPS.positions(dd)
% % 
% %             F_L(p) = aRPS .ROIs(nr).max_Fstim(d,p);
% %             if ~isempty(aRPS2)
% %                 F_R(p) = aRPS2.ROIs(nr).max_Fstim(d,p);
% %             end
% % 
% %         end
% %         
% %         F_L_int = interp1(X, F_L, Xi, 'pchip');
% %         F_L_int = smooth(F_L_int, smooth_factor);
% %     %     figure; plot(F_L);
% %     %     figure; plot(F_L_int);
% %         [~, Ix_L] = max(F_L_int);
% %         StripePositions_int_deg_L(i) = StripePositions_int_deg(Ix_L);
% %         if ~isempty(aRPS2)
% %             F_R_int = interp1(X, F_R, Xi, 'pchip');
% %             F_R_int = smooth(F_R_int, smooth_factor);
% %             [~, Ix_R] = max(F_R_int);
% %             StripePositions_int_deg_R(i) = StripePositions_int_deg(Ix_R);
% %             %
% %             F_LorR_int = max([F_L_int, F_R_int],[],2);
% %             [~, Ix_LorR] = max(F_LorR_int);
% %             StripePositions_int_deg_LorR(i) = StripePositions_int_deg(Ix_LorR);
% %             %
% %             StripePos_delta_deg(i) = StripePositions_int_deg_L(i) - StripePositions_int_deg_R(i);
% %         end
% %         
% %         
% %         
% % %         F_zint = zscore(F_int);
% % %         % Apply threshold:
% % %         F_zint_thr = F_zint;
% % %         F_zint_thr( F_zint_thr < thr_zscore ) = 0;
% %     
% %         switch Plot_RF case 1
% %             figure; hold on;
% %             plot( StripePositions_int_deg, F_L_int );
% %             plot( StripePositions_int_deg, F_R_int );
% %             title([ 'ROI # ' num2str(nr) ]);
% %         end    
% %     end
% % end


% function Param = LoadParam_fromInfo(Info) %#ok<STOUT>
%     auxdata_dir = ['E:\MyData\' Info.Mouse '\Data\' Info.Date filesep];
%     % Stimulation info:
%     % first, get the proper index to load:
%     InfoFiles = dir([auxdata_dir 'InfoExp*' '.mat']);
%     if isempty(InfoFiles)
%         auxdata_dir = ['I:\Alessandro La Chioma\RawData\alessandro\' Info.Mouse '\Data\' Info.Date filesep];
%         InfoFiles = dir([auxdata_dir 'InfoExp*' '.mat']);
%     end
%     if isempty(InfoFiles)
%         error('No InfoExp file found')
%     end
%     for f = 1:length(InfoFiles)
%         if ~isempty( strfind( InfoFiles(f).name, Info.Exp ) )
%             ind = f;
%             break
%         end
%     end
% 
% % catch
% %     auxdata_dir = ['I:\Alessandro La Chioma\RawData\alessandro\' Info.Mouse '\Data\' Info.Date filesep];
% %     load([auxdata_dir InfoStimFile(ind).name]);
% % end
% 
% InfoStimFile = dir([auxdata_dir 'InfoStim*.mat']);
% load([auxdata_dir InfoStimFile(ind).name], 'Param');
% 
% end




