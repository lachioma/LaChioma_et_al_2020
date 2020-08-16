function [StripePositions_deg, StripeStep_deg, StripeWidth_deg] =...
    GetPositionsRPS(RPS, Dir, Info, Param)


% Either Info or Param must be provided
if nargin<4 || isempty(Param)
%     Info = aRPS.Info;
    Param = LoadParam_fromInfo(Info);
end

% mouse's eyes position on the screen in [azimuth elevation], from the
% [left top] border of the screen:
mouse_center = [20 15]; % (cm) 

if strcmp(Dir, 'Elevation')
    dd = 2; %index used for RPS.positions(dd) to get figure size and in Plot_Traces
    cc = 'Y';
    if RPS.positions(1)==0
        d = 1;
    else
        d = 2;
    end
elseif strcmp(Dir, 'Azimuth')
    dd = 1; %index used for RPS.positions(dd) to get figure size and in Plot_Traces
    d = 1;
    cc = 'X';
end
stripeWidth_px  = Param.screenRes(dd) / RPS.tot_positions(dd);
StripeWidth_deg = stripeWidth_px / Param.pixperdeg(1);
mouse_dist_cm    = Param.mouse_dist_cm(1);
mouse_offset_deg = round(atan(mouse_center(dd)/mouse_dist_cm)*180/pi);

% StripePositions_deg = [0 : stripeWidth_deg/(RPS.OverlapStripes+1) : stripeWidth_deg/(RPS.OverlapStripes+1)*(RPS.positions(dd)-1) ];
% StripePositions_tot_deg = [0 : stripeWidth_deg : stripeWidth_deg*(RPS.tot_positions(dd)-1) ]';
StripePositions_tot_deg = [-mouse_offset_deg : StripeWidth_deg : StripeWidth_deg*(RPS.tot_positions(dd)-1) ]';
marg1_ix = 1+RPS.margin.(cc)(1); marg2_ix = RPS.tot_positions(dd) - RPS.margin.(cc)(2); 
StripePositions_deg = StripePositions_tot_deg(marg1_ix : marg2_ix);
if RPS.OverlapStripes
    StripePositions_deg = sort([StripePositions_deg; StripePositions_deg(1:end-1)+StripeWidth_deg/2]);
end
% the following calculates the actual center of position of the center of the
% stripe:
StripePositions_deg =  StripePositions_deg + StripeWidth_deg/2;
StripePositions_deg = -StripePositions_deg;

StripeStep_deg = abs(StripePositions_deg(1)-StripePositions_deg(2));


