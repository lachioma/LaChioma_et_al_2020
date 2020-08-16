function areah = AddErrorArea(x,y, erry, ColorArea, FaceAlpha, erry_bottom, EdgeColor )
% areah = AddErrorArea(x,y, erry, ColorArea, FaceAlpha, erry_bottom, EdgeColor )
% Add area for error (e.g. std or sem or quartiles) around y (mean or median).
% For asymmetric error area, use erry_bottom.
% Adapted from http://www.mathworks.com/matlabcentral/fileexchange/13188-shade-area-between-two-curves):

if isempty(x)
    x = 1:length(y);
end
if all(erry==0)
    return
end
if ~exist('ColorArea','var') || isempty(ColorArea)
%     ColorArea = [0 0 1]*0.7 ;
    ploth = findobj(gca,'type','line');
    % Color of last plotted line (it's the line with index 1):
    ColorArea = ploth(1).Color;
end
if ~exist('FaceAlpha','var') || isempty(FaceAlpha)
    FaceAlpha = 0.3;
end
if ~exist('erry_bottom','var') || isempty(erry_bottom)
    erry_bottom = erry;
end
if ~exist('EdgeColor','var') || isempty(EdgeColor)
    EdgeColor = 'none';
end
x            = reshape(x,1,[]);
y            = reshape(y,1,[]);
erry         = reshape(erry ,1,[]);
erry_bottom  = reshape(erry_bottom ,1,[]);

if any(isnan(y))
    x = x(~isnan(y));
    erry = erry(~isnan(y));
    erry_bottom = erry_bottom(~isnan(y));
    y = y(~isnan(y));
end


% Fill area for e.g. standard error of the mean
hold on
areah = fill( [ x fliplr([ x ]) ] , ...      % x values
                [ y + erry ...  % upper part
                fliplr(y - erry_bottom ) ] ,... % lower part
               ColorArea ) ; % color of the area
% areah = patch( [ x fliplr([ x ]) ] , ...      % x values
%                 [ y + erry ...  % upper part
%                 fliplr(y - erry_bottom ) ] ,... % lower part
%                ColorArea ) ; % color of the area
% Default area settings
% set edge color and transparency (Alpha) of the area
% set(areah, 'EdgeColor','none' , 'FaceAlpha',FaceAlpha) ; 
set(areah, 'EdgeColor',EdgeColor , 'FaceAlpha',FaceAlpha) ; 
% Bring error area behind meany:
uistack(areah, 'down');

