function C = load_DGD_colors(type)
% Output variable C is in the format
% C.(Area){i}
% with Area being {'V1','LM','RL'};


[scriptPath, scriptName] = fileparts(mfilename('fullpath')); %#ok<*ASGLU>
% load('C:\Code\Alessandro\analysis for DGD\ColorSchemes.mat');
load(fullfile(scriptPath, 'ColorSchemes.mat'));

if nargin < 1 || isempty(type) || type==1
    C = CC;
else
    C = CC2;
end

% CC=[];
% nc1 = 30; cix1 = [26,21,18]; RdBu = cbrewer2('RdBu',nc1); RdBu = mat2cell(RdBu, ones(nc1,1));
% nc2 = 60; cix2 = [52,42,37]; PiYG = cbrewer2('PiYG',nc2); PiYG = mat2cell(PiYG, ones(nc2,1));
% nc3 = 18; cix3 = [14,8,4];   Oranges = cbrewer2('Oranges',nc3);  Oranges = mat2cell(Oranges, ones(nc3,1));
% CC.V1 = RdBu(cix1);    %blue
% CC.LM = PiYG(cix2);    %green
% CC.RL = Oranges(cix3);  %red
% figure; hold on; cnt=0;
% Areas = {'V1','LM','RL'};
% for ax = 1 : numel(Areas)
%     for i = 1 : 3
%         cnt=cnt+1;
%         bar(cnt, 1, 'FaceColor', CC.(Areas{ax}){i},'LineStyle','none');
%     end
% end
% 
% CC=[];
% nc1 = 30; cix1 = [28,24,21]; RdBu = cbrewer2('RdBu',nc1); RdBu = mat2cell(RdBu, ones(nc1,1));
% nc2 = 60; cix2 = [55,45,40]; PiYG = cbrewer2('PiYG',nc2); PiYG = mat2cell(PiYG, ones(nc2,1));
% nc3 = 18; cix3 = [16,11,7];  Oranges = cbrewer2('Oranges',nc3);  Oranges = mat2cell(Oranges, ones(nc3,1));
% CC.V1 = RdBu(cix1);    %blue
% CC.LM = PiYG(cix2);    %green
% CC.RL = Oranges(cix3);  %red
% figure; hold on; cnt=0;
% Areas = {'V1','LM','RL'};
% for ax = 1 : numel(Areas)
%     for i = 1 : 3
%         cnt=cnt+1;
%         bar(cnt, 1, 'FaceColor', CC.(Areas{ax}){i},'LineStyle','none');
%     end
% end