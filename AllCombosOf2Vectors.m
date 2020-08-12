function d = AllCombosOf2Vectors(a,b)
% a = 1:3;
% d = AllCombosOf2Vectors(a,a)
% 
% >> d =
% 
%      1     1
%      2     1
%      3     1
%      1     2
%      2     2
%      3     2
%      1     3
%      2     3
%      3     3
%
% Look also at function combnk:
% c = combnk(1:4,2)
% >> c =
%    3   4
%    2   4
%    2   3
%    1   4
%    1   3
%    1   2
     
[A,B] = meshgrid(a,b);
c = cat(2,A',B');
d = reshape(c,[],2);

