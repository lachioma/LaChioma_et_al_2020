function [q2, q1, q3] = circ_quartile(alpha, AngleUnit, AngWrap)
% It returns median (q2), 25th (q1) and 75th (q2) percentile of a sample of
% circular data.
% Alessandro La Chioma, 2019-01-29

if nargin < 2
    AngleUnit = 'deg';
end
if nargin < 3
    AngWrap = 'angwrap';
end
if ~isvector(alpha)
    error('Input variable ''alpha'' must be a vector');
else
    alpha = reshape(alpha,[],1);
end

alpha = alpha(~isnan(alpha));

if strcmp(AngleUnit,'deg')
    alpha = alpha*pi/180;
end

q2 = circ_median(alpha);

dd = circ_dist(alpha, q2);

q1 = circ_median( alpha(dd <= 0) );
q3 = circ_median( alpha(dd >= 0) );

if strcmpi(AngleUnit,'deg')
%     q1 = angwrap(q1*180/pi);
%     q2 = angwrap(q2*180/pi);
%     q3 = angwrap(q3*180/pi);
    q1 = (q1*180/pi);
    q2 = (q2*180/pi);
    q3 = (q3*180/pi);
end
if strcmpi(AngWrap,'angwrap')
    q1 = angwrap(q1);
    q2 = angwrap(q2);
    q3 = angwrap(q3);
end