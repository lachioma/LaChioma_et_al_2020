function A2 = angwrap(A_deg, WrapInterval)

% angwrap returns the input angles A_deg in the interval -180 < A <= 180 or
% -90 < A <= 90.
% Input A_deg can be an array and must be in degrees.
%
% Alessandro La Chioma - 2017-08-15 - created
% Alessandro La Chioma - 2017-10-22 - added A2 = mod(A2, sign(A2)*360);

if nargin < 2
    WrapInterval = 180;
end
if WrapInterval ~= 180 && WrapInterval ~= 90
    error('Wrap interval must be either 180 or 90')
end

A2 = A_deg;

% The following is to have all angles in the interval [-360 360]:
% (e.g. 400 -> 40, -400 -> -40)
A2 = mod(A2, sign(A2)*360);

A2( A2  >  WrapInterval ) = A2( A2 >   WrapInterval ) - WrapInterval*2;
A2( A2 <= -WrapInterval ) = A2( A2 <= -WrapInterval ) + WrapInterval*2;