function [ang] = GetPhase(a,b)

% PROGRAM "GetPhase"
% Tool to calculate phase of semi-annual or annual cycle.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. a: estimated parameter 'a' from LS method 
%       2. b: estimated parameter 'b' from LS method 
%       
% Output: 
%       1. phase, unit radian

ang = atan(a/b);
c = a;
if c>0 && b >0
    ang = ang;
elseif c>0 && b<0
    ang = ang + pi;
elseif c<0 && b<0
    ang = ang + pi;
elseif c<0 && b>0
    ang = ang + 2*pi;
end


