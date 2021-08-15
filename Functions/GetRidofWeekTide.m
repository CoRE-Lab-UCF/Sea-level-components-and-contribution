function [name] = GetRidofWeekTide(name)

% PROGRAM "GetRidofWeekTide"
% Tool to exclude SA (annual cycle) and SSA (semi-annual cycle) from tidal constituents. 
% This is becuase the UTide package also estimate seasonal cycle but we do
% not need that in our analysis.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. name: name of consituents
%            
% Output: 
%       1. name: name of constituents withour SA and SSA. 

BOOL = NaN*ones(6,1);
m = length(name);
CARD = {'SSA';'SA';}; % 'MSM';'MM';'MSF';'MF'
for k = 1:length(CARD)
    for q = 1:length(name)
        if (length(name{q})==length(CARD{k})) && sum(name{q}==CARD{k})==length(name{q})
            BOOL(k) = q;
        end     
    end 
end

s = isnan(BOOL);
BOOL(s) = [];
name(BOOL) = [];

