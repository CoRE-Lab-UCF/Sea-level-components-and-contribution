function [coef] = Write2Coef(coef,CARD,INDEX)

% PROGRAM "Write2Coef"
% Tool to exclude SA (annual cycle) and SSA (semi-annual cycle) from coef, 
% coef is a struct containing info about tidal constituents, and will be
% used to reconstruct tides. This is becuase the UTide package also estimate seasonal cycle but we do
% not need that in our analysis.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. coef:  a struct containing info on tidal constituents, and will be
%         used to reconstruct tides.
%       2. CARD: the name of tidal constituents used to reconstruct tides. 
%          CARD will be wirtten into coef   
%       3. INDEX: an switch whether consider SLR in reconstruction process or not . If NOT, use
%       'NoMSL'. Otherwise, use any number or string.
% Output: 
%       1. coef 

name = coef.name;
% CARD = {'SSA';'MSM';'MM';'MSF';'MF'}; % 
m1 = length(name);
A = zeros(m1,1);
A_ci = zeros(m1,1);
g = zeros(m1,1);
g_ci = zeros(m1,1);    
BOOL = [];
for k = 1:length(CARD)
    for q = 1:m1
        if (length(name{q})==length(CARD{k})) && sum(name{q}==CARD{k})==length(name{q})
            BOOL(k) = q;
        end     
    end 
end

if ~isempty(BOOL)
    A(BOOL) = coef.A(BOOL);
    A_ci(BOOL) = coef.A_ci(BOOL);
    g(BOOL) = coef.g(BOOL);
    g_ci(BOOL) = coef.g_ci(BOOL);
end

coef.A = A;
coef.A_ci = A_ci;
coef.g = g;
coef.g_ci = g_ci;

if INDEX == 'NoMSL'
    coef.slope = 0;
    coef.mean = 0;
end




