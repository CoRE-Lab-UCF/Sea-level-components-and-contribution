function [OUT] = InterMontoHour(MSL)


% PROGRAM "InterMontoHour"
% Tool to interpolate monthly mean sea level(MSL) to hourly MSL.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. MSL: monthly MSL. A matrix with three columns: [year,month,MSL]  
%       
% Output: 
%       1. OUT: hourly MSL. A matrix with six columns: [year, month, day, hour, time, values]

load time
TEMP = time;
TEMP(:,5) = TEMP(:,1) + TEMP(:,2)/100;
TEMP(:,6) = 99999;
MSL(:,4) = MSL(:,1) + MSL(:,2)/100;
for j = 1:length(MSL)

    BOOL = find(TEMP(:,5)==MSL(j,4));
    TEMP(BOOL,6) = MSL(j,3);

end

BOOL = find(TEMP(:,6)==99999);
TEMP(BOOL,:) = [];

OUT = TEMP;








