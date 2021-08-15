function [SLH] = ShapeRecord2Mon(signal,time)

% PROGRAM "ShapeRecord2Mon"
% Tool to reshpae hourly water level to monthly format.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. signal: hourly water level. 
%       2. time: A matrix with four columns: [year,month,day,hour] ;
%       
% Output: 
%       1. SLH: A matrix with columns: 
%          (1). the first two columns: [year, month]
%          (2). the third column: total days in this month 
%          (3). the forth column:  total number of water level observations in this month
%          (4). fith to last column: water level observation in this month
if length(signal)~= length(time)
   error([' not equal length ']);
end

m = length(time);

monthold = 0;   cont = 0;   cont1 = 0;

for j = 1:m
%     year = fix(time(j));
%     doy = (time(j) - year)*1000-1 +0.000000001;
%     [month,day,allday] = Doy2YMD(year,doy);
    year = time(j,1);
    month = time(j,2);
    day = time(j,3);

    if month~= monthold
        cont = cont + 1;
        cont1 = 1;
        

        MON = [31 28 31 30 31 30 31 31 30 31 30 31];
        if mod(year,4)~=0
            MON(2) = 28;
        elseif mod(year,100)==0 && mod(year,400)~=0
            MON(2) = 28;
        else
            MON(2) = 29;
        end


        allday =  MON(month);
        
        SLH(cont,1:4) = [year,month,allday,allday*24] ;
        SLH(cont,4+cont1) = signal(j);
    else   
        cont1 = cont1 + 1;
        SLH(cont,4+cont1) =  signal(j);        
    end
    monthold = month;
end








