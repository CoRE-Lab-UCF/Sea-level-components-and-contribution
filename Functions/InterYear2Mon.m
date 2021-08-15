function [T_BOX] = InterYear2Mon(MSL_YEAR,MSL_MON)
  
% PROGRAM "InterYear2Mon"
% Tool to interpolate Yearly mean sea level(MSL) to monthly MSL.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. MSL_YEAR: yearly MSL. A matrix with two columns: [time,MSL]  
%       2. MSL_MON:  Time of monthly MSL, A vector with format like
%       year + month/12 - 1/24. i.e., 1991.0417 (Jan 1991 ) = 1991 + 1/12 - 1/24
%       
% Output: 
%       1. T_BOX: monthly MSL. A matrix with three columns: [year,
%       month,MSL].

TIME = MSL_YEAR(:,1);
MSL = MSL_YEAR(:,2);
FIT_MON = [];
m = length(TIME);
for j = 1:m-1
    TM = [];
    if j==1
        SLOPE = (MSL(j+1) - MSL(j));
        TEMP = [-6:11]/12*SLOPE;
        TM(1:12,1) = TIME(j);
        TM(1:12,2) = [1:12]';
        TM(13:18,1) = TIME(j+1);
        TM(13:18,2) = [1:6]';
        TM(:,3) = TEMP'+MSL(j);
        FIT_MON = [FIT_MON;TM];
        
    elseif j==m-1
        SLOPE = (MSL(j+1) - MSL(j));
        TEMP = [0:17]/12*SLOPE;
        TM(1:6,1) = TIME(j);
        TM(1:6,2) = [7:12]';
        TM(7:18,1) = TIME(j+1);
        TM(7:18,2) = [1:12]';
        TM(:,3) = TEMP'+MSL(j);
        FIT_MON = [FIT_MON;TM];
        
    else
        SLOPE = (MSL(j+1) - MSL(j));
        TEMP = [0:11]/12*SLOPE;
        TM(1:6,1) = TIME(j);
        TM(1:6,2) = [7:12]';
        TM(7:12,1) = TIME(j+1);
        TM(7:12,2) = [1:6]';
        TM(:,3) = TEMP'+MSL(j);
        FIT_MON = [FIT_MON;TM];
    end

end

TIME = [];
TIME_MON = MSL_MON(:,1);
TIME(:,1) = fix(TIME_MON);
TIME(:,2) = round((TIME_MON - TIME(:,1))*12 + 1/24);

t = TIME(:,1) + TIME(:,2)/100;

tmon = FIT_MON(:,1) + FIT_MON(:,2)/100;
[a b c] = intersect(tmon,t);

T_BOX = TIME(:,1:2);
T_BOX(c,1:3) = FIT_MON(b,:);

BOOL = (T_BOX(:,3)==0)+0;
DIFF = diff(BOOL);
st = find(DIFF==-1);
en = find(DIFF==1);
if isempty(st) && ~isempty(en) 
    SLOPE = MSL(end) - MSL(end-1);
    TEMP = [18:18+length(T_BOX)-en-1]/12*SLOPE+MSL(end);
    T_BOX(en+1:end,3) = TEMP;    
elseif ~isempty(st) && isempty(en) 
    SLOPE = MSL(2) - MSL(1);
    TEMP = [-7-st+1:-7]/12*SLOPE+MSL(1);
    T_BOX(1:st,3) = TEMP;    
elseif ~isempty(st) && ~isempty(en)     
    SLOPE = MSL(2) - MSL(1);
    TEMP = [-7-st+1:-7]/12*SLOPE+MSL(1);
    T_BOX(1:st,3) = TEMP;
    
    SLOPE = MSL(end) - MSL(end-1);
    TEMP = [18:18+length(T_BOX)-en-1]/12*SLOPE+MSL(end);
    T_BOX(en+1:end,3) = TEMP;    
   
    
end






