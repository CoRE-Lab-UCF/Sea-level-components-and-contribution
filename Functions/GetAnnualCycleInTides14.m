function [HC] = GetAnnualCycleInTides14(WL,Lat,name,WinLen,CARD2,HC,Percent)

% PROGRAM "GetAnnualCycleInTides14"
% Tool to do todal harmonic analysis.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. WL: water level. A matrix with five columns: [year, month, day, hour, water level]
%       2. Lat: latitude of the gauge
%       3. name: name of the station.
%       4. WinLen: window length for tidal harmonic analysis, unit: month.
%          usually we choose '12' months.
%       5. CARD2: the name of tidal consitutents, usually we use 143 constituents.
%       6, HC: a struct save each station's results.    
%       7. Percent: a threshold. If data completeness is less that (1-Percent),  
%          we will not conduct tidal harmonic analysis in this calendar year.
% Output: 
%       1. HC: a struct save each station's results. HC contains 5 matrixs and 1 cell:
%             For 'Amp' , a matrix saves amplitudes of components over time.
%             For 'Pha' , a matrix saves phases of components over time.
%             For 'AmpErr' , a matrix saves amplitude errors of components over time.
%             For 'PhaErr' , a matrix saves phase errors of components over time.
%             For 'consname' , a cell saves constituents analyzed.
%             For 'Tide_weekly', a the reconstructed tides.
%
% Reference:  http://www.po.gso.uri.edu/~codiga/utide/utide.html
% 
T0 = WL(:,1) + WL(:,2)/100 ;
[SLH] = ShapeRecord2Mon(WL(:,5),WL(:,1:4));
SLH(1,:) = [];
% SLH(end,:) = [];
[m,n] = size(SLH);
FIT = NaN*ones(length(WL),1);
FIT2 = NaN*ones(length(WL),1);
FIT3 = NaN*ones(length(WL),length(CARD2));
YEAR = unique(WL(:,1));
YEAR(1) = [];
b0 = find(WL(:,1)==YEAR(1));
m0 = length(YEAR);
BAR_BOX = NaN*ones(m0,1);

for j = 1:length(YEAR)
    TEMP_H = [];
    TEMP_T = [];
    CONT = 0;    
    BASE = WinLen-12;
    STEP = round(BASE/2);
    BE = [];
    Start_time = [];UT_tide = [];
    BOOL0 = find(WL(:,1)==YEAR(j));
    BOOL = find(SLH(:,1)==YEAR(j));
    if (j-1)*12<STEP
        T2 = find(T0==SLH(BOOL(end) + STEP,1) +SLH(BOOL(end) + STEP,2)/100);
        TEMP_1 = [WL(b0:T2(end),1:4),zeros(length([b0:T2(end)]),2)];
        Start_time = datenum(TEMP_1);
        UT_tide = WL(b0:T2(end),5);
    elseif (j-1)*12>=STEP &&( m0-j)*12>=STEP
        T3 = find(T0==SLH(BOOL(1) - STEP,1) +SLH(BOOL(1) - STEP,2)/100);
        T2 = find(T0==SLH(BOOL(end) + STEP,1) +SLH(BOOL(end) + STEP,2)/100);
        TEMP_1 = [WL(T3(1):T2(end),1:4),zeros(length([T3(1):T2(end)]),2)];        
        Start_time = datenum(TEMP_1);
        UT_tide = [WL(T3(1):T2(end),5);];  
    elseif (m0-j)*12<STEP
        T3 = find(T0==SLH(BOOL(1) - STEP,1) +SLH(BOOL(1) - STEP,2)/100);
        TEMP_1 = WL(T3(1):end,1:4);
        Start_time = datenum([TEMP_1,zeros(length(TEMP_1),2)]);
        UT_tide = WL(T3(1):end,5);

    end
    
    latitude = Lat;
    if sum(isnan(UT_tide))<=length(UT_tide)*Percent
        [ coef ] = ut_solv(Start_time,UT_tide,[],latitude,'auto','Huber');
        HC.(name).Amp(1:length(coef.A(:,1)),j) = coef.A(:,1);
        HC.(name).Pha(1:length(coef.A(:,1)),j) = coef.g(:,1);
        HC.(name).AmpErr(1:length(coef.A(:,1)),j) = coef.A_ci(:,1);
        HC.(name).PhaErr(1:length(coef.A(:,1)),j) = coef.g_ci(:,1);
        HC.(name).consname(1:length(coef.name),j) = coef.name;          
        coef.mean = 0; coef.slope = 0;
        if j==1
            PreStart_time = datenum([WL(1:BOOL0(end),1:4) zeros(length([1:BOOL0(end)]),2)]); 
            [ sl_fit] = ut_reconstr(PreStart_time,coef); 
            FIT(1:BOOL0(end),1) = sl_fit';
            
            [CARD] = GetRidofWeekTide(coef.name);
            [COEF] = Write2Coef(coef,CARD,'NoMSL');
            [ sl_fit] = ut_reconstr(PreStart_time,COEF); 
            FIT2(1:BOOL0(end),1) = sl_fit';  
            
%             for jj = 1:length(CARD2)
%                 [COEF] = Write2Coef(coef,CARD2(jj),'NoMSL');
%                 [ sl_fit, ~,BAR ] = ut_reconstr_BAR(PreStart_time,COEF); 
%                 FIT3(1:BOOL0(end),jj) = sl_fit';         
%             end
            
        elseif j==m0
            PreStart_time = datenum([WL(BOOL0,1:4) zeros(length(BOOL0),2)]);
            [ sl_fit] = ut_reconstr(PreStart_time,coef); 
            FIT(BOOL0,1) = sl_fit';
            
            [CARD] = GetRidofWeekTide(coef.name);            
            [COEF] = Write2Coef(coef,CARD,'NoMSL');
            [ sl_fit] = ut_reconstr(PreStart_time,COEF); 
            FIT2(BOOL0,1) = sl_fit';    
                        
%             for jj = 1:length(CARD2)
%                 [COEF] = Write2Coef(coef,CARD2(jj),'NoMSL');
%                 [ sl_fit, ~,BAR ] = ut_reconstr_BAR(PreStart_time,COEF); 
%                 FIT3(BOOL0,jj) = sl_fit';         
%             end            
            
        else
            PreStart_time = datenum([WL(BOOL0,1:4) zeros(length(BOOL0),2)]);
            [ sl_fit] = ut_reconstr(PreStart_time,coef);  
            FIT(BOOL0,1) = sl_fit';
            
            [CARD] = GetRidofWeekTide(coef.name);
            [COEF] = Write2Coef(coef,CARD,'NoMSL');
            [ sl_fit] = ut_reconstr(PreStart_time,COEF); 
            FIT2(BOOL0,1) = sl_fit';               
            
%             for jj = 1:length(CARD2)
%                 [COEF] = Write2Coef(coef,CARD2(jj),'NoMSL');
%                 [ sl_fit, ~,BAR ] = ut_reconstr_BAR(PreStart_time,COEF); 
%                 FIT3(BOOL0,jj) = sl_fit';         
%             end                
            
        end    
        

    else
        HC.(name).Amp(1:80,j) = NaN;
        HC.(name).Pha(1:80,j) = NaN;
        HC.(name).AmpErr(1:80,j) = NaN;
        HC.(name).PhaErr(1:80,j) = NaN;    
        HC.(name).consname(1:80,j) = {'NaN'}; 
    end
    
    HC.(name).CenTime(j,1:2) = [SLH(j+round(WinLen/2)-1,1:2)];
    
end
% HC.(name).BAR = BAR_BOX;
% HC.(name).Tide_All = FIT;
HC.(name).Tide_weekly = FIT2;
% HC.(name).Tide = FIT - FIT2;
% HC.(name).WL = WL;
% HC.(name).Tide_Each = FIT3;




