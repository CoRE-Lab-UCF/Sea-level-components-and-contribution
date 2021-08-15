clc
clear


% load('F:\ReCalculation\Download_MonthlyMSL_NOAA\NOAA_MSL_ALL');
% 
% load('F:\ReCalculation\Download_MonthlyMSL_NOAA\NOAA_WL');
% NOAA_WL_SHORT = NOAA_WL;
% load('F:\ReCalculation\SeasonalSL\SENTIVITY_ALL');

% NAMEb = NOAA_MSL_SHORT.NAME_before1990;
% 
% MSL_M = NOAA_MSL_SHORT.(NAMEb{32}).Monthly;
% MSL_Y = NOAA_MSL_SHORT.(NAMEb{32}).Yearly;
% 
% Lat = NOAA_WL_SHORT.(NAMEb{32}).Lat;
% WL = NOAA_WL_SHORT.(NAMEb{32}).WL;
% MSL0 = NOAA_WL_SHORT.(NAMEb{32}).MSL;
% MHHW = NOAA_WL_SHORT.(NAMEb{32}).MHHW;
% NFL = NOAA_WL_SHORT.(NAMEb{32}).NFL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data preparation
load Data_San_Francisco
% MSL_M : monthly mean sea level
% MSL_Y : yearly mean sea level
% Lat : latitude at San Francisco
% WL : water level. WL is a matrix with five columns. First to froth column: [year month day hour](datum: GMT)
% The fifth column: water level observation (Datum: STND).  
% MSL0 : mean sea level. This mean sea level is the
% average of water level from 1983 to 2001.
% MHHW : mean higher high water
% NFL = nuisance flooding threshold. 

BOOL = find(MSL_M(:,1)==1950);  % To select data from 1950-now.
if ~isempty(BOOL)
    MSL_M = MSL_M(BOOL(1):end,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% design for Sea level rise
[ts_yr,BestRC]=TS_Stat_4_yr(MSL_Y(:,[1,4]),15); % Singular spectrum analysis
[FIT] = InterYear2Mon([MSL_Y(:,1),BestRC],MSL_M(:,1)+MSL_M(:,2)/12-1/24); %Interplote Yearly to monthly
BOOL = find(FIT(:,1)>2019);
FIT2 = FIT;
FIT2(BOOL,:) = [];
[NonLin_Hour] = InterMontoHour(FIT2(:,:)); % interpolate monthly to hourly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% design for Seasonal cycle
[TotalPreMSL] = GetMovingWindow_Monthly(MSL_M(:,3),MSL_M(:,1)+MSL_M(:,2)/12-1/24,5);      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% design for interannual to decadal variarity

MSL_M(:,5) = MSL_M(:,4) - FIT(:,3) - sum(TotalPreMSL.Monthly(:,3:4),2);
[RecSig_1] = WavePacketFilter(MSL_M(:,5),12,1/13,1/2,1); % wavelet package filter.
[RecSig_3] = MSL_M(:,5) - RecSig_1;
RecSig_2 = zeros(length(RecSig_1),1);
RecSig_4 = zeros(length(RecSig_1),1);
RecSig_Non = FIT(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% design for Tide
load ALL_NAME  % 143 tidal constituents
HC = [];
[HC] = GetAnnualCycleInTides14(WL,Lat,'San_Francisco__CA',12,ALL_NAME,HC,0.25);
TIDE = HC.San_Francisco__CA.Tide_weekly;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% design for non-tidal residual
NTR = WL;
T_WL = WL(:,1) + WL(:,2)/12-1/24;
T_MSL = MSL_M(:,1)+MSL_M(:,2)/12-1/24;
BOOL0 = find(T_WL==T_MSL(1));
WL = WL(BOOL0(1):end,:);
TIDE = TIDE(BOOL0(1):end,:);
NON = RecSig_Non - MSL0;
ASA = [TotalPreMSL.Monthly(:,1),sum(TotalPreMSL.Monthly(:,3:4),2),...
       RecSig_4,RecSig_2,RecSig_3,RecSig_4,zeros(length(RecSig_1),1),RecSig_4];
for k = 1:length(ASA)       
    BOOL1 = find((WL(:,1)+WL(:,2)/12-1/24)==ASA(k));
    NTR(BOOL1,5) = WL(BOOL1,5) - sum(ASA(k,2:6)) - NON(k) - MSL0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate contribution

BOOL = find(WL(:,5)>NFL);
Contri = WL(BOOL,1:5);
Contri(:,6) = NFL;
Contri(:,7) = MSL0;

for k = 1:length(BOOL)

    T = (Contri(k,1) + Contri(k,2)/12-1/24);
    a = find(ASA(:,1)==T);

    if isempty(a)
        Contri(k,9:24) = NaN*ones(1,16);
    else
        Contri(k,8) = NON(a);
        Contri(k,9:15) = ASA(a,2:end);
        Contri(k,16) = sum(Contri(k,8:13));
        Contri(k,17) = Contri(k,15) - Contri(k,7);
        Contri(k,18) = MHHW;
        Contri(k,19) = TIDE(BOOL(k),1) + MSL0; %  MSL0 + tides

        a = Contri(k,19) - Contri(k,18); %above MHHW   
        c = sum(Contri(k,10:13)); % inter decadal
        ee = Contri(k,5) - Contri(k,19)  - sum(Contri(k,8:13)); % residual
        d = Contri(k,9); % seasonal
        b = Contri(k,8); %  nonlinear
        ABS = abs(a) +abs(c)+abs(d) + abs(ee)+ abs(b);

        Contri(k,21) = 100*b/ABS; %non linear 
        Contri(k,22) = 100*a/ABS; %above MHHW   
        Contri(k,23) = 100*d/ABS; % seasonal
        Contri(k,24) = 100*c/ABS; % inter decadal
        Contri(k,25) = 100*ee/ABS; % residual

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot water level decomposition

load time_plot
T1 = time_plot(:,1) + time_plot(:,2)/100 +time_plot(:,3)/10000 +time_plot(:,4)/1000000;
T2 = WL(:,1) + WL(:,2)/100 +WL(:,3)/10000 +WL(:,4)/1000000;
[a1,b1,c1] = intersect(T1,T2);       


figure(1);set(gcf,'unit','normalized ','position',[0.01,0.01,0.7,0.8]);
subplot(5,2,5)
plot(time_plot(b1,5),WL(c1,5),'linewidth',2);
title('Water level');ylabel('meter')
xlim([1950 2020]);grid on;set(gca,'fontsize',12)

subplot(5,2,2)
%         plot(ASA(:,1),NON+NOAAData.(NAME{j}).MSL,'linewidth',2);     
plot(ASA(:,1),NON,'linewidth',2);     

title('SLR')
xlim([1950 2020]);grid on;set(gca,'fontsize',12)

subplot(5,2,4)
plot(ASA(:,1),sum(ASA(:,3:6),2),'linewidth',2);         
title('ID')
xlim([1950 2020]);grid on;set(gca,'fontsize',12)

subplot(5,2,6)
plot(ASA(:,1),ASA(:,2),'linewidth',2);         
title('SC')
xlim([1950 2020]);grid on;set(gca,'fontsize',12)

subplot(5,2,8)
plot(time_plot(b1,5),TIDE(c1)+MSL0 - MHHW,'linewidth',2);
title('TA')
xlim([1950 2020]);grid on;set(gca,'fontsize',12)

subplot(5,2,10)

plot(time_plot(b1,5),NTR(c1,5)-TIDE(c1),'linewidth',2);
title('NTR')
xlim([1950 2020])
grid on;set(gca,'fontsize',12)



