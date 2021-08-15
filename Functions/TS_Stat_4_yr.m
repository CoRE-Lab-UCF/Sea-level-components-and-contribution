function [ts_yr,BestRC]=TS_Stat_4_yr(Time_Monatswerte,window); 
% _________________________________________________________________________
%
% PROGRAM "TS_Stat_4"
%
% Statistics tool for trend analysis of water level time series
%
% Version 1.1
%
% © fwu, Thomas Wahl
%
% Research Institute for Water and Environment
% University of Siegen
%
% Date: 02/11/2009
%
% Created as part of the AMSeL project (WBL199D)
%
% _________________________________________________________________________
%
% INPUT DATA:
%
% Water level time series
%
% Matrix "Time_Monatswerte", result of TideAnalyseMonat
%
%
% ANALYZES:
%
% (1) In the case of monthly values ??as input data, these
% initially determined ANNUAL VALUES, with the monthly values ??corresponding to
% Weighted number of days from which they were originally generated
% will.
%
% References: http://www.pol.ac.uk/psmsl/manuals/ioc_14ii.pdf (p.32)
%
%
% (2) The time series is EXTRAPOLATED beyond the ends to be at the
% later trend analyzes also indicate reliable values ??for the margins
% can. Various methods are used here
% - extrapolation via linear trend (trend of the first / last
% Window, serves here as an aid)
% - Through MC 100 AR1 realizations of those with the right to trend
% Generate initial time series, which are longer
% and their ends are appended to the output time series.
% This is based on the assumption that the "long-term
% Trend "of the extended time series is similar to the
% of the original time series

% OPTIONAL:
% - extrapolation over long-term mean (mean value of
% first / last time window)
% - extrapolation via linear trend (entire time series)
% - Last / first part (1 window length) of the output time series
% hang at the end mirrored horizontally
% - Last / first part (1 window length) of the output time series
% hang at the end mirrored horizontally and vertically


% References:
% MANN M.E. (2004) - On smoothing potentially non-stationary time series
% MOORE C.J. et al. (2005) New Tools for Analyzing Time Series
% Relationships and Trends
% BOX G. & JENKINS G. (1976) - Time series analysis - forecasting and
% control

% (3) STATIONARITY TESTS (KS test, Mann-Kendall test, sliding window test)
% Perform% for all EOFs to get those with trend information
% Extract%. Definition: EOF is considered unsteady if it is from
% at least two of the three tests are defined as such. The
% Stationarity tests are carried out for the extrapolated time series,
% at which the edge window is extrapolated via the linear trend. The
% EOFs that are recognized as unsteady are used for all others
% Reconstructions used.
%
% Credentials:
% VAN GELDER P. et al. (2008) - Data management of extreme marine and
% coastal hydro-meteorological events
% CHEN H.L. at al. (2002) - Testing Hydrologic Time Series for Stationarity
% MAN H.B. (1945) - Nonparametric Test Against Trend
%
%
% (4) Perform SINGULAR SYSTEM ANALYZES (SSA) for all time series
%
% References:
% GHIL M. et al. (2001) - Advanced spectral methods for climatic time series
% GOLYANDIDA N. et al. (2001) - Analisys of Time Series Structure - SSA and
% Related Techniques
%
% Time series of the NONLINEAR TRENDS with the relevant EOFs
% Reconstruct%, mean square error (MSE) in each case
% Determine% fit, the reconstruction with the smallest MSE
% Select% and determine 95% confidence intervals.
%
% Credentials:
% MANN M.E. (2004) - On smoothing potentially non-stationary time series
% JEVREJEVA S. et al. (2006) - Nonlinear trends and multiyear cycles in sea
% level records
%
%
% (5) MOVING LINEAR TRENDS (GLT) AND MOVING AVERAGE (GLM) over a
% Determine% or more time windows to be defined at the beginning. The
% Extrapolation of the ends is done using the method used in the case of the
% non-linear trends produced the best result. 95%
% Determine% confidence intervals for each linear trend.
%
% Credentials:
% HOLGATE S. (2007) - On the decadal rates of sea level change during the
% twentieth century
% HOLGATE S. et al. (2004) - Evidence for enhanced coastal sea level rise
% durcing the 1990s
%
%
% (6) LINEAR, EXPONENTIAL AND SQUARE FUNCTION to the whole
% Adjust% time series and determine 95% confidence intervals.
%
% Credentials:
% Jevrejeva et al. (2008) - Recent globa
% _________________________________________________________________________

%% ABFRAGEN:

% tic
% 
% gaugename=input('Bitte geben Sie den Pegelnamen ein:  ','s');
% disp(' ');
% % resolution=input('Handelt es sich um Monatswerte (im Moment gehen nur Monatswerte)? (ja: 0 / nein: 1)  ');
% % disp(' ');
% % window=input('Geben Sie eine Zeitfensterlänge X an (Empehlung:  X=N/5, N = Anzahl der Jahre):  ');
% % disp(' ');
% MC=input('Wieviele Monte-Carlo-Simulationen sollen durchgeführt werden?  ');

% EIGENSCHAFTEN DER EINGANGSZEITREIHE


% INPUT: nput should be annual mean sea level in this case as a matrix 
% with date in column 1 and water level in column 2


%% ANALYSESCHRITT (1)
% From the monthly values (from tide analysis month.m)
% annual values (weighted) are determined. An annual value is only then
% determined if at least 11 monthly values are available.
ts_yr=Time_Monatswerte;
MC = 10000;
gaugename = 'San Francisco';
% timeseries2=Time_Monatswerte2;

% ts=[datevec(timeseries(:,1)),timeseries(:,[2,3])];
% b=unique(ts(:,1));
% Length=length(b);
% ts_yr=zeros(Length,2);
% 
%     for i=1:Length;
%         y=find(ts==b(i,1));
%         if length(y)>=11
%         days=sum(ts([y(1,1):y(end)],8));
%         value_yr=sum(ts([y(1,1):y(end)],8).*ts([y(1,1):y(end)],7))./days;
%         a(:,:)=[b(i,1),value_yr];
%         ts_yr(i,:)=a;
%         end
%     end
%     
% gap=find(ts_yr(:,1)==0);
% ts_yr(gap,:)=[];
% gap=diff(ts_yr(:,1));
% gap=find(gap>1);
% gap=[1:max(gap)];
% ts_yr(gap,:)=[];
% 
% gap_mon=find(timeseries2(:,2)<min(ts_yr(:,1)));
% timeseries2(gap_mon,:)=[];

Length=length(ts_yr);
% window=round(Length/5);

% window=30;

% end

clear b y days value_yr a

% ANNUAL VALUES: here the time column has to be converted into serial date !!!!!
% if resolution==1
% ts_yr=timeseries;
% Length=length(ts_yr(:,1));
% end

%% ANALYSESCHRITT (2): Padding

% % 
% Extrapolation by a time window length with the mean of the last
% or first time window
% 
% Mean1_value=mean(ts_yr([1:window],2));
% Mean1_time=[ts_yr(1,1)-window:ts_yr(1,1)-1]';
% Mean1=Mean1_time;
% Mean1(:,2)=Mean1_value;
% Mean2_value=mean(ts_yr([Length-window+1:Length],2));
% Mean2_time=[ts_yr(Length,1)+1:ts_yr(Length,1)+floor(window)]';
% Mean2=Mean2_time;
% Mean2(:,2)=Mean2_value;
% ts_exp1=[Mean1;ts_yr;Mean2];
% 
% clear Mean1_value Mean1_time Mean2_value Mean2_time 

% Extrapolation by a time window length with a linear trend of the
% first / last window

Mean1_value=mean(ts_yr([1:window],2));
Mean1_time=[ts_yr(1,1)-window:ts_yr(1,1)-1]';
Mean1=Mean1_time;
Mean1(:,2)=Mean1_value;

Mean2_value=mean(ts_yr([Length-window+1:Length],2));
Mean2_time=[ts_yr(Length,1)+1:ts_yr(Length,1)+floor(window)]';
Mean2=Mean2_time;
Mean2(:,2)=Mean2_value;

trend_para=polyfit(ts_yr([1:window],1),ts_yr([1:window],2),1);

for i=1:window
    Mean1(i,2)=Mean1(i,1)*trend_para(1,1)+trend_para(1,2);
end

trend_para=polyfit...
(ts_yr([Length-window+1:Length],1),ts_yr([Length-window+1:Length],2),1);

for i=1:window
    Mean2(i,2)=Mean2(i,1)*trend_para(1,1)+trend_para(1,2);
end
ts_exp2=[Mean1;ts_yr;Mean2];

% LL=length(ts_exp2);

% Determine the linear trend of the time series
trend_para=polyfit(ts_yr([1:Length],1),ts_yr([1:Length],2),1);

% 
% Mean1(:,2)=Mean1(:,1)*trend_para(1,1)+trend_para(1,2);
% Mean2(:,2)=Mean2(:,1)*trend_para(1,1)+trend_para(1,2);
% 
% ts_exp3=[Mean1;ts_yr;Mean2];
% 
% % Hang the last / first part mirrored horizontally at the end
% 
% Mean1(:,2)=flipud(ts_yr([2:window+1],2));
% Mean2(:,2)=flipud(ts_yr([Length-window:Length-1],2));
% 
% ts_exp4=[Mean1;ts_yr;Mean2];
% 
% % Hang the last / first part at the end, mirrored horizontally and vertically
% 
% Mean1(:,2)=-flipud(ts_yr([2:window+1],2))+2*ts_yr(1,2);
% Mean2(:,2)=-flipud(ts_yr([Length-window:Length-1],2))+2*ts_yr(Length,2);
% 
% ts_exp5=[Mean1;ts_yr;Mean2];
% 
% clear Mean1 Mean2

% Through Monte Carlo simulation 100 AR1 realizations of
% Generate starting time series and adopt ends

% Form differences in order to generate stationary time series

% Find the best parametric fit

% linear model

[cfun,gof]=fit(ts_yr(:,1),ts_yr(:,2),'poly1');
pConf = predint(cfun,ts_yr(:,1),0.95,'functional','off');
cfun_ex=cfun(ts_exp2(:,1));
cfun=cfun(ts_yr(:,1));



% square model

[cfun1,gof1]=fit(ts_yr(:,1),ts_yr(:,2),'poly2');
pConf1 = predint(cfun1,ts_yr(:,1),0.95,'functional','off');
cfun_ex1=cfun1(ts_exp2(:,1));
cfun1=cfun1(ts_yr(:,1));

% exponential (s) Modell (e)

[cfun2,gof2]=fit(ts_yr(:,1),ts_yr(:,2)+abs(min(ts_yr(:,2))),'exp1');
pConf2 = predint(cfun2,ts_yr(:,1),0.95,'functional','off');
cfun_ex2=cfun2(ts_exp2(:,1));
cfun_ex2=cfun_ex2-abs(min(ts_exp2(:,2)));
cfun2=cfun2(ts_yr(:,1));
cfun2=cfun2-abs(min(ts_yr(:,2)));
pConf2=pConf2-abs(min(ts_yr(:,2)));

Errors=[gof.rmse;gof1.rmse;gof2.rmse];

BestParaFit=find(Errors==min(Errors));


if BestParaFit==1;
    cfun_best=cfun;
    cfun_ex_best=cfun_ex;
    pConf_best=pConf;
    Str='Polynom 1. Ordnung';
end

if BestParaFit==2;
    cfun_best=cfun1;
    cfun_ex_best=cfun_ex1;
    pConf_best=pConf1;
    Str='Polynom 2. Ordnung';
end

if BestParaFit==3;
    cfun_best=cfun2;
    cfun_ex_best=cfun_ex2;
    pConf_best=pConf2;
    Str='Exponentialfunktion';
end

D=ts_yr(:,2)-cfun_best;

% Determine parameters for AR1 model

% [a,e]=arcov(D,1);
% e=sqrt(e);

[g,a,mu2]=ar1(D);

% Generate 10,000 realizations

% X=zeros(Length+2*window,MC);
% X(1,:)=sqrt(a^2/(1-g^2))*randn(1,MC);
% z=e*randn(Length+2*window,MC);

Noise=ar1noise(Length+2*window,MC,g,a);
X=Noise;

% for i=2:Length+2*window
%   X(i,:)=g*X(i-1,:)+Noise(i,:);
% end

clear a e z

% Integration / intersection with lin. Trend of the initial time series

for i=1:MC
    X([window+1:Length+window],i)=D;
    X(:,i)=X(:,i)+cfun_ex_best;
end

X(:,:)=X-(X(window+1)-ts_yr(1,2));

% Write all "synthetically" extended time series into a matrix

% X(:,[MC+1:MC+5])=...
%     [ts_exp1(:,2),ts_exp2(:,2),ts_exp3(:,2),ts_exp4(:,2),ts_exp5(:,2)];

%% ANALYSESCHRITT (3): Stationarity tests

% Determine which EOFs are to be accepted as decisive later

X1=ts_yr(:,2); 
    
% SSA-Schritt 1: Embeddding
N=length(X1); 
K=N-window+1; 
Y=zeros(window,K);  
for i=1:K;
    Y(1:window,i)=X1(i:window+i-1);   
end
    
% SSA-step 2: Singular Value Decomposition (SVD)
S=Y*Y'; 
[U,lamda,V]=svd(S);
[d,i]=sort(-diag(lamda));     
d=-d;
Sum=sum(d); 

V=Y'*U; 

% I1=zeros(MC+5,5);

I1=zeros(1,5);

% 2-D KOLMOGOROV-SMIRNOV TEST

% Creation of a matrix with part time series of the length window, always at
% dt = 1 moved (sliding window)

for b=1:5;
    nw1 = window;
    nts = length(V) - (window-1); % Number of possible part-time series
    ts = ones(window,nts);     % Definition of a matrix for part-time series
    D=V(:,b);
    
for i=1:nts;
    ts(:,i) = D(i:nw1);
    nw1 = nw1 + 1;
end

% The first part time series is set as the reference time series. The others
% Time series are always compared with this, whether this is a common one
% Have density function (two-sample K-S test)

% if H = 0: same pdf; if H = 1: different pdf

for i=1:nts;
    [H(i,1), P(i,1), ksstat(i,1)]= kstest2(ts(:,1),ts(:,i),0.05);
end

n1 = sum(H(:,1));  % Determining the number of '1'
proz(b,1) = n1/nts;     % Percentage of '1'

Sign=find(proz(:,1)>=0.3); % Value from sensitivity analyzes
proz(:,1)=0;
proz(Sign,1)=1;


% SLIDING-WINDOW-TEST

x = window;  
D=V(:,b);

% Calculation of the mean value within an x-year time window,
% which is always shifted by 1 year

j=1;
for i=x:length(V);
    Mean = mean(D(j:i,:));
    j = j+1;
    
    m19v(i,1) =Mean;
    
end

% The calculated value is always based on the last position of the time window
% so that the first (x-1) positions are filled with the xth value
% will

m19v(1:(x-1),1) = m19v(x,1);

% Calculation of the 95% confidence range for the mean
% of the first time window of length x

meanm19v = mean (D((1:x),:));
% varm19v  = var  (D((1:x),:));
stdm19v =  std (D((1:x),:));
a = ones(length(V),1);

cm19v = 1.960;  % Table value for 95% percentile from standard normal distribution
km19v = cm19v * stdm19v / sqrt(x);

m19ciupm   = (meanm19v + km19v) .* a;
m19cidownm = (meanm19v - km19v) .* a;

u1 = find(m19v > m19ciupm | m19v < m19cidownm);
u1 = length(u1);

proz(b,2) = u1/nts;

% Limit value function is the result of sensitivity analyzes

Sign=find(proz(:,2)>=(0.6881*window^-0.0713));
proz(:,2)=0;
proz(Sign,2)=1;

% MANN-KENDALL-TEST

D=V(:,b);
t=[1:length(D)]';
datain = [t,D];
[taub h sig Z S sigma sen] = ktaub(datain, 0.05);

proz(b,3) = h;

end

proz(:,4)=sum(proz');

I=find(proz(:,4)==3);
I=I';

Anz=length(I);
I1(1,[1:Anz])=[I];
gg=find(I1==0);
I1(gg)=[];
I=I1;

if length(I)==0;
    I=[1];
end
clear proz

%% Analyseschritt 4: SSA

RMSE_start=zeros(MC,1);
RMSE_end=zeros(MC,1);

Rec=zeros(Length,length(a));
Deviation=zeros(Length,length(a));

% I=[1];

for a=1:MC;

    X1=X(:,a);

    % SSA-Schritt 1: Embeddding
    N=length(X1); 
    K=N-window+1; 
    Y=zeros(window,K);  
    for i=1:K;
        Y(1:window,i)=X1(i:window+i-1);   
    end

    % SSA-Schritt 2: Singular Value Decomposition (SVD)
    S=Y*Y'; 
    [U,lamda,V]=svd(S);
    [d,i]=sort(-diag(lamda));     
    d=-d;
    Sum=sum(d); 

    V=Y'*U; 

    % SSA-Schritt 4: Rekonstruktion

    I=[1];

    Vt=V';
    rca=U(:,I)*Vt(I,:);

    y=zeros(N,1);  
    Lp=min(window,K);
    Kp=max(window,K);

    for k=0:Lp-2;
        for m=1:k+1;
            y(k+1)=y(k+1)+(1/(k+1))*rca(m,k-m+2);
        end
    end

    for k=Lp-1:Kp-1;
        for m=1:Lp;
            y(k+1)=y(k+1)+(1/(Lp))*rca(m,k-m+2);
        end
    end

    for k=Kp:N;
        for m=k-Kp+2:N-Kp+1;
            y(k+1)=y(k+1)+(1/(N-k))*rca(m,k-m+2);
        end
    end

    r=X(:,a)-y;
       %vr=(sum(d(I))/sev)*100;

    % Rec1(:,a)=y;

    L1=length(X1);   
    y([L1-window+1:L1],:)=[];    
    y([1:window],:)=[];   

    Rec(:,a)=y;
    Deviation(2:Length,a)=diff(y)*10;
    clear y

    L=length(ts_yr);

    s_end=sqrt(mean((ts_yr([L-window+1:L],2)-Rec([L-window+1:L],a)).^2));
    s_start=sqrt(mean((ts_yr([1:window],2)-Rec([1:window],a)).^2));
    RMSE_start(a,1)=s_start;
    RMSE_end(a,1)=s_end;

end
% end

% Since it happens that additional EOFs as
% are considered authoritative, there may be several reconstructions at the end
% evaluated

% Mode=mode(I1);
% 
% NM2=find(I1(:,2)~=Mode(:,2));
% NM3=find(I1(:,3)~=Mode(:,3));
% NM4=find(I1(:,4)~=Mode(:,4));
% NM5=find(I1(:,5)~=Mode(:,5));
% 
% NotMode=[NM2;NM3;NM3;NM4;NM5]
% NotMode=unique(NotMode);
% 
% RMSE2=RMSE(NotMode,:);
% RMSE(NotMode,:)=[];
%     
% I2=I1(NotMode,:);
% I1(NotMode,:)=[];
%    
% X2=X(:,NotMode);
% X(:,NotMode)=[];
%     
% Rec2(:,:)=Rec(:,NotMode);
% Rec(:,NotMode)=[];

% Best Fit

MinError_start=min(RMSE_start);
bfit_start=find(RMSE_start==MinError_start);

MinError_end=min(RMSE_end);
bfit_end=find(RMSE_end==MinError_end);

% Verwendete EOFs zur Rekonstruktion mit kleinstem RMSE

% EOFs=find(I1(bfit,:)>0);

% bestExp_start=X(1:L1-window,bfit_start);
% bestExp_end=X(L1-window+1:L,bfit_start);
% bestExp=[bestExp_start;bestExp_end];

Rec_start=Rec(1:L-window,bfit_start);
Rec_end=Rec(L-window+1:L,bfit_end);

BestRC=[Rec_start;Rec_end];

BestDeviate=diff(BestRC)*10;

% EOFs=I1(bfit,[1:length(EOFs)]);

Error_SSA_Ges=sqrt(mean((ts_yr(:,2)-BestRC(:,1)).^2));

Error_SSA_Teil=sqrt(mean((ts_yr(11:L-10,2)-BestRC(11:L-10,1)).^2));

% possibly second reconstruction, depending on the type of extrapolation
% Output time series additional EOFs can be decisive

% if I2(:,1)~=0;
% MinError2=min(RMSE2);
% bfit2=find(RMSE2==MinError2);
% 
% EOFs2=find(I2(bfit2,:)>0);
% bestExp2=X2(:,bfit2);
% BestRC2=Rec2(:,bfit2);
% EOFs2=I2(bfit2,[1:length(EOFs2)]);
% MinErrorGes1=sqrt(mean((ts_yr(:,2)-BestRC2(:,1)).^2));
% end

% % 95% Konfidenzintervalle für Rekonstruktion(en)
%  
%  for i=1:L
%     [mue,sigma]=normfit(detrend(bestExp([i:i+2*window-1],1)));
%     confup=BestRC(i,1)+((sigma/sqrt(2*window))*1.96); % Tabellenwert (s. Papula III S. 517) 
%     confdown=BestRC(i,1)-((sigma/sqrt(2*window))*1.96);
%     conf(i,1)=confup;
%     conf(i,2)=confdown;
%     confup1=BestRC(i,1)+((sigma/sqrt(2*window))); % Tabellenwert (s. Papula III S. 517) 
%     confdown1=BestRC(i,1)-((sigma/sqrt(2*window)));
%     conf1(i,1)=confup1;
%     conf1(i,2)=confdown1;
%  end
% 
% if I2(:,1)~=0;
% 
% 
%  for i=1:L
%     [mue,sigma]=normfit(detrend(bestExp2([i:i+2*window-1],1)));
%     confup2=BestRC2(i,1)+((sigma/sqrt(2*window))*1.96); % Tabellenwert (s. Papula III S. 517) 
%     confdown2=BestRC2(i,1)-((sigma/sqrt(2*window))*1.96);
%     conf2(i,1)=confup2;
%     conf2(i,2)=confdown2;
%     confup3=BestRC2(i,1)+((sigma/sqrt(2*window))); % Tabellenwert (s. Papula III S. 517) 
%     confdown3=BestRC2(i,1)-((sigma/sqrt(2*window)));
%     conf3(i,1)=confup3;
%     conf3(i,2)=confdown3;
%  end  
% end

%% ANALYSESCHRITT (5): GLT und GLM

% sliding linear trends (20 year olds !!!)

window1=20;

% % mit Extrapolation
% 
% glintrend=zeros(L,2);
% 
% for i=round(window1/2)+1:round(L+window1/2)
% cfun = fit(ts_exp2([i:i+window1-1],1),bestExp([i:i+window1-1]),'poly1');    
% confBounds=confint(cfun);
% trend_year_MSL=cfun.p1*10;
% glintrend(i-round(window1/2),1)=confBounds(2,1)*10-trend_year_MSL;
% glintrend(i-round(window1/2),2)=trend_year_MSL;
% end

% ohne Extrapolation

glintrend=zeros(L-window1,2);

for i=1:L-window1
cfun = fit(ts_yr(i:i+window1-1,1),ts_yr(i:i+window1-1,2),'poly1');    
confBounds=confint(cfun);
trend_year_MSL=cfun.p1*10;
glintrend(i,1)=confBounds(2,1)*10-trend_year_MSL;
glintrend(i,2)=trend_year_MSL;
end
    
    
% gleitende Mittel

% % mit Extrapolation
% 
% glmean=zeros(L,1);
% 
% for i=round(window1/2)+1:round(L+window1/2)
% Mean = mean(bestExp([i:i+window1-1],1));    
% glmean(i-round(window1/2),1)=Mean;
% end    

% ohne Extrapolation

glmean=zeros(L-window1,1);

for i=1:L-window1
Mean = mean(ts_yr([i:i+window1-1],2));    
glmean(i,1)=Mean;
end 

Error_GLM=sqrt(mean((ts_yr(window1/2+1:L-window1/2,2)-glmean(:,1)).^2));
 

%% OUTPUT

% GRAPHIC REPRESENTATION (S)
%__________________________________________________________________________
%
% SSA & 1. Derivation
%__________________________________________________________________________
%
s=zeros(window,2);

for i=1:window
    Min_start=min(Rec(i,:));
    Max_start=max(Rec(i,:));
    s(i,1)=Min_start;
    s(i,2)=Max_start;
end

s_Dev=zeros(window-1,2);

for i=1:window-1
    Min_start=min(Deviation(i+1,:));
    Max_start=max(Deviation(i+1,:));
    s_Dev(i,1)=Min_start;
    s_Dev(i,2)=Max_start;
end

s_Dev(window,:)=BestDeviate(window);

% s_Dev=diff(s)*10;

e=zeros(L,2);

for i=L-window+1:L
    Min_end=min(Rec(i,:));
    Max_end=max(Rec(i,:));
    e(i,1)=Min_end;
    e(i,2)=Max_end;
end

e([1:L-window],:)=[];

e_Dev=zeros(L,2);

for i=L-window+2:L
    Min_end=min(Deviation(i,:));
    Max_end=max(Deviation(i,:));
    e_Dev(i,1)=Min_end;
    e_Dev(i,2)=Max_end;
end

e_Dev([1:L-window+1],:)=[];

e_1(1,1)=BestDeviate(L-window);
e_1(1,2)=BestDeviate(L-window);
e_Dev=[e_1;e_Dev];

% e_Dev=diff(e)*10;


% figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
% 
% subplot(2,1,1,'Parent',figure1,'XTickLabel','','FontSize',14);
% % Create axes
% % axes1 = axes('Parent',figure1,'XTickLabel','',...
% %     'Position',[0.13 0.3876 0.775 0.5374],...
% %     'FontSize',14);
% % Uncomment the following line to preserve the X-limits of the axes
% xlim([1840 2010]);
% % Uncomment the following line to preserve the Y-limits of the axes
% miny=min(ts_yr(:,2))-5;
% maxy=max(ts_yr(:,2))+5;
% 
% ylim([miny maxy]);
% 
% box('on');
% hold('all');
% 
% L3=fill([ts_yr([1:window],1)' fliplr(ts_yr([1:window],1)')],[s(:,1)' fliplr(s(:,2)')],[0.502 0.502 0.502]);
% hold all
% fill([ts_yr([L-window+1:L],1)' fliplr(ts_yr([L-window+1:L],1)')],[e(:,2)' fliplr(e(:,1)')],[0.502 0.502 0.502])
% hold all
% plot(ts_yr([1:window],1),s(:,1),'LineWidth',3,'Color',[0.502 0.502 0.502])
% plot(ts_yr([1:window],1),s(:,2),'LineWidth',3,'Color',[0.502 0.502 0.502])
% plot(ts_yr([L-window+1:L],1),e(:,1),'LineWidth',3,'Color',[0.502 0.502 0.502])
% plot(ts_yr([L-window+1:L],1),e(:,2),'LineWidth',3,'Color',[0.502 0.502 0.502])
% 
% L1=plot(ts_yr(:,1),ts_yr(:,2),'LineWidth',3,'Color',[0 0 0]);
% hold all
% L2=plot(ts_yr(:,1),BestRC,'LineWidth',3,'Color',[0 0 1]);
% 
% 
% % Create ylabel
% ylabel({'Water level [cmPN]'},'FontSize',16);
% 
% legend1 = legend([L1,L2,L3],'MSL-Annual valueS','SSA-Reconstruction with the smallest RMSE',...
%     'Uncertainty  (from 10,000 MCAP simulations)');
% set(legend1,'Position',[0.146 0.8201 0.4031 0.08712]);
% 
% subplot(2,1,2,'Parent',figure1,'FontSize',14);
% % Create axes
% % axes2 = axes('Parent',figure1,'Position',[0.13 0.11 0.775 0.261],...
% %     'FontSize',14);
% % Uncomment the following line to preserve the X-limits of the axes
% xlim([1840 2010]);
% % Uncomment the following line to preserve the Y-limits of the axes
% ylim([-4.99 11.99]);
% box('on');
% hold('all');
% 
% % Create plot
% L5=fill([ts_yr([2:window+1],1)' fliplr(ts_yr([2:window+1],1)')],[s_Dev(:,1)' fliplr(s_Dev(:,2)')],[0.8 0.8 0.8]);
% hold all
% fill([ts_yr([L-window+1:L],1)' fliplr(ts_yr([L-window+1:L],1)')],[e_Dev(:,2)' fliplr(e_Dev(:,1)')],[0.8 0.8 0.8])
% hold all
% plot(ts_yr([2:window+1],1),s_Dev(:,1),'LineWidth',3,'Color',[0.8 0.8 0.8])
% plot(ts_yr([2:window+1],1),s_Dev(:,2),'LineWidth',3,'Color',[0.8 0.8 0.8])
% plot(ts_yr([L-window+1:L],1),e_Dev(:,1),'LineWidth',3,'Color',[0.8 0.8 0.8])
% plot(ts_yr([L-window+1:L],1),e_Dev(:,2),'LineWidth',3,'Color',[0.8 0.8 0.8])
% 
% L4=plot(ts_yr(2:L,1),BestDeviate,'LineWidth',3,'Color',[0 0 1]);
% plot([1840;2010],[0;0],'LineWidth',1,'Color',[0 0 0]);
% 
% % Create xlabel
% xlabel({'Jahr'},'FontSize',16);
% 
% % Create ylabel
% ylabel({'Anstiegsrate MSL aus','SSA-Rekonstruktion [mm/a]'},'FontSize',16);
% 
% legend2 = legend([L4,L5],'Anstiegsrate MSL (1. Ableitung der SSA-Rekonstruktion mit dem kleinsten RMSE)','Unsicherheitsbereich (aus 10000 MCAP-Simulationen)');
% set(legend2,'Position',[0.1459 0.3698 0.582 0.0598]);
% 
% saveas(gcf,sprintf('SSA-MCAP %s ', gaugename));
% 
% %__________________________________________________________________________
% %
% % GLT & GLM
% %__________________________________________________________________________
% 
% figure2 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
% 
% % GLT
% 
% % Create axes
% axes2 = axes('Parent',figure2,'XTickLabel','',...
%     'Position',[0.13 0.5515 0.775 0.4418],...
%     'FontSize',14);
% % Uncomment the following line to preserve the X-limits of the axes
% xlim([1840 2010]);
% % Uncomment the following line to preserve the Y-limits of the axes
% ylim([-11.9 11.9]);
% box('on');
% hold('all');
% 
% % Create plot
% L1=plot(ts_yr(window1/2+1:L-window1/2,1),glintrend(:,2),'Parent',axes2,'LineWidth',3,...
%     'DisplayName','20-jähriger gleitender linearer Trend');
% 
% % Create errorbar
% L2=errorbar(ts_yr(window1/2+1:L-window1/2,1),glintrend(:,2),glintrend(:,1),'Marker','+','LineStyle','none',...
%     'DisplayName','95%-Konfidenzintervalle',...
%     'Color',[0.3137 0.3137 0.3137],...
%     'Parent',axes2);
% 
% plot([1840;2010],[0;0],'Parent',axes2,'LineWidth',1,'Color',[0 0 0]);
% 
% % Create ylabel
% ylabel({'Anstiegsrate MSL aus','Jahreswerten [mm/a]'},'FontSize',16);
% 
% % Create legend
% legend1 = legend([L1,L2],'20-jähriger gleitender linearer Trend','95%-Konfidenzintervalle');
% set(legend1,'Position',[0.1471 0.9141 0.2805 0.0598]);
% 
% % GLM
% 
% % Create axes
% axes3 = axes('Parent',figure2,'Position',[0.13 0.07309 0.775 0.4673],...
%     'FontSize',14);
% % Uncomment the following line to preserve the X-limits of the axes
% xlim([1840 2010]);
% % Uncomment the following line to preserve the Y-limits of the axes
% ylim([miny maxy]);
% box('on');
% hold('all');
% 
% % Create plot
% plot(ts_yr(:,1),ts_yr(:,2),'Parent',axes3,'LineWidth',3,'DisplayName','MSL-Jahreswerte',...
%     'Color',[0 0 0]);
% 
% % Create plot
% plot(ts_yr(window1/2+1:L-window1/2,1),glmean,'Parent',axes3,'LineWidth',3,'Color',[0 0 1],...
%     'DisplayName','20-jähriges gleitendes Mittel');
% 
% % Create xlabel
% xlabel({'Jahr'},...
%     'FontSize',16);
% 
% % Create ylabel
% ylabel({'Wasserstand [cmPN]'},'FontSize',16);
% 
% % Create legend
% legend2 = legend(axes3,'show');
% set(legend2,'Position',[0.1488 0.4551 0.2297 0.0598]);
% 
% saveas(gcf,sprintf('GLT-GLM %s ', gaugename));
% 
% %__________________________________________________________________________
% %
% % ParaFit
% %__________________________________________________________________________
% 
% % t_mon=timeseries2(:,2)+timeseries2(:,1)*1/12;
% 
% figure3 = figure('PaperType','a4letter','PaperSize',[20.98 29.68]);
% 
% % Create axes
% axes1 = axes('Parent',figure3,'FontSize',14);
% % % Uncomment the following line to preserve the X-limits of the axes
% xlim([1840 2010]);
% % % Uncomment the following line to preserve the Y-limits of the axes
% 
% % miny=min(timeseries2(:,4))-5;
% % maxy=max(timeseries2(:,4))+20;
% % 
% % ylim([miny maxy]);
% 
% box('on');
% hold('all');
% 
% % % Create plot
% % L1=plot(t_mon(:,1),timeseries2(:,4),'Parent',axes1,'LineWidth',2,'Color',[0.5 0.5 0.5],...
% %     'DisplayName','MSL-Monatswerte');
% 
% L4=fill([ts_yr(:,1)' fliplr(ts_yr(:,1)')],[(pConf_best(:,1))' fliplr((pConf_best(:,2))')],[0.8 0.8 0.8]);
% plot(ts_yr(:,1),pConf_best(:,1),'Parent',axes1,'LineWidth',2,'Color',[0.8 0.8 0.8])
% plot(ts_yr(:,1),pConf_best(:,2),'Parent',axes1,'LineWidth',2,'Color',[0.8 0.8 0.8])
% 
% % Create multiple lines using matrix input to plot
% L2=plot(ts_yr(:,1),ts_yr(:,2),'Parent',axes1,'LineWidth',3,'Color',[0 0 0]);
% 
% L3=plot(ts_yr(:,1),cfun_best,'Parent',axes1,'LineWidth',3,'Color',[0 0 1]);
% 
% % Create xlabel
% xlabel({'Jahr'},'FontSize',16);
% 
% % Create ylabel
% ylabel({'Wasserstand [cmPN]'},'FontSize',16);
% 
% % Create legend
% legend1 = legend([L2,L3,L4],'MSL-Jahreswerte',Str,'95%-Konfidenzintervall');
% set(legend1,'Position',[0.1563 0.7863 0.1961 0.1144]);
% 
% Trend=num2str(trend_para(1)*10,'%1.2f');
% Start=num2str(ts_yr(1,1));
% End=num2str(ts_yr(length(ts_yr),1));
% 
% % Create textbox
% annotation(figure3,'textbox',...
%     'String',{['Lin. Trend (',Start,' - ',End,') = ',Trend, ' mm/a']},...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'Position',[0.1565 0.7331 0.2607 0.03876]);
% 
% saveas(gcf,sprintf('Para-Fit %s ', gaugename));
% 
% Errors_Lin_Quad_Exp=Errors;
% 
% savefile1= sprintf ('Ergebnis_Auswertung %s ', gaugename);
% save (savefile1,'ts_yr','Errors_Lin_Quad_Exp','Error_SSA_Ges','Error_SSA_Teil','Error_GLM','e','s','BestRC','glmean','glintrend','I','Rec');
% 
% % toc
% 
% %end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
