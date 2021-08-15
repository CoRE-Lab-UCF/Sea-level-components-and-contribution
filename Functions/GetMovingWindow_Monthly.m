function [TotalPreMSL] = GetMovingWindow_Monthly(MSL,TIME,Window)

% PROGRAM "GetMovingWindow_Monthly"
% Tool to estimate monthly amplitude and phase of seasonal cycle in MSL with moving least squares. 
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. MSL: monthly MSL. A vector.   
%       3. TIME: time. A vector saves [year+month/12-1/24].
%       2. Window: moving window length for least squares. unit: year
%
% Output: 
%       1. TotalPreMSL: a struct saves two matrixs: 'Monthly' and 
%          'AmpPha'. 
%          For 'Monthly', it has 10 columns: 
%             (1) time
%             (2) reconstructed trend
%             (3) reconstruct semi-annual cycle
%             (4) reconstruct annual cycle
%             (5) semi-annual amplitude
%             (6) annual amplitude
%             (7) std of semi-annual amplitude
%             (8) std of annual amplitude
%             (9) std of semi-annual cycle
%             (10) std of annual cycle
%
%          For 'AmpPha', it has 5 columns: 
%             (1) time
%             (2) semi-annual amplitude
%             (3) annual amplitude
%             (4) semi-annual phase
%             (5) annual phase
%
% Reference: Wahl et al., 2014. Rapid changes in the seasonal sea level cycle 
%            along the US Gulf coast from the late 20th century. 


Window = 5; % window: scalar (years) , like:  5 

m = length(MSL);
CONT = 0;
PreMSL = [];
PreSASSA = [];
for j = 1:m-Window*12+1
    
    L = MSL(j:j+Window*12-1);   % build matrix
%%  Sida's Method - Time datum & Phase are always changing      
%     t = [1:Window*12]';
%     A(1:Window*12,1) = ones(Window*12,1);
%     A(1:Window*12,2) = t/Window*12
%     A(1:Window*12,3) = sin(2*pi/6*t);    % semi-annual cycle
%     A(1:Window*12,4) = cos(2*pi/6*t);    %             amplitude
%     A(1:Window*12,5) = sin(2*pi/12*t);   % annual      cycle
%     A(1:Window*12,6) = cos(2*pi/12*t);   %             amplitude
   
%%  Time datum based on real time       
    t = [TIME(j):1/12:TIME(j+Window*12-1)]';
    A(1:Window*12,1) = ones(Window*12,1);
    A(1:Window*12,2) = t/Window*12;
    A(1:Window*12,3) = sin(2*pi/0.5*t);
    A(1:Window*12,4) = cos(2*pi/0.5*t);
    A(1:Window*12,5) = sin(2*pi/1*t);
    A(1:Window*12,6) = cos(2*pi/1*t);
    
    AA = A;
    BOOL = isnan(L);
    if sum(BOOL)>=12*2  % set threshold stop fitting (if missing data is too much) 
        
        if j == 1  % the beginning of the fitting
            Beginning(1:Window*6-1,1) = TIME(1:Window*6-1,1);
            Beginning(1:Window*6-1,2:4) = NaN*ones(Window*6-1,3);
            X = [1:6]';PreErr = [1 1];
            PreMSL(j,1) = TIME(Window*6+j-1);
            PreMSL(j,2) = NaN; 
        elseif j==m-Window*12+1   % the end of the fitting 
            Ending(1:Window*6+1,1) = TIME(end-Window*6:end,1);
            Ending(1:Window*6+1,2:4) = NaN*ones(Window*6+1,3);
            X = [1:6]';PreErr = [1 1];           
        else 
            X = [1:6]';PreErr = [1 1];       
            PreMSL(j,1) = TIME(Window*6+j-1);
            PreMSL(j,2) = NaN; 
        end
        AmpPha(j,1) = TIME(Window*6+j-1,1);
        AmpPha(j,2:5) = NaN;                
    else  
        A(BOOL,:) = [];
        L(BOOL) = [];        
        X = inv(A'*A)*A'*L;  % least squares analysis
        V = A*X-L;
        sigma = sqrt(V'*V/(Window*12-sum(BOOL)-6));  % std
        DX = sqrt(diag(sigma^2*inv(A'*A))); % std of each parameter
        t = abs(X)./DX; % significance test (one sigma level)
        TEMP = AA(:,1:2)*X(1:2);  % reconstruct trend
        PreMSL(j,1) = TIME(Window*6+j-1,1);% + TIME(Window*6+j-1,2)/12-1/12;
        PreMSL(j,2) = TEMP(Window*6);        
        TEMP = AA(:,3:4)*X(3:4); % reconstruct semi-annual cycle
        PreMSL(j,3) = TEMP(Window*6);               
        TEMP = AA(:,5:6)*X(5:6); % reconstruct annual cycle
        PreMSL(j,4) = TEMP(Window*6);
        PreMSL(j,5) = norm(X(3:4));   %  semi-annual amplitude
        PreMSL(j,6) = norm(X(5:6));   %  annual amplitude
        
%         PreMSL(j,11) = (1/(1+X(3)^2/X(4)^2)*1/X(4))^2*DX(3)^2 + (1/(1+X(3)^2/X(4)^2)*X(3)/X(4)^2)^2*DX(4)^2;
%         PreMSL(j,12) = (1/(1+X(5)^2/X(6)^2)*1/X(6))^2*DX(5)^2 + (1/(1+X(5)^2/X(6)^2)*X(5)/X(6)^2)^2*DX(6)^2;
        
%%%%%%%%%%%%%%%%%%%%% replaced j,7-10 below %%%%%%%%%%%%%%%%%%%%%%
%       PreMSL(j,7) = sqrt((X(3)^2*DX(3)^2 +  X(4)^2*DX(4)^2)/(X(3)^2 + X(4)^2));  % std of semi-annual phase
%       PreMSL(j,8) = sqrt((X(5)^2*DX(5)^2 +  X(6)^2*DX(6)^2)/(X(5)^2 + X(6)^2));   % std of annual phase
                                   
%       PreMSL(j,9) = sqrt(AA(Window*6,3)^2*DX(3)^2  + AA(Window*6,4)^2*DX(4)^2);
%       PreMSL(j,10) = sqrt(AA(Window*6,5)^2*DX(5)^2  + AA(Window*6,6)^2*DX(6)^2);
%%%%%%%%%%%%%%%%%%%%%% REPLACED lines of PreMSL(j,7-10) with: %%%%%%%%%
        DXX = sigma^2*inv(A'*A);
        TEMP1 = [X(3)/sqrt(X(3)^2 + X(4)^2), X(4)/sqrt(X(3)^2 + X(4)^2)];
 %%%    TEMP2 = [X(5)/sqrt(X(3)^2 + X(4)^2), X(6)/sqrt(X(5)^2 + X(6)^2)];
        TEMP2 = [X(5)/sqrt(X(5)^2 + X(6)^2), X(6)/sqrt(X(5)^2 + X(6)^2)];                                   
        %PreMSL(j,7) = TEMP1*DXX(3:4,3:4)*TEMP1';  % std of semi-annual amplitude
        %PreMSL(j,8) = TEMP2*DXX(5:6,5:6)*TEMP2';  % std of annual amplitude
        %PreMSL(j,9) = AA(Window*6,3:4)*DXX(3:4,3:4)*AA(Window*6,3:4)'; % std of semi-annual cycle
        %PreMSL(j,10) = AA(Window*6,5:6)*DXX(5:6,5:6)*AA(Window*6,5:6)';  % std of annual cycle"
        %%Replaced Lines 97-100 (j7-10) with: %%%%%%%%%%%%%%%%%          
        PreMSL(j,7) = sqrt(TEMP1*DXX(3:4,3:4)*TEMP1');  % std of semi-annual amplitude
        PreMSL(j,8) = sqrt(TEMP2*DXX(5:6,5:6)*TEMP2');  % std of annual amplitude
        PreMSL(j,9) = sqrt(AA(Window*6,3:4)*DXX(3:4,3:4)*AA(Window*6,3:4)'); % std of semi-annual cycle
        PreMSL(j,10) = sqrt(AA(Window*6,5:6)*DXX(5:6,5:6)*AA(Window*6,5:6)');  % std of annual cycle"
        PreErr = [PreMSL(j,9) PreMSL(j,10)];
        if j==1
            Beginning(1:Window*6-1,1) = TIME(1:Window*6-1,1);
            TEMP = AA(:,1:2)*X(1:2);
            Beginning(1:Window*6-1,2) = TEMP(1:Window*6-1,1);
            TEMP = AA(:,3:4)*X(3:4);
            Beginning(1:Window*6-1,3) = TEMP(1:Window*6-1,1);
            TEMP = AA(:,5:6)*X(5:6);
            Beginning(1:Window*6-1,4) = TEMP(1:Window*6-1,1);
        elseif j==m-Window*12+1
            Ending(1:Window*6,1) = TIME(end-Window*6+1:end,1);
            TEMP = AA(:,1:2)*X(1:2);
            Ending(1:Window*6,2) = TEMP(Window*6+1:Window*6*2,1);
            TEMP = AA(:,3:4)*X(3:4);
            Ending(1:Window*6,3) = TEMP(Window*6+1:Window*6*2,1);
            TEMP = AA(:,5:6)*X(5:6);
            Ending(1:Window*6,4) = TEMP(Window*6+1:Window*6*2,1);     
        end        
        %%
        AmpPha(j,1) = TIME(Window*6+j-1,1);
        AmpPha(j,2) = norm(X(3:4));
        AmpPha(j,3) = norm(X(5:6));
        [ang] = GetPhase(X(3),X(4));
        AmpPha(j,4) = ang/pi*180;       %annual phase
        [ang] = GetPhase(X(5),X(6)) ;    
        AmpPha(j,5) = ang/pi*180;       %semi-annual phase         
    end 
end
TotalPreMSL.Monthly = [Beginning,zeros(Window*6-1,6);PreMSL;Ending,zeros(Window*6,6)];

TotalPreMSL.AmpPha = AmpPha;





