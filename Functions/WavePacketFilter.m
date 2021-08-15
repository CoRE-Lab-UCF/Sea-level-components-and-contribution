function [RecSig] = WavePacketFilter(signal,level,low_freq,high_freq,Hz)

% PROGRAM "WavePacketFilter"
% Tool to filter signal with wavelet decomposition.
% Written by Sida Li
% Date: 14/8/2021
%
% Input:
%       1. signal: a vector 
%       2. level:  the levels we want to decompose the signal
%       3. low_freq: low frequency cutoff
%       4. high_freq: high frequency cutoff
%       5. Hz: The sampling frequency， usuaally "1".
%            
% Output: 
%       1. RecSig: filtered signal. 

% signal 输入信号（向量）
% level 小波包分解层数
% low_freq  低频截止频率
% low_freq  高频截止频率
% Hz  采用频率
% RecSig 是滤波后的信号

% for example :  有一个信号，信号的采样 周期是一个月采样一个值，
% 想提取2-7年左右的信号，则：
%  [RecSig] = WavePacketFilter(signal,6,1/84,1/16,1) ;
% 这里，6代表小波包分解层数；最后一个参数 hz=1（我们心里明白，这个hz对应的单位就是 “月”）
% 1/84 表示低频截止频率， 12*7年 = 84 月； 1/18表示高频截止频率 12*1.5年 = 18月；
% 这里需要自己设计的是 level 、 low_freq、 low_freq
% 关于 level 的设计 ，对于具有11周期的太阳黑子信号，建议level取 7 或取8，或取9，均可试试
% 高低频的截止频率建议  low_freq = 1/(14年*12) = 1/168；   high_freq = 1/(8年*12) = 1/106
% 或者基于上述截止频率微调

wpt= wpdec(signal,level,'db45');  % coif5   dmey
nodes = get(wpt,'Tn');
N_cfs=length(nodes);%小波包系数个数
ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
nodes_ord=nodes(ord);%重排后的小波系数

step = 1:2^level;
freq = 0:Hz/2/2^level:Hz/2;
stop_num = [];start_num = [];
for j = 1:2^level
    if (freq(j) - low_freq)<=0 & (freq(j+1) - low_freq)>=0
        start_num = j;
    end
    
    if (freq(j) - high_freq)<=0 & (freq(j+1) - high_freq)>=0
        stop_num = j;
    end    
    
end
 
for i=1:length(freq)-1  
    rex3(:,i)=wprcoef(wpt,nodes_ord(i));  %?????????????        
end
RecSig = sum(rex3(:,start_num:stop_num),2);


