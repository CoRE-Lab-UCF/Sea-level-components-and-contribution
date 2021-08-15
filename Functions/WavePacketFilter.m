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
%       5. Hz: The sampling frequency�� usuaally "1".
%            
% Output: 
%       1. RecSig: filtered signal. 

% signal �����źţ�������
% level С�����ֽ����
% low_freq  ��Ƶ��ֹƵ��
% low_freq  ��Ƶ��ֹƵ��
% Hz  ����Ƶ��
% RecSig ���˲�����ź�

% for example :  ��һ���źţ��źŵĲ��� ������һ���²���һ��ֵ��
% ����ȡ2-7�����ҵ��źţ���
%  [RecSig] = WavePacketFilter(signal,6,1/84,1/16,1) ;
% ���6����С�����ֽ���������һ������ hz=1�������������ף����hz��Ӧ�ĵ�λ���� ���¡���
% 1/84 ��ʾ��Ƶ��ֹƵ�ʣ� 12*7�� = 84 �£� 1/18��ʾ��Ƶ��ֹƵ�� 12*1.5�� = 18�£�
% ������Ҫ�Լ���Ƶ��� level �� low_freq�� low_freq
% ���� level ����� �����ھ���11���ڵ�̫�������źţ�����levelȡ 7 ��ȡ8����ȡ9����������
% �ߵ�Ƶ�Ľ�ֹƵ�ʽ���  low_freq = 1/(14��*12) = 1/168��   high_freq = 1/(8��*12) = 1/106
% ���߻���������ֹƵ��΢��

wpt= wpdec(signal,level,'db45');  % coif5   dmey
nodes = get(wpt,'Tn');
N_cfs=length(nodes);%С����ϵ������
ord=wpfrqord(nodes);  %С����ϵ�����ţ�ord�����ź�С����ϵ���������ɵľ�����3��ֽ��[1;2;4;3;7;8;6;5]
nodes_ord=nodes(ord);%���ź��С��ϵ��

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


