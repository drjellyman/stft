%% Check getting the power in the time and frequency domains are equal!

close all; clear all; 

% Import target audio
[s,fs] = audioread('273177__xserra__la-vaca-cega-eva.wav'); 
s = s(:,1);
S = fft(s);
N = length(s);

% Calculate equivalent ENERGY in both domains
sEn = (s'*s) 
SEn = (1/length(S))*(S'*S) 

% We can also normalize by N or sqrt(N)
sEn = (s'*s) / sqrt(N)
SEn = (1/length(S))*(S'*S) / sqrt(N)

% Power is Energy/time
sPow = sEn/(N/fs)
SPow = SEn/(N/fs)