close all; clear all;
Fs = 8e3; 
t = [0:1/Fs:1]';
f = 73;
s = sin(2*pi*f*t);

S = fft(s);
freq = [0:(1/(length(s)-1)):1]';

tShift = 5; % Shift in samples
SShifted = S.*exp(-j*2*pi*tShift*freq);
sShifted = ifft(SShifted);

figure; scatter(t,s,'.'); hold on; scatter(t,sShifted,'.'); legend('s','sShifted'); grid on;
  

% Ok what about upsampling then downsampling. Too simplistic? 
% This is a delay from my freqBF.m script -> 0.023597307572534, which with
% Fs = 8e3 means a delay of 377.55692116 samples.


