close all; clear all;
Fs = 16e3;
t = [0:1/Fs:1-1/Fs]';
n = length(t);
s = sin(2*pi*0.3e3*t);
noise = 0.1*randn(n,1);
s = s + noise;

% Take fft of s
S = fft(s);
S = abs(S(1:n/2+1))/Fs;
S(2:end-1) = 2*S(2:end-1);
f = [0:Fs/2]';

% Take ifft of S

sound(s,Fs)

plot(f,S)