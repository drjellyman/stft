close all
clear all

% check fft

Fs = 16e3;
t = [0 : 1/Fs : 1-1/Fs]';
n = length(t);
f1 = 1e3; 
f2 = 4.3e3;
s1 = sin(2*pi*f1*t);
s2 = 0.6*cos(2*pi*f2*t); 
s3 = 0.1*randn(n,1);
s = s1 + s2 + s3; 
S = fft(s);
S = abs(S/n); 
S = 2*S(1:n/2);

% stft 




figure; plot(t,s); legend('s(t)'); xlabel('t (seconds)'); grid on; 
figure; plot(S); legend('S(f)'); xlabel('f (Hz)'); grid on; 