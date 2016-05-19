close all; clear all;

Fs = 24e3; % (Hz)
t = [0:1/Fs:2-1/Fs]'; % (s)
sf = 933; % (Hz)
s = 0.01 * sin(2*pi*sf*t);
N = length(s);

S2 = fft(s) / N; % 2 sided spectrum of s
f2 = [0:1/(N-1):1]'; % [0:Fs/(N-1):Fs]'; % 2 sided frequency domain
figure; plot(f2,abs(S2)); legend('magnitude(S2)');
figure; plot(f2,angle(S2)); legend('phase(S2)');

sHat2 = ifft(S2);

sau = 0;
Z = exp(-j*2*pi*[0:N-1]'*sau/N) .* S2; 
figure; plot(f2,abs(Z)); legend('magnitude(Z)');
figure; plot(f2,angle(Z)); legend('phase(Z)');
