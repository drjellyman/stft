close all; 
clear all;
[y,Fs] = audioread('/home/jellymdayl/Documents/phd/Matlab/2000_20140121-1125/50_male_speech_english_ch10_orth_2Y.flac');
y = y - mean(y); % Remove any dc offset from y
y = y .* (1/max(y)); % Normalize y to peak at 1, in order to coincide with the window function
sig_duration = length(y)/Fs; % Length in seconds of imported sound file
t = [0:1/Fs:sig_duration-1/Fs]'; % Time domain variable t
N = 2^16; % Window length
n = [1:length(y)]'; % Time axis by sample number
sigma = 0.3; % Window variance? 
nWindowIter = floor(2*(length(y)/N) - 1)
for aaa = 1:nWindowIter
    w = exp(-0.5*((n-aaa*(N-1)/2)/(sigma*(N-1)/2)).^2); % Shifted Gaussian window function 
    X = fft(y(1:N).*w(1:N)); % Spectrum of the windowed signal
    X2 = abs(X/N); % 2 sided distribution of X
    X1(:,aaa) = 2*X2(1:N/2+1); % Chop off the other side of the spectrum
end
Xsum = sum(X1')';
Xtrue = abs(fft(y)/N);
Xtrue = 2*Xtrue(1:N/2+1);
norm(Xtrue-Xsum)

plot(t,y,t,w); legend('s(t)', 'w(t)'); xlabel('t (seconds)'); grid on;
figure; plot(Xsum); xlim([0, Fs/2]); hold on; plot(Xtrue); legend('Xsum','Xtrue');
% spectrogram(y, 2^7, 'yaxis');
