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
wSum = zeros(length(t),1);
for aaa = 1:nWindowIter
    w = exp(-0.5*((n-aaa*(N-1)/2)/(sigma*(N-1)/2)).^2); % Shifted Gaussian window function 
    X = fft(y(1:N).*w(1:N)); % Spectrum of the windowed signal
    X2 = abs(X/N); % 2 sided distribution of X
    X1(:,aaa) = 2*X2(1:N/2+1); % Chop off the other side of the spectrum
    wSum = wSum + w;
end
Xsum = sum(X1')';
Xtrue = abs(fft(y)/length(y));
Xtrue = Xtrue(1:length(y)/2+1);
Xtrue(2:end-1) = 2*Xtrue(2:end-1);


plot(t,y,t,w); legend('s(t)', 'w(t)'); xlabel('t (seconds)'); grid on;
figure; plot(Xsum);  hold on; plot(Xtrue); legend('Xsum','Xtrue');
% spectrogram(y, 2^7, 'yaxis');
figure; plot(wSum); legend('wSum');

% So it would appear that the sum of the windows is not 1, which would seem
% to be a problem for reconstruction? Would it? Try anyways???

% Lets check out the spectrum of the window
w = exp(-0.5*((n(1:N)-(N-1)/2)/(sigma*(N-1)/2)).^2); % Shifted Gaussian window function 
W = fft(w);
W = abs(W/N);
W = W(1:N/2+1);
W(2:end-1) = 2*W(2:end-1);
figure; plot(W); legend('W');
figure; plot(w); legend('w');
% So it pretty much has only very low content. Which leads to the idea of
% low pass filtering the signal in order to window ? f

% One thing I can do is scale the summed spectrum to compare with the true
% spectrum:


% Let's check that (sum(window coefficients))^2 = 1
checkw = (sum(w))^2
