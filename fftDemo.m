close all; clear all; 

% Proper use of fft, ifft, fftshift etc

% Import audio
[s,fs] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency
s = s(1:2^16); % truncate audio
N = length(s); % total number of samples
ts = 1/fs; % sample period
tmax = (N-1)*ts; % length of signal in time
t = 0:ts:tmax; % time domain 

% Plot time domain signal
fig01 = figure; plot(t,s);

% Find dft of s
f = 0:fs/(N-1):fs;
S = fft(s);

% Plot 2 sided frequency domain with Hz on x axis
fig02 = figure; subplot(2,1,1); plot(f,abs(S)); title('Two sided magnitude spectrum');
subplot(2,1,2); plot(f,angle(S)); title('Two sided phase spectrum'); xlabel('Frequency (Hz)'); 

% Plot 2 sided frequency domain with normalized frequency axis
fn = 0:1/(N-1):1;
fig03 = figure; subplot(2,1,1); plot(fn,abs(S)); title('Normalized two sided magnitude spectrum');
subplot(2,1,2); plot(fn,angle(S)); title('Normalized two sided phase spectrum'); xlabel('Frequency'); 

% Plot 1 sided frequency domain with normalized frequency axis
fnOneSide = fn(1:N/2+1);
SOneSide = S(1:N/2+1);
fig03 = figure; subplot(2,1,1); plot(fnOneSide,abs(SOneSide)); title('Normalized one sided magnitude spectrum');
subplot(2,1,2); plot(fnOneSide,angle(SOneSide)); title('Normalized one sided phase spectrum'); xlabel('Frequency'); 

% By chopping the spectrum in half, we have reduced the energy represented
% so it is common to double the magnitude. This excludes the DC and nyquist
% values. 
SOneSide2(2:end-1) = 2*SOneSide2(2:end-1);