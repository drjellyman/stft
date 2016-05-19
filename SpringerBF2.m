close all; clear all;

% Closed form LCMV BF from Springer (eq 47.7)
% W_LCMV(k,l) = (Z(k,l)Z^H(k,l))^-1 A(k) / (A^H(k)(Z(k,l)Z^H(k,l))^-1 A(k))

% Import audio files
[s1,Fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency

% Shorten the signals
NSig = 2^12;
s1 = s1(1:NSig); s2 = s2(1:NSig);

%% STFT
K = 2^8+1; % Window length in samples
[S1, L] = stft(s1,K);
S2 = stft(s2,K);

% Have a look
mySpectrogram(S1);

%% Recreate signal
sHat = myOverlapAdd(S1);

