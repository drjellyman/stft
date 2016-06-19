%% adjust input audio levels

close all; clear all; 

% Import audio files
[s1,fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s1 = source, Fs1 = sampling frequency
[s2,fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % s2 = noise stationary, Fs2 = sampling frequency
[s3,fs3] = audioread('317354__speedenza__shall-i-compare-thee-voice.wav'); % s1 = source, Fs1 = sampling frequency
[s4,fs4] = audioread('273177__xserra__la-vaca-cega-eva.wav'); % s2 = noise stationary, Fs2 = sampling frequency
% Make them longer for longer adaptation
s1 = [s1];
s2 = [s2];
s3 = [s3];
s4 = [s4(:,1)];

% Downsample to 16k
fsd = 16e3; % fsd = desired sampling frequency
s1 = resample(s1,1,fs1/fsd);
s2 = resample(s2,1,fs2/fsd);
s3 = resample(s3,1,fs3/fsd);
s4 = resample(s4,1,fs4/fsd);
fs = fsd;

% Shorten the signals
K = 2^9+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.
lmin = min([length(s1),length(s2),length(s3),length(s4)]);
ld = 10*fs-mod(10*fs,K-1); % ld = desired length = seconds * seconds^-1
if lmin < ld
    s1 = s1(1:lmin); s2 = s2(1:lmin); s3 = s3(1:lmin); s4 = s4(1:lmin);
else
    s1 = s1(1:ld); s2 = s2(1:ld); s3 = s3(1:ld); s4 = s4(1:ld);
end
s = [s1,s2,s3,s4];
srms = rms(s)
srmsinv = 1./srms;
s = s * diag(srmsinv);
srmsPostNorm = rms(s)

figure; plot(s(:,1)); hold on; plot(s(:,2)); plot(s(:,3)); plot(s(:,4));