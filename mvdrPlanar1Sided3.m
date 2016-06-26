close all; clear all;

% Import audio files
[s1,fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency

s1 = resample(s1,1,3);
s2 = resample(s2,1,3);
fs = fs/3;

s1 = zeros(length(s1),1);
s2 = zeros(length(s2),1);

% Shorten the signals
NSig = 2^16; % Truncate signals to be a y*0.5 multiple of the window length, where y is an integer.
start = 1;
s1 = s1(start:NSig-1+start); s2 = s2(start:NSig-1+start);
s = [s1,s2]; 

%% STFT 
K = 2^8+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.

% pad the source signals so the 1st half window doesn't distort the data
s1Padded = [zeros((K-1)/2,1);s(:,1);zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s(:,2);zeros((K-1)/2,1)];

% Take stft of both sources
[S1,L] = stft(s1Padded,K);
S2 = stft(s2Padded,K);

%% set up geometry
doa = [-0.2, 0];
d = 0.04; 
c = 343; 
dt = d*sin(doa)/c;
