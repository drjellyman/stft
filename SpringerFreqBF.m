clear all; close all; 

% Import audio files
[s1,Fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac');
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac');

% Truncate one of the files so they have equivalent lengths
if (length(s1) > length(s2)) s1 = s1(1:length(s1));
elseif (length(s1) < length(s2)) s2 = s2(1:length(s1));
end
s = [s1, s2]; clear s1; clear s2;

% Initialize vars
sLenT = length(s1)/Fs1; % sLenT = source duration in seconds
N = length(s1); % N = number of samples in source .wav files
velSnd = 334; % The velocity of sound in air (ms^-1)
A = zeros(4,4); % freq domain ATF
K = 4; 
M = 4; 
tau = 2.92e-5;
for k = 1:K
    for m = 0:M
        A(k,m+1) = exp(-j*2*pi*tau*m/k);
    end
end
