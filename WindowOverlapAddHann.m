close all; clear all;

% Import audio files
[s1,Fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency

% Shorten the signals
NSig = 2^10;
s1 = s1(1:NSig); s2 = s2(1:NSig);

K = 2^8+1; % Window length in samples
w = 0.5*(1-cos((2*pi*[0:K-1]')/(K-1))); % Hann window for stft
L = floor(length(s1(:,1))/(K/2)) -1; % The number of windows used across the full length of the input signal

% Check window
wSum = zeros(NSig,1);
sp = 1 + floor(K/2) * [0:L]'; % sp = start points; the indexes that each consecutive window shouuld start at
fig1 = figure; hold on;
for l = 1:L
    wSum(sp(l):sp(l)+K-1) = wSum(sp(l):sp(l)+K-1) + w;
    
    % Plot the offset window
    wCheck(:,l) = zeros(NSig,1); 
    wCheck(sp(l) : sp(l)+K-1,l) = w;
    plot(wCheck,'--');
end

% Plot the result of overlap-add
plot(wSum); set(gca,'fontsize',18); xlim([0,1000]); ylim([-0.05,1.05]); grid on;
