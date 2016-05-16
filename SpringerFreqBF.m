clear all; close all; 

% Import audio files
[s,Fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[ns,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency

% Truncate one of the files so they have equivalent lengths
if (length(s) > length(ns)) s = s(1:length(s));
elseif (length(s) < length(ns)) ns = ns(1:length(s));
end

% STFT
K = 2^10; % Window length ~ 11 ms @ 48 kHz. How long is too long? 
w = sqrt(0.5*(1-cos((2*pi*[1:K]')/(K-1)))); % Sqrt hann window for stft
L = floor(length(s(:,1))/(K/2)) -1; % The number of windows used across the full length of the input signal

% Calculate the stft of the source and the stationary intereferer
for l = 1:L
    S(:,l) = fft(s(0.5*K*(l-1)+1 : 0.5*K*(l+1)) .* w ); 
    NS(:,l) = fft(ns(0.5*K*(l-1)+1 : 0.5*K*(l+1)) .* w ); 
end

% Initialize vars
sDur = length(s(:,1))/Fs1; % sDur = source duration in seconds
sLen = length(s(:,1)); % N = number of samples in (truncated) source .wav files
velSnd = 334; % The velocity of sound in air (ms^-1)
M = 4; % The number of sensors
sTau = 0; % Delay between mics for s
A = exp(-j*2*pi*sTau*[1:M]'*(ones(K,1)./[1:K]')'); % A is meant to be my freq domain ATF, made of delays/phase shifts only
% A = exp(-j*2*pi*sTau*[0:M-1]/K);


% Phase shift the noise to locate in a different spatial region
NS4 = NS; NS4(:,:,2) = NS; NS4(:,:,3) = NS; NS4(:,:,4) = NS; 
nTau = 2.92e-5; % Delay between mics for s
for m = 1:M
    for l = 1:L
        Z(:,l,m) = S(:,l) + (exp(-j*2*pi*nTau*(m-1)*[1:K]') .* NS(:,l)); % Z is the freq domain input to the bf
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency domain BF

% Precompute simplifications based on A (assumed delay only ATF)
P = eye(4,4) - A*A'/(norm(A)^2);
F = A'/(norm(A)^2);
W = F; % Initialize weight vector
mu = 0.5; % Arbtrary, smaller is slower? 

% Iterate over time and improve W as you go
for l = 2:L
    for m = 1:M 
        ZW = Z(:,l,m).*W(:,m); % Component-wise multiply each freq domain input bin with its filter weight
    end
    Y(:,l) = sum(ZW,3);
    
    % Update weights
    for k = 1:K
        W(k,:) = P * (W(k,:)-mu*squeeze(Z(k,l,:))'*conj(Y(k,l)))' + F(k,:)';
    end
end

% istft using overlap add to get bf output y(t)
% Take ifft of each block
yBlocks = ifft(Y);
y = zeros(length(s(:,1)),1); % Initialize y, length is same as original signals

% Now overlap add each block 
for l = 1: L
    y(0.5*K*(l-1)+1 : 0.5*K*(l+1)) = y(0.5*K*(l-1)+1 : 0.5*K*(l+1)) + yBlocks(:,l).*w;
end


