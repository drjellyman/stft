% BF conjugate symmetric
close all; clear all;

% Import audio files
[s1,fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency
clear Fs2; 
% s1 = resample(s1,1,3);
% s2 = resample(s2,1,3);
% fs = fs/3;

s1 = zeros(length(s1),1);
% s2 = zeros(length(s2),1);

% Shorten the signals
NSig = 2^17; % Truncate signals to be a y*0.5 multiple of the window length, where y is an integer.
start = 29500;
s1 = s1(start:NSig-1+start); s2 = s2(start:NSig-1+start);
clear start;
s = [s1,s2]; clear s1; clear s2;

%% STFT 
K = 2^8+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.

% pad the source signals so the 1st half window doesn't distort the data
s1Padded = [zeros((K-1)/2,1);s(:,1);zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s(:,2);zeros((K-1)/2,1)];

% Take stft of both sources
[S1,L] = stft(s1Padded,K);
[S2,L2] = stft(s2Padded,K);
clear L2;

% Discard top half of the spectrum, dc, and fs/2.
% S1 = S1(2:(K-1)/2,:);
% S2 = S2(2:(K-1)/2,:);
% K1 = length(S1(:,1)); % K1 = the truncated number of frequency bins

%% Place sensors
NSources = length(s(1,:));
M = 8; % M = number of sensors
dz = 0.9; % ds = distance between sensors
zPos = ones(3,M);
zPos(1,:) = zPos(1,:).*([0:M-1]*dz+1);
sPos = [1,2,1 ; 8,2.5,1]';
figure; plot3(zPos(1,:),zPos(2,:),zPos(3,:),'o',sPos(1,:),sPos(2,:),sPos(3,:),'*'); grid on; 
xlabel('x'); ylabel('y'); zlabel('z'); xlim([0, 9]); ylim([0,5]); zlim([0,3]); 

% Find the minimum distance from each source to the array
for m = 1:M
    for ns = 1:NSources
        d(m,ns) = norm(sPos(:,ns)-zPos(:,m));
    end
end
[d,i] = min(d); % d contains the minimum distance, i contains the index of the closest sensor
zmean = mean(zPos');
for ns = 1:NSources
    a = sPos(:,ns) - zmean';
    b = zPos(:,ns) - zmean';
    phi(ns) = atan2(norm(cross(a,b)),dot(a,b));
    dzphi(:,ns) = [0:M-1]'*dz*cos(phi(ns))+d(ns);
end


%% Create atf for both sources
c = 343; % Speed of sound in m.s^-1
% kdom = [1:K]'*fs/(K-1);
% kdom = [0:K-1]'*fs/(K-1);
kdom = [0:(K-1)/2 , -(K-1)/2:-1]';
for m = 1:M
    A(:,m) = exp(-j*2*pi*kdom'*dzphi(m,1)/c)/dzphi(m,1);
    A2(:,m) = exp(-j*2*pi*kdom'*dzphi(m,2)/c)/dzphi(m,2);
end

%% Create observations
for l = 1:L
    for m = 1:M
        Z(:,l,m) = A(:,m).*S1(:,l)+A2(:,m).*S2(:,l);
    end
end

%% ifft of Z1 (sensor 1)
Z1 = squeeze(Z(:,:,1));
Z2 = squeeze(Z(:,:,2));
% Z1 = [zeros(1,L) ; Z1 ; zeros(1,L) ; (flipud(Z1)) ; zeros(1,L)];
% mySpectrogram(Z1)
z1 = myOverlapAdd(Z1);
z2 = myOverlapAdd(Z2);
figure; plot(myNormalize(z1)); hold on; plot(myNormalize(z2));

%% Filter and sum fixed bf
% W0 = A/||A||^2
for l = 1:L
    for m = 1:M
        W0(:,l,m) = A(:,m)/(norm(A(:,m))^2);
    end
end
% Yfbf = W0^H * Z
for l = 1:L
    for m = 1:M
        W0Z(:,l,m) = squeeze(W0(:,l,m)).*squeeze(Z(:,l,m));
    end
end
Yfbf = sum(W0Z,3);
%% ifft of Yfbf 
yfbf = myOverlapAdd(Yfbf);










