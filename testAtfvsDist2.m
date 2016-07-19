close all; clear;
Nsrcs = 26; % Nsrcs = number of sources

s = audioread(strcat('/audio/','422-122949-0013.flac')); 
fs = 16e3;
K = 2^9+1; % K = window length in samples, and the number of frequency bins
Khalf = (K-1)/2-1;
tls = 5; % tls = target length in seconds
tl = tls*fs-mod(tls*fs,K-1); % tl = target length in samples, adjusted for window length and sampling frequency
s = s(1:tl); % truncation
ss = s(1:end-1);
s =  [zeros((K-1)/2,1);s;zeros((K-1)/2,1)]; % zero padding
[S,L] = stft(s,K);
S = S(2:(K-1)/2,:); % half spectrum

dist = 11;% 32.8;%27.4;

fdom = (fs/(K-1)) * (1:Khalf)';
c = 343; % c = speed of sound in m.s^-1
MaxAbsCorr = zeros(length(dist),1);
ExPwr = zeros(Khalf,length(dist));
for d=1:length(dist)
    A = exp(-1i*2*pi*fdom*dist(d)/c) ;%/ (4*pi*dist(d)); 
%     ExPwr(:,d) = fdom*dist(d)/c;
    X = zeros(Khalf,L); 
    for l = 1:L
        X(:,l) = X(:,l) + A .* S(:,l);
    end
    X = [zeros(1,L);X;zeros(2,L);conj(flipud(X))];
    x = myOverlapAdd(X);
    MaxAbsXcorr(d) = max(abs(xcorr(x,s)));
end
% figure; semilogy(dist,MaxAbsXcorr); grid on;

%% 
% figure; hold on; grid on; plot(ExPwr(end,:))

%% Go again but use full length fft for freq dom phase shift
KK = length(ss);
SS = fft(ss); 
SS = SS(2:(KK-1)/2); % Half spectrum
fdelta = fs/(KK-1);
ffdom = fdelta*(1:(KK-1)/2-1);
dist = 34.481999999999;%[0.1:0.001:100];%5.6;% 32.8;%27.4;
for d=1:length(dist)
    AA = exp(-1i*2*pi*ffdom*dist(d)/c) ;%/ (4*pi*dist(d)); 
%     ExPwr(:,d) = fdom*dist(d)/c;
    XX = zeros(KK,1); 
    XX = AA.' .* SS;
    
    XX = [0;XX;0;0;conj(flipud(XX))];
    xx = ifft(XX);
    MaxAbsXcorr2(d) = max(abs(xcorr(xx,ss)));
end
figure; semilogy(dist,MaxAbsXcorr2); grid on;

%% 
figure; plot(xx,'.'); hold on ; grid on; plot(ss,'o');
