%% Example code for dealing with conjugate symmetry. It seems easier to just 
% use the fftshift functionality to get it right. 

%% Working version
% clear all; close all;
% [s,fs] = audioread('273177__xserra__la-vaca-cega-eva.wav'); 
% s = s(1:15,1);
% K = length(s);
% S = fft(s);
% 
% Sshift = fftshift(S);
% Sshift = [0;Sshift(2:(K-1)/2);0;Sshift((K+1)/2+1:end-1);0];
% ssHat = ifft(ifftshift(Sshift))
% Shalf = Sshift(2:(K-1)/2);
% SHat = [0;Shalf;0;conj(Shalf(end:-1:1));0];
% 
% figure; plot(abs(Sshift));hold on; plot(abs(SHat),'*')
% figure; plot(angle(Sshift));hold on; plot(angle(SHat),'*')
% 
% shat = ifft(ifftshift(SHat))
% figure; plot(shat); hold on; plot(ssHat,'*')


%% Experimenting version
clear all; close all;
[s,fs] = audioread('273177__xserra__la-vaca-cega-eva.wav'); 
s = s(1:15,1);
K = length(s);
S = fft(s);
figure; plot(abs(S));
% sHat = ifft(S);

figure; plot(abs(S(2:(K+1)/2-1))); hold on; 
plot(abs(S(end:-1:(K+1)/2+2)),'*');
SHat = [0;S(2:(K+1)/2-1);0;0;conj(S((K+1)/2-1:-1:2))];
figure; plot(abs(SHat));
figure; plot(angle(SHat));
sHat = ifft(SHat);

% Sshift = fftshift(S);
% Sshift = [0;Sshift(2:(K-1)/2);0;Sshift((K+1)/2+1:end-1);0];
% ssHat = ifft(ifftshift(Sshift))
% Shalf = Sshift((K+1)/2+1:end-1);
% SHat = [0;Shalf;0;conj(Shalf(end:-1:1));0];
% 
% figure; plot(abs(ifftshift(Sshift)));hold on; plot(abs(SHat),'*')
% figure; plot(angle(ifftshift(Sshift)));hold on; plot(angle(SHat),'*')

% shat = ifft(ifftshift(SHat))
% figure; plot(shat); hold on; plot(ssHat,'*')
