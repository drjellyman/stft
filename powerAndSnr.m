close all; clear

% Import target audio
[s1,fs] = audioread('317354__speedenza__shall-i-compare-thee-voice.wav'); % s1 = source, Fs1 = sampling frequency
[s2,fs] = audioread('273177__xserra__la-vaca-cega-eva.wav'); 
K = 2^9+1; % K = window length in samples, and the number of frequency bins
tls = 10; % tls = target length in seconds
tl = tls*fs-mod(tls*fs,K-1); % tl = target length in samples, adjusted for window length and sampling frequency
s1 = s1(1:tl,1);
s2 = s2(1:tl,1);
N = length(s1);

% Time domain
sHs1 = s1'*s1; 
sHs2 = s2'*s2;
energy_s1a = sHs1;
energy_s1b = norm(s1)^2;
energy_s2a = sHs2;
energy_s2b = norm(s2)^2;
power_s1 = (1/N)*(sHs1);
power_s2 = (1/N)*(sHs2);

% Fourier domain
S1a = fft(s1); 
S2a = fft(s2);
SHS1a = S1a'*S1a;
SHS2a = S2a'*S2a;
energy_S1a = (1/N)*SHS1a;
energy_S2a = (1/N)*SHS2a;
power_S1a = ((1/N)^2)*SHS1a;
power_S2a = ((1/N)^2)*SHS2a;

% Time-frequency domain

s1Padded = [zeros((K-1)/2,1);s1;zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s2;zeros((K-1)/2,1)];
[S1b,L] = stft(s1Padded,K);
[S2b,L] = stft(s2Padded,K);
% [S1b,L] = stft(s1,K);
% [S2b,L] = stft(s2,K);
S1bhalf = S1b(2:(K+1)/2-1,:); % Truncate to half spectrum
S2bhalf = S2b(2:(K+1)/2-1,:); % Truncate to half spectrum

energy_S1b = 0;
energy_S2b = 0;
power_S1b = 0;
power_S2b = 0;
for l=1:L
    tmp1 = S1bhalf(:,l)'*S1bhalf(:,l);
    tmp2 = S2bhalf(:,l)'*S2bhalf(:,l);
    energy_S1b = energy_S1b + (2/K)*(tmp1);
    energy_S2b = energy_S2b + (2/K)*(tmp2);
    power_S1b = power_S1b + (1/(L-1))*(2/K)^2*(tmp1);
    power_S2b = power_S2b + (1/(L-1))*(2/K)^2*(tmp2);
end

% istft
S1bhalfrecon = [zeros(1,L);S1bhalf;zeros(2,L);conj(flipud(S1bhalf))];
S2bhalfrecon = [zeros(1,L);S2bhalf;zeros(2,L);conj(flipud(S2bhalf))];
s1recon = myOverlapAdd(S1bhalfrecon);
s2recon = myOverlapAdd(S2bhalfrecon);
s1recon = s1recon((K-1)/2+1:end-((K-1)/2));
s2recon = s2recon((K-1)/2+1:end-((K-1)/2));
sHsrecon1 = s1recon'*s1recon; 
sHsrecon2 = s2recon'*s2recon;
energy_s1c = sHsrecon1;
energy_s2c = sHsrecon2;
power_s1c = (1/N)*(sHsrecon1);
power_s2c = (1/N)*(sHsrecon2);



EandP = [energy_s1a,power_s1,energy_s2a,power_s2;
        energy_S1a,power_S1a,energy_S2a,power_S2a;
        energy_S1b,power_S1b,energy_S2b,power_S2b;
        energy_s1c,power_s1c,energy_s2c,power_s2c]
error = EandP(1,:)-EandP(3,:)

% fid = fopen('stftenergyerror.txt', 'at');
% fprintf(fid, 'No padding, L, error = %3.8f, %3.8f, %3.8f, %3.8f\n', error);
% fclose(fid);