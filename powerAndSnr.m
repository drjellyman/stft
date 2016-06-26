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
energy_s1 = sHs1;
energy_s2 = sHs2;
power_s1 = (1/N)*(sHs1);
power_s2 = (1/N)*(sHs2);
% energy_s2b = norm(s2)^2;
% energy_s1b = norm(s1)^2;
SNRdBt = 10*log10(power_s1/power_s2);

% Fourier domain
S1fft = fft(s1); 
S2fft = fft(s2);
SHS1fft = S1fft'*S1fft;
SHS2fft = S2fft'*S2fft;
energy_S1fft = (1/N)*SHS1fft;
energy_S2fft = (1/N)*SHS2fft;
power_S1fft = ((1/N)^2)*SHS1fft;
power_S2fft = ((1/N)^2)*SHS2fft;
SNRdBfft = 10*log10(power_S1fft/power_S2fft);

% Time-frequency domain
s1Padded = [zeros((K-1)/2,1);s1;zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s2;zeros((K-1)/2,1)];
[S1stft,L] = stft(s1Padded,K);
[S2stft,L] = stft(s2Padded,K);
% [S1b,L] = stft(s1,K);
% [S2b,L] = stft(s2,K);
S1stftHalf = S1stft(2:(K+1)/2-1,:); % Truncate to half spectrum
S2stftHalf = S2stft(2:(K+1)/2-1,:); % Truncate to half spectrum
S1stftHalfInclDCFS2 = S1stft(1:(K+1)/2,:); % Truncate to half spectrum
S2stftHalfInclDCFS2 = S2stft(1:(K+1)/2,:); % Truncate to half spectrum

% istft with two sided spectrum
s1stft = myOverlapAdd(S1stft);
s2stft = myOverlapAdd(S2stft);
s1stft = s1stft((K-1)/2+1:end-((K-1)/2));
s2stft = s2stft((K-1)/2+1:end-((K-1)/2));
sHs1stft = s1stft'*s1stft; 
sHs2stft = s2stft'*s2stft;
energy_s1stft = sHs1stft;
energy_s2stft = sHs2stft;
power_s1stft = (1/N)*(sHs1stft);
power_s2stft = (1/N)*(sHs2stft);
norm(s1stft-s1)
SNRdBistft = 10*log10(power_s1stft/power_s2stft);

% sum the energy/power in freq domain for 2sided stftwqqqqqqqqqaQ

energy_S1stftHalf = 0;
energy_S2stftHalf = 0;
power_S1stftHalf = 0;
power_S2stftHalf = 0;
for l=1:L
    tmp1 = S1stftHalf(:,l)'*S1stftHalf(:,l);
    tmp2 = S2stftHalf(:,l)'*S2stftHalf(:,l);
    energy_S1stftHalf = energy_S1stftHalf + (2/K)*(tmp1);
    energy_S2stftHalf = energy_S2stftHalf + (2/K)*(tmp2);
    power_S1stftHalf = power_S1stftHalf + (1/(L-1))*(2/K)^2*(tmp1);
    power_S2stftHalf = power_S2stftHalf + (1/(L-1))*(2/K)^2*(tmp2);
end
SNRdBstftHalf = 10*log10(power_S1stftHalf/power_S2stftHalf);

% istft with half spectrum
S1stftHalfRecon = [zeros(1,L);S1stftHalf;zeros(2,L);conj(flipud(S1stftHalf))];
S2stftHalfRecon = [zeros(1,L);S2stftHalf;zeros(2,L);conj(flipud(S2stftHalf))];
s1stftHalfRecon = myOverlapAdd(S1stftHalfRecon);
s2stftHalfRecon = myOverlapAdd(S2stftHalfRecon);
s1stftHalfRecon = s1stftHalfRecon((K-1)/2+1:end-((K-1)/2));
s2stftHalfRecon = s2stftHalfRecon((K-1)/2+1:end-((K-1)/2));
sHs1stftHalfRecon = s1stftHalfRecon'*s1stftHalfRecon; 
sHs2stftHalfRecon = s2stftHalfRecon'*s2stftHalfRecon;
energy_s1stftHalfRecon = sHs1stftHalfRecon;
energy_s2stftHalfRecon = sHs2stftHalfRecon;
power_s1stftHalfRecon = (1/N)*(sHs1stftHalfRecon);
power_s2stftHalfRecon = (1/N)*(sHs2stftHalfRecon);
norm(s1stftHalfRecon-s1)
SNRdBistftHalf = 10*log10(power_s1stftHalfRecon/power_s2stftHalfRecon);


% Try again with half spectrum including DC and fs/2
energy_S1stftHalfInclDCFS2 = 0;
energy_S2stftHalfInclDCFS2 = 0;
power_S1stftHalfInclDCFS2 = 0;
power_S2stftHalfInclDCFS2 = 0;
for l=1:L
    tmp1 = S1stftHalfInclDCFS2(:,l)'*S1stftHalfInclDCFS2(:,l);
    tmp2 = S2stftHalfInclDCFS2(:,l)'*S2stftHalfInclDCFS2(:,l);
%     energy_S1stftHalfInclDCFS2 = energy_S1stftHalfInclDCFS2 + (2/K)*(tmp1);
%     energy_S2stftHalfInclDCFS2 = energy_S2stftHalfInclDCFS2 + (2/K)*(tmp2);
    energy_S1stftHalfInclDCFS2 = energy_S1stftHalfInclDCFS2 + (2/(K+1))*(tmp1);
    energy_S2stftHalfInclDCFS2 = energy_S2stftHalfInclDCFS2 + (2/(K+1))*(tmp2);
    power_S1stftHalfInclDCFS2 = power_S1stftHalfInclDCFS2 + (1/(L-1))*(2/(K+1))^2*(tmp1);
    power_S2stftHalfInclDCFS2 = power_S2stftHalfInclDCFS2 + (1/(L-1))*(2/(K+1))^2*(tmp2);
end
SNRdBstftHalfInclDCFS2 = 10*log10(power_S1stftHalfInclDCFS2/power_S2stftHalfInclDCFS2);

% istft with half spectrum
S1stftHalfInclDCFS2Recon = [S1stftHalfInclDCFS2;conj(flipud(S1stftHalfInclDCFS2(2:end,:)))];
S2stftHalfInclDCFS2Recon = [S2stftHalfInclDCFS2;conj(flipud(S2stftHalfInclDCFS2(2:end,:)))];
s1stftHalfInclDCFS2Recon = myOverlapAdd(S1stftHalfInclDCFS2Recon);
s2stftHalfInclDCFS2Recon = myOverlapAdd(S2stftHalfInclDCFS2Recon);
s1stftHalfInclDCFS2Recon = s1stftHalfInclDCFS2Recon((K-1)/2+1:end-((K-1)/2));
s2stftHalfInclDCFS2Recon = s2stftHalfInclDCFS2Recon((K-1)/2+1:end-((K-1)/2));
sHs1stftHalfInclDCFS2Recon = s1stftHalfInclDCFS2Recon'*s1stftHalfInclDCFS2Recon; 
sHs2stftHalfInclDCFS2Recon = s2stftHalfInclDCFS2Recon'*s2stftHalfInclDCFS2Recon;
energy_s1stftHalfInclDCFS2Recon = sHs1stftHalfInclDCFS2Recon;
energy_s2stftHalfInclDCFS2Recon = sHs2stftHalfInclDCFS2Recon;
power_s1stftHalfInclDCFS2Recon = (1/N)*(sHs1stftHalfInclDCFS2Recon);
power_s2stftHalfInclDCFS2Recon = (1/N)*(sHs2stftHalfInclDCFS2Recon);
norm(s1stftHalfInclDCFS2Recon-s1)
SNRdBistftHalfInclDCFS2 = 10*log10(power_s1stftHalfInclDCFS2Recon/power_s2stftHalfInclDCFS2Recon);


EandP = [energy_s1,power_s1,energy_s2,power_s2;
        energy_S1fft,power_S1fft,energy_S2fft,power_S2fft;
        energy_s1stft,power_s1stft,energy_s2stft,power_s2stft;
        energy_S1stftHalf,power_S1stftHalf,energy_S2stftHalf,power_S2stftHalf;               
        energy_s1stftHalfRecon,power_s1stftHalfRecon,energy_s2stftHalfRecon,power_s2stftHalfRecon;
        energy_S1stftHalfInclDCFS2, power_S1stftHalfInclDCFS2, energy_S2stftHalfInclDCFS2, power_S2stftHalfInclDCFS2;
        energy_s1stftHalfInclDCFS2Recon, power_s1stftHalfInclDCFS2Recon, energy_s2stftHalfInclDCFS2Recon, power_s2stftHalfInclDCFS2Recon]
error = [EandP(1,:)-EandP(3,:);
        EandP(1,:)-EandP(4,:);
        EandP(1,:)-EandP(5,:);
        EandP(1,:)-EandP(6,:);
        EandP(1,:)-EandP(7,:)]
SNR = [SNRdBt;
       SNRdBfft;
       SNRdBistft;
       SNRdBstftHalf;
       SNRdBistftHalf;
       SNRdBstftHalfInclDCFS2;
       SNRdBistftHalfInclDCFS2]

% fid = fopen('stftenergyerror.txt', 'at');
% fprintf(fid, 'No padding, L, error = %3.8f, %3.8f, %3.8f, %3.8f\n', error);
% fclose(fid);