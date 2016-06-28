clc

fs = 16000;
x = randn(fs*30,1); % 20 s of noise

K = 256;   % frame length          
H = 256/4; % frame hop
N = K;     % DFT size

scale = 2*sqrt(H)/sqrt((4*0.54^2+2*0.46^2)*(K+1)); % modified hamming window
w = scale*(0.54-0.46*cos(2*pi*(1:K+1)/(K+1)))';
w(end)=[];

overlap = K-H;
samples = floor((length(x)-overlap)/H); % number of samples

frames = zeros(K,samples);
for s=1:samples
   frames(:,s) = w.*x((1:K)+(s-1)*H); 
end

stft_x = fft(frames,N); % 2 sided stft

X = (K/H)*abs(stft_x/N).^2; % power normalized spectrogram (K/H accounts for overlap)

% total power
power_t = mean(x.^2)
power_f = mean(sum(X))

error = power_t - power_f

% the error exists because the sum of shifted windows is not constant at
% the start and end of the signal. In these regions the window reduces the 
% power, so power_t > power_f. The longer the signal, the less this effect.