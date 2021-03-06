close all
clear all

% check fft

Fs = 16e3;
t = [0 : 1/Fs : 1-1/Fs]';
n = length(t);
f1 = 1e3; 
f2 = 4.3e3;
s1 = sin(2*pi*f1*t);
s2 = 0.6*cos(2*pi*f2*t); 
noise = 0.1*randn(n,1);
s = s1 + s2;
snrTrue = snr(s,noise)
% sPowerRMS = rms(s)^2
sPowerRMS = ((mean(s.^2))^(0.5))^2 % power = (rms(s))^2
% noisePowerRMS = rms(noise)^2
noisePowerRMS = ((mean(noise.^2))^(0.5))^2
snrDayle = 10 * log10(sPowerRMS / noisePowerRMS)
s = s + noise;
S = fft(s);
S = abs(S/n); 
S = S(1:n/2+1);
S(2:end-1) = 2*S(2:end-1);

% stft 
sig_duration = length(s)/Fs; % Length in seconds of imported sound file
t = [0:1/Fs:sig_duration-1/Fs]'; % Time domain variable t
N = 2^9; % Window length
sample_n = [1:length(s)]'; % Time axis by sample number
sigma = 0.3; % Window variance? 
nWindowIter = floor(2*(length(s)/N) - 1)
for aaa = 1:nWindowIter
    w = exp(-0.5*((sample_n-aaa*(N-1)/2)/(sigma*(N-1)/2)).^2); % Shifted Gaussian window function 
    X = fft(s(1:N).*w(1:N)); % Spectrum of the windowed signal
    X2 = abs(X/N); % 2 sided distribution of X
    X1temp = X2(1:N/2+1); % Chop off the other side of the spectrum
    X1temp(2:end-1) = 2*X1temp(2:end-1);
    X1(:,aaa) = X1temp; % This conserves total power, 0 Hz and Nyquist freq are not doubled as they only occur once in the two sided spectrum
end
XSum = sum(X1')';
fScaled = [0:(Fs/2)/(N/2):length(S)]';

% Approximate snr from mixed signal
notch1 = designfilt('bandstopfir','PassbandFrequency1', 0.115, ...
     'StopbandFrequency1',0.1249, ...
     'StopbandFrequency2', 0.1251, ...
    'PassbandFrequency2', 0.135, ...
        'PassbandRipple1', 1, ...
    'StopbandAttenuation', 60, ...
        'PassbandRipple2', 1, ...
           'DesignMethod', 'equiripple');
% fvtool(notch1)
notch2 = designfilt('bandstopfir','PassbandFrequency1', 0.53, ...
     'StopbandFrequency1', 0.5373, ...
     'StopbandFrequency2', 0.5375, ...
    'PassbandFrequency2', 0.5425, ...
        'PassbandRipple1', 1, ...
    'StopbandAttenuation', 60, ...
        'PassbandRipple2', 1, ...
           'DesignMethod', 'equiripple');
%        fvtool(notch2)
% sNotched = filter(); 
sFilt = filter(notch1,s);
sFilt = filter(notch2, sFilt);
SFilt = fft(sFilt);
SFilt2 = abs(SFilt/n);
SFilt1 = SFilt2(1:n/2+1);
SFilt1(2:end-1) = 2*SFilt1(2:end-1);

sPlusNoisePowerRMS = ((mean(s.^2))^(0.5))^2
NoisePowerRMSApprox = ((mean(SFilt1.^2))^0.5)^2
SNRApprox = 10 * log10((sPlusNoisePowerRMS-NoisePowerRMSApprox)/(NoisePowerRMSApprox))
SNRApprox = 10 * log10((sPlusNoisePowerRMS)/(NoisePowerRMSApprox))


W = fft(w);
W = W(1:N/2+1);
W = (1/(Fs*N))*(abs(W).^2);
W(2:end-1) = 2*W(2:end-1);

figure; plot(t,s); legend('s(t)'); xlabel('t (seconds)'); grid on; 
figure; plot(S); legend('S(f)'); xlabel('f (Hz)'); grid on; 
figure; plot(fScaled,XSum); legend('XSum'); xlabel('f (Hz)'); grid on; 
figure; plot(SFilt1); legend('SFilt1'); xlabel('f (Hz)'); grid on; 
figure; plot(fScaled,10*log10(W)); legend('W'); xlabel('f (Hz)'); grid on;