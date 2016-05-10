close all; clear all;
Fs = 8e3; 
t = [0:1/Fs:1]';
f = 73;
s = sin(2*pi*f*t);

S = fft(s);
freq = [0:(1/(length(s)-1)):1]';

tShift = 5; % Shift in samples
SShifted = S.*exp(-j*2*pi*tShift*freq);
sShifted = ifft(SShifted);

figure; scatter(t,s,'.'); hold on; scatter(t,sShifted,'.'); legend('s','sShifted'); grid on;

% There is a problem with using non integer delays, where the delay becomes
% incorrect and the amplitude of the output is reduced. Something to do
% with circular convolution?

% Ok what about upsampling then downsampling. Too simplistic? 
% This is a delay from my freqBF.m script -> 0.023597307572534, which with
% Fs = 8e3 means a delay of 377.55692116 samples.
td = 0.023597307572534; % td = time delay
sd = 1.5%td/(1/Fs); % sd = sample delay

N = length(S);
for k = 1:N
    if k < N-k
        SShifted2(k) = S(k) * exp(2*pi*i * sd*k/N);
    elseif k > N-k
        SShifted2(k) = S(k) * exp(2*pi*i * sd*(k-N)/N);
    elseif k == N/2
        SShifted2(k) = S(k) * cos(2*pi * sd/2);
    end
end
sShifted2 = ifft(SShifted2);

figure; scatter(t,s,'.'); hold on; scatter(t,sShifted2,'.'); legend('s','sShifted'); grid on;



sPower = sum(s.^2)
sShiftedPower = sum(sShifted.^2)
sShifted2Power = sum(sShifted2.^2)



% 
%                 exp(2*pi*i * s*k/N) for k < N - k,
% by
%                 exp(2*pi*i * s*(k-N)/N) for k > N - k,
% and by
%                 cos(2*pi * s/2) for k = N - k = N/2.

