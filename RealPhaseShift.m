close all; clear all; 

% How to phase shift a signal and get a real signal out
% ref: http://dsp.stackexchange.com/questions/10250/how-to-circularly-shift-a-signal-by-a-fraction-of-a-sample

a = [1,2,3,3,4,4,3,3,2]';
N = length(a);
% k = [0:N-1];

% If you want the shifted output of the IFFT to be real, the phase 
% twist/rotation in the frequency domain has to be conjugate 
% symmetric, as well as the data. This can be accomplished by 
% adding an appropriate offset to your complex exp()'s exponent, 
% for the given phase slope, so that the phase of the upper (or 
% negative) half, modulo 2 Pi, mirrors the lower half in the FFT 
% aperture. The complex exponential shift function can also be made 
% conjugate symmetric by indexing it from -N/2 to N/2 with a phase 
% of zero at index 0.
k = [0:4,-4:-1];
A = fft(a);
p = exp(-1j*2*pi*2.21*k/N)';
AShift = A .* p;
aHat = ifft(AShift);
plot(real(aHat)); hold on; plot(a,'--');

% lets look at our exp for both cases:
k = [0:N-1];
figure; plot([-4:4],real(fftshift(exp(-1j * 2*pi * 0.5/N * k))));

k = [0:4,-4:-1];
figure; plot([-4:4],real(fftshift(exp(-1j * 2*pi * 0.5/N * k))));

