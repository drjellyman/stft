clear all;
close all;

T = 1;  %length of sampling sequence in s
N = 64; %number of samples
M = 3; % number of periods per sequence
ts = T/N; %sample interval
fs = 1/ts %sampling frequency
tmax = (N-1)*ts;
t = 0:ts:tmax;
y = sin(2*pi*M*t);

fig01 = figure;
plot(t,y);
grid on;

%% We do the FT
Y = fft(y);

%% We create a frequency vector in natural order
% -fs/2, ..., 0, ... +fs/2
f =fftshift(( 0:(fs-1)) - fs/2);

%% Show Magnitude spectrum
% There shold be only two lines at -M and +M
figure;
plot(f,abs(Y),'o');
grid on;

%% Attempt at phase shift
% y/t) -> Y(w)
% y(t-t0) -> Y(w) * exp(-i*w*t0)
% Phase shift of pi/2 in frequncy domain is equavalent to as time shift
% of T/4 in time domain

Y = Y.*exp(-i*2*pi*f*T/4);

% Inverse FT
u = ifft(Y);

figure
hold on;
plot(t,real(u),'b-');
plot(t,real(y),'r-');
hold off;
grid;