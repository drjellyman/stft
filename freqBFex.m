close all; clear all;

f = 1e3
d = 3;
c = 343; 

fs = 3*f; % sampling frequency
Om = 2*pi*f/fs; % angular digital frequency
N = 50;
n = (0:N-1)'; % discrete time samples
del_s = d*fs/c; % some constants
mv = 0:M-1;
a = ones(M, 1)/M; % rectangular window
s_psi = sin((-90:90)*pi/180); % sine of incoming plane wave angle

jfi_b = -j*2*pi*del_s*sin(dir*pi/180)/N ; % beam direction parameter

% Rotational ( Fourier Transform ) operator
XFT = exp(jfi_b*(0:(N/2)-1)'*(0:M-1)) ;
for id = 1:length(s_psi) % FOR each INPUT direction
q = del_s*s_psi(id) ;
x = sin(Om*( n(:, ones(1, M)) + q*mv(ones(N, 1), :) ));
X = fft(x); % get the discrete fourier transform
X = X(1:N/2,:); % only interested in the single sided spectrum
Y_b = (X.*XFT)*a ; % Fourier transform of the beam "dir"
y_b = abs(ifft(Y_b)); % get back y_b(n)
P(id) = sum(y_b.^2); % get average power value of y_b
end;