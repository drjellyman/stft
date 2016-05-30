close all;
clear all; 
% Haykin MVDR example
M = 5; % M = number of sensors
NSources = 2;
f = 200; % f = frequency of received signal (Hz)
fs = 16e3; 
ts = 1/fs;
N = 32e3;
tmax = (N-1)*ts;
t = [0:ts:tmax].';
soi = sin(2*pi*f*t);
interf = soi; 
c = 343; % c = speed of sound in air (m.s^-1)
lambda = c/f; % lambda = wavelength of received signal (m)
d = lambda/2; % d = distance between sensors (m)
phit = asin(0.2); % phit = phi target = source of interest angle
phii = asin(0); % phii = phi interferer = angle of interference signal
theta = pi*sin(phit); % theta is the angle of arrival
sTheta = exp(-j*[0:M-1].'*theta); % sTheta = the steering vector
mu = 1e-8; % mu = step size parameter

% create observations
% Calculate observation signals x_m, as a mixture of the two sources based
% on the distance between the source and the sensor
x = zeros(N, M);
NOrder = 3;
x = repmat(interf,1,M);
    for m = 1:M
        dphit = d*sin(phit);
        dt = dphit/c;
        ds = m*dt/ts;
        dsInt = floor(ds) - 1; % Split up the delay into an integer and a rational between 1 and 2
        dsSmall = ds - dsInt;         
        sTemp = [soi(dsInt+1:end) ; zeros(dsInt,1)]; % Delay the source by ssdsInt
        x(:,m) = x(:,m) + filter(lagrange(NOrder,dsSmall),1,sTemp) + 0.02*randn(N,1); % Delay the source by ssdsSmall, then weight and sum into the obervation signal
    end
% 
% figure; plot(x(:,1)); hold on; 
% plot(x(:,2));
% plot(x(:,3));
% plot(x(:,4));
% plot(x(:,5));

