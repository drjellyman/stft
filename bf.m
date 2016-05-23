% close all
% clear all
% clc
% m = 8;
% signal=(sin(20* pi/180));
% ad = exp (-1j * pi *[0: m-1]'*signal);% array response vector in the direction of desired signal. expect the direction of the array response vector
% wop = ad;
% thetas = [-90:90];
% tm = thetas * pi/180;
% am = exp (-1j * pi * [0: m-1]'* (sin (tm)));
% A = abs (wop'* am);% array response array response
% A = A / max (A);
% figure, polar (tm, A)
% A = 10 * log10 (A);% log figure log plot
% title ('Normalized magnitude response array polar diagram, eight array elements')
% figure, plot (thetas, A);
% title ('eight microphones');
% xlabel ('angle[degrees]');
% ylabel ('Normalized Beam Power[dB]');
% grid on

close all; clear all;
t = (0:999)'/1000;
s = sin(2*pi*t);
ad = exp(-1i*pi*(0:7)*sin(50*pi/180));
x = s*ad;
y = x*ad';
figure; subplot(211),plot(t,s),subplot(212),plot(t,real(y));
figure; plot(real(x(:,1))); hold on; plot(real(x(:,2)));plot(real(x(:,3)));