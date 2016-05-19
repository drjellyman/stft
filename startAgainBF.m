close all; clear all;
    
Fs = 48e3; % sampling rate
t = [0:1/Fs:1-1/Fs]';
s = 0.01*sin(2*pi*1.37e3*t); % source signal
as = 0.05; % as = array spacing (m)
sv = 343; % sv = sound velocity (m.s^-1)


aoa = 5*(2*pi/360);
tau = (as*sin(aoa))/(sv); % tau is the delay between microphones for s (s)
sau = tau * Fs; % sau is tau converted to samples
N = length(s);
S = fft(s);
Smag = abs(S);
Smag = Smag(1:N/2+1);
Smag(2:end-1) = 2 * Smag(2:end-1);


% A = exp(-j*2*pi*[0:3]'*sau);
% Z = (A*S')';
% for m = 1:4
%     z(:,m) = ifft(Z(:,m));
% end
% plot(real(z(:,1))) ; hold on ; plot(real(z(:,2))); hold on ;
% plot(real(z(:,3))); hold on; plot(real(z(:,4))); legend('1','2','3','4');

k = [0:1/(length(S)):1]';
kmag = [0:1/(length(Smag)-1):1]';
% figure; plot(kmag,Smag);
% figure; plot(kmag*Fs/2,Smag);

for aaa = 1:length(k) -1
    Z(aaa) = exp(j*2*pi*k(aaa)*24000/N) * S(aaa);
end
% figure; plot(abs(Z));
% figure; plot(abs(fftshift(Z)));
zHat = ifft(Z);
figure; plot(t,s); hold on; plot(t,imag(zHat)); legend('s','zHat');

% SS = fft(s);
% ss = ifft(SS);
% figure; plot(ss);