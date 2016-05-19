close all ; clear all ; 

Fs = 1;
s = cos(2*pi*0.01*[0:99]');
N = length(s);
Ts = 1/Fs;
Tmax = (N-1)*Ts;
t = 0:Ts:Tmax;
figure; plot(t,s);
f = -Fs/2:Fs/(N-1):Fs/2;
S = fft(s)/N;
S2 = fft(s)/N;
delay = [1,2,3,3.2,3.5]';
k = [0:N-1]';
for c = 1:5    
    A = exp(-j*2*pi*k*delay(c)/N);
    Sshifted(:,c) = A .* S;
end
figure; plot(t,ifft(Sshifted(:,1)),t,ifft(Sshifted(:,2)),t,ifft(Sshifted(:,3)),t,ifft(Sshifted(:,4)),t,ifft(Sshifted(:,5))); legend('delay=1','delay=2','delay=3','delay=3.2','delay=3.5');
set(gca,'fontsize',18)