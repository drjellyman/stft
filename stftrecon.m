close all; clear all;

% Import audio file
[y,Fs] = audioread('/home/jellymdayl/Documents/phd/Matlab/2000_20140121-1125/50_male_speech_english_ch10_orth_2Y.flac');
n = length(y);
n_dom = [1:n]';


% Take full length fft
ydft = fft(y); % Two sided dft of y
ypsd = (1/n)*(abs(ydft(1:n/2+1)).^2); % Make one sided, take magnitude, square, and scale
f = (Fs/(2*(length(ypsd))))*[0:length(ypsd)-1]'; % Freq variable for plot
figure; plot(f,ypsd); legend('ypsd(f)'); xlabel('f (Hz)'); grid on; 


% Reconstruct signal
yhat = ifft(ydft);
e = norm(y-yhat);
% figure; plot(n_dom,y,n_dom,yhat,'--'); legend('y(t)','yhat(t)'); xlabel('f (Hz)'); grid on; 
% It appears to have perfect reconstruction, as expected.


% Now the same with STFT
N = 2^9; % Window length
nWin = [1:N]';
w = 0.5*(1-cos((2*pi*nWin)/(N-1)));
% w_sum = zeros(n,1); % Use w_sum to check the window sums to 1 over time
for k = 1:floor(n/(N/2)) -1
    ywDft(:,k) = fft(y(0.5*N*(k-1)+1 : 0.5*N*(k+1)) .* w); 
end

ywDftSum = sum(ywDft');
ywSumPsd = (1/N)*(abs(ywDftSum(1:N/2+1)).^2);
figure; plot(ywSumPsd); legend('ywSumPsd'); xlabel('n (samples)'); grid on; 

ywPsd = (1/N)*(abs(ywDft(1:N/2+1,:)).^2);
figure; imagesc(10*log10(flip(ywPsd,1))); xlabel('n (samples)'); 