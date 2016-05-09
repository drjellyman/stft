close all; clear all;

% Import audio file
[y,Fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac');
n = length(y);
n_dom = [1:n]';
t_length = length(y)/Fs;


% Take full length fft
ydft = fft(y); % Two sided dft of y
ypsd = (1/n)*(abs(ydft(1:n/2+1)).^2); % Make one sided, take magnitude, square, and scale
f = (Fs/(2*(length(ypsd))))*[0:length(ypsd)-1]'; % Freq variable for plot
figure; plot(f,ypsd); legend('ypsd(f)'); xlabel('f (Hz)'); grid on; 


% Reconstruct signal
yHat1 = ifft(ydft);
e1 = norm(y-yHat1)
% figure; plot(n_dom,y,n_dom,yhat,'--'); legend('y(t)','yhat(t)'); xlabel('f (Hz)'); grid on; 
% It appears to have perfect reconstruction, as expected.


% Now the same with STFT
N = 2^9; % Window length
nWin = [1:N]';
w = sqrt(0.5*(1-cos((2*pi*nWin)/(N-1))));
w_sum = zeros(n,1); % Use w_sum to check the window sums to 1 over time
num_win = floor(n/(N/2)) -1;
for k = 1:num_win
    ywDft(:,k) = fft(y(0.5*N*(k-1)+1 : 0.5*N*(k+1)) .* w); 
    w_sum(0.5*N*(k-1)+1 : 0.5*N*(k+1)) = w_sum(0.5*N*(k-1)+1 : 0.5*N*(k+1)) + w.^2;
end

% Plot the summed window 
figure; plot(w_sum); legend('w\_sum'); xlabel('n (samples)'); grid on; 

% Take ifft of each block
ywIft = ifft(ywDft);

% Now overlap add each block 
yHat2 = zeros(n,1);
for k = 1: floor(n/(N/2)) -1
    yHat2(0.5*N*(k-1)+1 : 0.5*N*(k+1)) = yHat2(0.5*N*(k-1)+1 : 0.5*N*(k+1)) + ywIft(:,k).*w;
end

% Check difference between original y and reconstructed yHat
e2 = norm(y-yHat2)
e2_minus_e1 = e2-e1

% Plot the two signals y and yHat
figure; plot(y); hold on; plot(yHat2); legend('y','yHat2'); xlabel('n (samples)'); 

% Plot the spectrogram
ywPsd = (1/N)*(abs(ywDft(1:N/2+1,:)).^2);
time = linspace(0, t_length, length(ywPsd(1,:)))';
freq = linspace(0, Fs/2, length(ywPsd(:,1)))';
figure; %subplot(2,1,1); 
imagesc(time,freq,10*log10(ywPsd)); xlabel('time (s)'); ylabel('frequency (Hz)'); axis xy; colorbar;
%subplot(2,1,2); plot(n_dom/Fs,y); xlabel('time (s)'); ylabel('amplitude'); axis tight;

% Comparison with matlab stft
figure; spectrogram(y,w,'yaxis'); 