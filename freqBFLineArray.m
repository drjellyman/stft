% Randomly place an arbitrary number of sensors within a 10 x 5 x 3 m room.
% Randomly place two audio sources within the same room. 

close all; clear all;

% Variables
NSensors = 4; % The number of sensors
NSources = 2; % The number of sources
Fs = 16e3; % The sample rate at the sensors
roomDim = [10, 5, 3]';
velSnd = 334; % The velocity of sound in air (ms^-1)

% Import audio files
[s1,Fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac');
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac');

% Check sample rates of the audio files
if Fs1 ~= Fs2 fprintf('The audio files have different sample rates\n'); end

% Truncate one of the files so they have equivalent lengths
if (length(s1) > length(s2)) s1 = s1(1:length(s1));
elseif (length(s1) < length(s2)) s2 = s2(1:length(s1));
end
sDur = length(s1)/Fs1; % sDur = source duration in seconds
s = [s1, s2];
clear s1; clear s2;

% Randomly place sensors and sources
% sensorLocation = diag(roomDim)*rand(3,NSensors); % Each column of sensor_locations gives x,y,z coordinates of a sensor.
% sourceLocation = diag(roomDim)*rand(3,NSources);
sensorLocation = [5,1,1 ; 5.01,1,1 ; 5.02,1,1 ; 5.03,1,1]'; % linear array with < 0.5*wavelength spacing @ 16kHz
sourceLocation = [5.015,3,1 ; 7,4.5,2]';

% Calculate observation signals x_m, as a mixture of the two sources based
% on the distance between the source and the sensor
NSamples = sDur/(1/Fs);
x = zeros(NSamples, NSensors);
NOrder = 3;
for a = 1:NSources
    for b = 1:NSensors
        ssd = norm(sensorLocation(:,b) - sourceLocation(:,a)); % ssd = source sensor distance (m)
        ssds = (ssd/velSnd) / (1/Fs1); % Delay between the sensor and source in samples
        ssdsInt = floor(ssds) - 1; % Split up the delay into an integer and a rational between 1 and 2
        ssdsSmall = ssds - ssdsInt; 
        weight = 1/(ssd+1)^2; % Calculate the weight based on the source to sensor distance
        sTemp = [s(ssdsInt+1:end,a) ; zeros(ssdsInt,1)]; % Delay the source by ssdsInt
        x(:,b) = x(:,b) + weight * resample(filter(lagrange(NOrder,ssdsSmall),1,sTemp),1,(Fs1/Fs)); % Delay the source by ssdsSmall, then weight and sum into the obervation signal
    end
end

% Now STFT
N = 2^9; % Window length
nWin = [1:N]';
w = sqrt(0.5*(1-cos((2*pi*nWin)/(N-1))));
% w_sum = zeros(n,1); % Use w_sum to check the window sums to 1 over time
numWin = floor(length(x(:,1))/(N/2)) -1;
for a = 1:NSensors
    for k = 1:numWin
        X(:,k,a) = fft(x(0.5*N*(k-1)+1 : 0.5*N*(k+1), a) .* w ); 
    end
end

% % Plot the spectrogram
% xwPsd = (1/N)*(abs(xwDft(1:N/2+1,:)).^2);
% time = linspace(0, sDur, length(xwPsd(1,:)))';
% freq = linspace(0, Fs/2, length(xwPsd(:,1)))';
% figure; %subplot(2,1,1); 
% imagesc(time,freq,10*log10(xwPsd)); xlabel('time (s)'); ylabel('frequency (Hz)'); axis xy; colorbar;

H = ones(N,numWin,NSensors);
Y = sum(X .* H, 3);

% Ok looks like I need to rewrite the inverse with the larger dimensions



% Calculate y[n] using istft
% Take ifft of each block
yBlocks = ifft(Y);

y = zeros(length(x(:,1)),1); % Initialize y, length is same as original signals

% Now overlap add each block 

for k = 1: numWin
    y(0.5*N*(k-1)+1 : 0.5*N*(k+1)) = y(0.5*N*(k-1)+1 : 0.5*N*(k+1)) + yBlocks(:,k).*w;
end

%http://download.springer.com/static/pdf/514/bok%253A978-3-540-49127-9.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Fbook%2F10.1007%2F978-3-540-49127-9&token2=exp=1463094477~acl=%2Fstatic%2Fpdf%2F514%2Fbok%25253A978-3-540-49127-9.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Fbook%252F10.1007%252F978-3-540-49127-9*~hmac=307acecabf2e784cc4b023fde23fcfffa5cd6ac6abf27fe8ddd0a4ce65691903

% Freq domain Frost algorithm ref springer speech processing handbook
% min W^H(k,l)X(k,l)X^H(k,l)W(k,l) s.t. W^H(k,l)A(k)=1
% where W is the bf weights, X is the signal post stft, and A is the atf.



