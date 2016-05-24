close all; clear all;

% Closed form LCMV BF from Springer (eq 47.7)
% W_LCMV(k,l) = (Z(k,l)Z^H(k,l))^-1 A(k) / (A^H(k)(Z(k,l)Z^H(k,l))^-1 A(k))

% Import audio files
[s1,Fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency

% Shorten the signals
NSig = 2^17; % Truncate signals to be a y*0.5 multiple of the window length, where y is an integer.
start = 29500;
s1 = s1(start:NSig-1+start); s2 = s2(start:NSig-1+start);

% For looking at narrowband sources...
% s1 = cos(2*pi*1113*[1:NSig]'/Fs1);
% s2 = cos(2*pi*541*[1:NSig]'/Fs1);

% For checking timing and ATF
% s1 = zeros(NSig,1); 
s2 = zeros(NSig,1);
% s1(500:549) = ones(50,1);
% s1(550:599) = -ones(50,1);

% Normalize the three signals to have a maximum amplitude just less than 1
% s1 = myNormalize(s1); s2 = myNormalize(s2); 

%% STFT on s, use if creating observation signals in f domain
K = 2^8+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.
% s1Padded = [zeros((K-1)/2,1);s1;zeros((K-1)/2,1)];
% [S1,L] = stft(s1Padded,K);
% S2 = stft(s2,K);
% mySpectrogram(S1);

%% Recreate signal
% sHat = myOverlapAdd(S1);

% Reconstructed source signal vs original source
% figure; plot(s1,'.'); hold on; plot(sHat,'o'); legend('s1','sHat');

%% Springer BF
% Place sensors and sources, note that s3 is diffuse noise so has no
% location

zPos = [3,1,1 ; 4,1,1 ]'; 
% zPos = [3,1,1 ; 3.6,1,1 ; 4.2,1,1 ; 4.8,1,1 ; 5.4,1,1 ; 6,1,1 ; 6.6,1,1 ; 7.2,1,1]'; 
% zPos = [2,1,1 ; 2.5,1,1 ; 3,1,1 ; 3.5,1,1 ; 3.1,2,1 ; 4.1,2,1 ; 5.3,2,1 ; 6.4,2,1]';
sPos = [3.5,3.5,1 ; 8,2,1]';
NSensors = length(zPos(1,:));
NSources = length(sPos(1,:));
figure; plot(zPos(1,:),zPos(2,:),'o'); hold on; plot(sPos(1,:),sPos(2,:),'o'); xlim([0 10]); ylim([0 5]); legend('Sensors','Sources');

% Construct observation signals using interpolation
Fs = Fs1; % Fs sets the sample rate of the observations
s = [s1,s2];
nsWt = 0; % nsWt = noise weight
[z, d] = myObservInterp(zPos,sPos,s,Fs1,Fs,nsWt); % Fs1 is the original sample rate, Fs is the desired sample rate. z is the observed signals, A is the acoustic transfer function (ATF)

% Calculate A, the acoustic transfer function (ATF)
c = 343; % c = velocity of sound in air (m.s^-1)
for k = 1:K%K-1
    for m = 1:NSensors
%         A(k,m) = exp(-j*2*pi*((k-1)*Fs/(K-1))*d(1,m)/c)/d(1,m);
        A(k,m) = exp(-j*2*pi*((k-1)*Fs/(K-1))*0.5)/d(1,m);
    end
end
figure; plot(real(fftshift(A(:,1)))); legend('A(:,1)');

ksym = [0:128,-128:-1];
for k = 1:K%K-1%[0:127,-127:-1]%1:K-1
    for m = 1:NSensors
%         Asym(k,m) = exp(-j*2*pi*((ksym(k))*Fs/(K-1))*d(1,m)/c)/d(1,m);
        Asym(k,m) = exp(-j*2*pi*((ksym(k))*Fs/(K))*0.5)/d(1,m);
    end
end
figure; plot(real(fftshift(Asym(:,1))));legend('Asym(:,1)');

% Use A to create the observations
sPadded = [zeros((K-1)/2,NSensors) ; s ; zeros((K-1)/2,NSensors)];
[S,LS] = stft(sPadded,K);
for l = 1:LS
    for m = 1:NSensors
        ZA(:,l,m) = A(1:end-1,m).*S(:,l);
    end
end
za = myNormalize(myOverlapAdd(ZA(:,:,1)));
za = za(65:end);

zn = myNormalize(z(:,1)); s1n = myNormalize(s1); 
figure; %plot(zn(358:end),'.'); 
hold on; plot(s1n); plot(real(za)); legend('s1n','za');
%plot(z(:,3));
% plot(z(:,4));plot(z(:,5));
% plot(z(:,6));plot(z(:,7));
% plot(z(:,8)); legend('1','2','3','4','5','6','7','8');

%% STFT on z, use if observation signals already created in t domain
zPadded = [zeros((K-1)/2,NSensors) ; z ; zeros((K-1)/2,NSensors)];
[Z,L] = stft(zPadded,K);

% figure; plot(abs(Z(:,34,1)),'.'); hold on; plot(abs(ZA(:,35,1)),'o'); legend('Z','ZA');

%% Delay and Sum
% What happens if I just find the weights that undo the delay? i.e. W^H*A =
% I
WH = pinv(A);
Y = zeros(K-1,L);
for l = 1:L
    for k = 1:K-1
        Y(k,l) = WH(:,k)' * squeeze(Z(k,l,:));
    end
end
y = myNormalize(myOverlapAdd(Y));
% plot(real(y(358:end)),'--')
% y = y(590:end);
% figure; plot(real(y)); hold on; %plot(s1); legend('y','s1');
% cov(100*real(y((K-1)/2:end-(K-1)/2-1)),s1);
% figure; plot(abs(Y));
% mySpectrogram(Y);
% mySpectrogram(Z(:,:,1));
% figure; plot(s1);hold on;  plot(real(y)); legend('s1','y')

%% Closed form optimal solution (Springer 47.7) - not suitable for time
% varying environments. W = ((ZZ^H)^-1 A)/(A^H (ZZ^H)^-1 A) -> (k x l x m)
% for k = 1:K-1
%     for l = 1:L
%         Zm = squeeze(Z(k,l,:));
%         ZZH = Zm * Zm';
%         Ak = A(k,:)';
%         W(k,l) = (inv(ZZH)*Ak)/(Ak'*inv(ZZH)*Ak); % inv(ZZH) results in a
% %         singular matrix warning from matlab
%     end
% end

%% Adaptive algorithm (Springer table 47.1)
% for k = 1:K-1
%     P(k,:,:) = eye(8) - A(k,:)*A(k,:)'/(norm(A(k,:))^2);
%     F(k,:) = A(k,:)/(norm(A(k,:))^2);    
% end
% 
% W = zeros(K-1,L,NSensors);
% W(:,1,:) = F(k);
% Y = zeros(K-1,L);
% mu = 0.5;
% for l = 1:L
%     for k = 1:K-1
%         Y(k,l) = squeeze(W(k,l,:))'*squeeze(Z(k,l,:));
%         W(k,l+1,:) = squeeze(P(k,:,:))*(squeeze(W(k,l,:))-mu*(squeeze(Z(k,l,:))*Y(k,l)'))+F(k);
%         Ypow(k,l) = Y(k,l)*Y(k,l)';
%     end
%     
% end
% y = myOverlapAdd(Y);
% figure; plot(real(y)); hold on; 
% plot(s1); legend('y','s1');

% % Ypow2 = Y*Y';
% for l = 1:L
%     Ypowt(l) = norm(Ypow(:,l));
% end
% figure; plot(Ypowt); 

%% Springer GSC
% Fixed bf


