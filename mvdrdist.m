close all; clear all; 

% Import audio files
% [s1,fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s1 = source, Fs1 = sampling frequency
% [s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % s2 = noise stationary, Fs2 = sampling frequency
[s1,fs] = audioread('317354__speedenza__shall-i-compare-thee-voice.wav'); % s1 = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('273177__xserra__la-vaca-cega-eva.wav'); % s2 = noise stationary, Fs2 = sampling frequency
% Make them longer for longer adaptation
s1 = [s1];%;s1;s1;s1;s1;s1;s1;s1;s1;s1;s1;s1;s1;s1;s1;s1];
s2 = [s2(:,1)];%;s2;s2;s2;s2;s2;s2;s2;s2;s2;s2;s2;s2;s2;s2;s2];

% Option to make s2 noise only
% lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 0.45,...
%     'StopbandFrequency', 0.55, 'PassbandRipple', 1, ...
%     'StopbandAttenuation', 60, 'DesignMethod', 'equiripple');
% s2 = filter(lpFilt,0.5*max(s1)*  randn(length(s2),1))  ;
% s2 = 0.5*max(s1)*  randn(length(s2),1)  ;

% Downsample
ds_by = 3;
s1 = resample(s1,1,ds_by);
s2 = resample(s2,1,ds_by);
fs = fs/ds_by;

s1save = s1; 
% s1 = zeros(length(s1),1);
% s2 = zeros(length(s2),1);

% Shorten the signals
K = 2^9+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.
len_s = min(length(s1), length(s2));
len_s = len_s - mod(len_s,K-1);
s1 = s1(1:len_s);
s2 = s2(1:len_s);
s1save = s1save(1:len_s);
s = [s1,s2];    

% figure; subplot(3,1,1); plot(s1);axis tight; ylim([-0.06, 0.06]);set(gca,'fontsize',14); legend('Target'); subplot(3,1,2); plot(s2); 
% legend('Interferer'); grid on; set(gca,'fontsize',14); axis tight;
% ylim([-0.06, 0.06]);

%% STFT 
% pad the source signals so the 1st half window doesn't distort the data
s1Padded = [zeros((K-1)/2,1);s(:,1);zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s(:,2);zeros((K-1)/2,1)];
s1savePadded = [zeros((K-1)/2,1);s1save;zeros((K-1)/2,1)];

% Take stft of both sources
[S1,L] = stft(s1Padded,K);
[S2,L2] = stft(s2Padded,K);
[S1save,L3] = stft(s1savePadded,K);
% 
% figure; mySpectrogram(S1);

%% Place sensors
NSources = length(s(1,:));
M = 32; % M = number of sensors
dz = 0.01; % ds = distance between sensors (m)
zPos = ones(3,M);
zPos(1,:) = zPos(1,:).*([0:M-1]*dz+1); % Set sensor position
sAng = [pi/3, 0]'; % Set source angle of arrival
c = 343; % Speed of sound (m/s)
dt = dz*sin(sAng)/c; % dt = time delay between sensors (s)
dd = dz*sin(sAng(1));

figure; hold on;
d = 0.7;
sPos = mean(zPos,2) + [d*sin(sAng(1)), d*cos(sAng(1)),0]';
plot(sPos(1),sPos(2),'r*');
nPos = mean(zPos,2) + [d*sin(sAng(2)), d*cos(sAng(2)),0]';
plot(nPos(1),nPos(2),'r^');
for m=1:M
    plot(zPos(1,m),zPos(2,m),'bo');
end
xlim([0.5, 2.5]); ylim([0.5, 2]);
legend('Target','Interference','Sensors');
grid on; 
set(gca, 'fontsize', 14);

% for m=1:M
%     dist(:,m) = [norm(sPos(1)-zPos(:,m)), norm(sPos(2)-zPos(:,m))];
% end
% figure; plot(dist(1,:),'*:')
% grid on; set(gca,'fontsize',14); xlim([1,8]); xlabel('Sensor'); ylabel('Distance (m)');

%% Create atf for both sources (planar)
fdom = (fs/(K-1)) * [0:(K-1)/2 , -(K-1)/2:-1]';
for m = 1:M
        A(:,m) = exp(-j*2*pi*fdom'*m*dt(1)) ;
        A2(:,m) = exp(-j*2*pi*fdom'*m*dt(2));
end

%% Create atf for both sources (point)
% fdom = (fs/(K-1)) * [0:(K-1)/2 , -(K-1)/2:-1]';
% 
% for m = 1:M
%         A(:,m) = exp(-j*2*pi*fdom'*dist(1,m)/c) ;
%         A2(:,m) = exp(-j*2*pi*fdom'*dist(2,m)/c);
% end

%% Create observations Z
Z = zeros(K,L,M);
for l = 1:L
    for m = 1:M
        Z(:,l,m) = A(:,m).*S1(:,l)+A2(:,m).*S2(:,l) ;
    end
end
z1 = myOverlapAdd(Z(:,:,1));
% subplot(3,1,3); plot(z1); set(gca,'fontsize',14); legend('Observation Mic #1');
% axis tight;ylim([-0.06, 0.06]);

%% Remove redundant frequency info

%% A1) Frequency domain Frost optimum weights

% for k=1:K
%     R = zeros(M,M); % R is the spatial covariance of the inputs
%     for l=1:L
%         Ztmp = squeeze(Z(k,l,:));
%         R = R + Ztmp*Ztmp'; % Sum the covariance over l
%     end
%     R = R + 1e-9*eye(M); % Diagonal loading
%     Rinv = inv(R);
%     Ak = A(k,:).';
%     W(k,:) = Rinv*Ak/(Ak'*Rinv*Ak); % Calculate optimum weights vector 
% end
% Wopt = W;
 
% Find output using optimum weights
% Y = zeros(K,L);
% for l=1:L
%     Ztmp = squeeze(Z(:,l,:));
%     Y(:,l) = sum(conj(W).*Ztmp,2).';
% end
% y = myOverlapAdd(Y);
% y = y((K-1)/2+1:end-(K-1)/2);
% figure; plot(s1); hold on; plot(y); legend('Target signal','BF output');
% set(gca,'fontsize',14); grid on; axis tight; ylim([-0.05, 0.05]);


%% A2) Adaptive Frost
% 
for k = 1:K
    Atmp = A(k,:).';
    P(k,:,:) = eye(M) - (Atmp*Atmp')/(norm(Atmp)^2);
    F(k,:) = Atmp/(norm(Atmp)^2);
end
W = F; % Initialize weight vector

% Iterate
mu = 0.0001 ; % mu = step size
lambda = 0;%1e-12; % lambda = diagonal loading % 
% muMaybe = zeros(K,L);
for l=1:L
    Ztmp = squeeze(Z(:,l,:));
    Y(:,l) = sum(W.*Ztmp,2); 
    for k = 1:K
        Ptmp = squeeze(P(k,:,:));
        Ztmp = squeeze(Z(k,l,:));
%         lambda = 2*Ztmp'*Ztmp; % optimal lambda is 2 * noise power?
%         doesn't work
        Ftmp = F(k,:).';
        Wtmp = W(k,:).';
        Rtmp = Ztmp*Ztmp'+lambda*eye(M);
%         PowScale = 0.01*Ztmp'*Ztmp;
%         W(k,:) = Ptmp*(Wtmp-mu*Ztmp*Y(k,l)')+Ftmp; % Update weights
        W(k,:) = Ptmp*(Wtmp-mu*Rtmp*Wtmp)+Ftmp; % Update weights

    end
%     Wdiff(l) = norm(Wopt-W);
end
y = myOverlapAdd(Y);
figure; plot(y);

%% Check final weight vector
st = linspace(-pi/2,pi/2,L); % Set up DOA space
dt = dz*sin(st)/c; % Set up time delay space
figure; hold on;
plot([sAng(1), sAng(1)],[-200, 50],'--', [sAng(2), sAng(2)],[-200, 50],'--');
Psave = zeros(126,1);
for k=2:127
    delays = exp(j*2*pi*fdom(k)*[1:M]'*dt);
    P = 20*log10(abs(W(k,:)*delays).^2);
    plot(st,P);
    Psave(k)=P(1745);
end

grid on; xlabel('Look direction (radians)');ylabel('Power (dB)');
legend('Target','Interferer'); axis tight; set(gca,'fontsize',14);

figure; plot(Psave,'*'); legend('Power at interferer direction by frequency bin');

%% Use final weight vector to process mixed z's
Zb = zeros(K,L,M);

for l = 1:L
    for m = 1:M
        Zb(:,l,m) = A(:,m).*S1(:,l)+A2(:,m).*S2(:,l) ; % S1save
    end
end
z1b = myOverlapAdd(Z(:,:,1));


for l=1:L
    Zbtmp = squeeze(Zb(:,l,:));
    Yb(:,l) = sum(W.*Zbtmp,2); 
end
yb = myOverlapAdd(Yb);
figure; plot(yb);
