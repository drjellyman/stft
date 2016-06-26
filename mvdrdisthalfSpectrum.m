close all; clear all; 

% Import audio files
% [s1,fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s1 = source, Fs1 = sampling frequency
% [s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % s2 = noise stationary, Fs2 = sampling frequency
[s1,fs1] = audioread('317354__speedenza__shall-i-compare-thee-voice.wav'); % s1 = source, Fs1 = sampling frequency
[s2,fs2] = audioread('273177__xserra__la-vaca-cega-eva.wav'); % s2 = noise stationary, Fs2 = sampling frequency
% Make them longer for longer adaptation
s1 = [s1];
s2 = [s2(:,1)];

% Make interferer noise
% s2 = 0.5*max(s1)*  randn(length(s2),1)  ;

% Downsample
fsd = 16e3; % fsd = desired sampling frequency
fs = fs1/round(fs1/fsd); % fs = actual sampling frequency post resample
s1 = resample(s1,1,round(fs1/fsd));
s2 = resample(s2,1,round(fs2/fsd));

% s1 = zeros(length(s1),1);
% s2 = zeros(length(s2),1);

% Shorten the signals
K = 2^9+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.
lmin = min([length(s1),length(s2)]);
ld = 30*fs-mod(30*fs,K-1); % ld = desired length = seconds * seconds^-1
if lmin < ld
    s1 = s1(1:lmin); s2 = s2(1:lmin); 
else
    s1 = s1(1:ld); s2 = s2(1:ld); 
end
s = [s1,s2];
% Normalize input files to have rms = 1
srms = rms(s);
srmsinv = 1./srms;
s = s * (0.1*diag(srmsinv));





% len_s = min(length(s1), length(s2));
% len_s = len_s - mod(len_s,K-1);
% s1 = s1(1:len_s);
% s2 = s2(1:len_s);
% s = [s1,s2];    
% srms = rms(s);
% s = [s(:,1).*0.01/srms(1), s(:,2).*0.01/srms(2)];

%% STFT 
% pad the source signals so the 1st half window doesn't distort the data
s1Padded = [zeros((K-1)/2,1);s(:,1);zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s(:,2);zeros((K-1)/2,1)];
% s1savePadded = [zeros((K-1)/2,1);s1save;zeros((K-1)/2,1)];

% Take stft of both sources and truncate to exclude negative frequencies
% as well as dc and fs/2.
[S1,L] = stft(s1Padded,K);
S1half = S1(2:(K+1)/2-1,:);
[S2,L2] = stft(s2Padded,K);
S2half = S2(2:(K+1)/2-1,:);

%% Place sensors
NSources = length(s(1,:));
M = 30; % M = number of sensors
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

%% Create atf for both sources (planar)
% fdom = (fs/(K-1)) * [0:(K-1)/2 , -(K-1)/2:-1]';
Khalf = (K-1)/2-1;
fdom = (fs/(K-1)) * [1:Khalf]';
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

Z = zeros(Khalf,L,M);
for l = 1:L
    for m = 1:M
        Z(:,l,m) = A(:,m).*S1half(:,l)+A2(:,m).*S2half(:,l);
    end
end
% z1 = myOverlapAdd(Z(:,:,1));
% subplot(3,1,3); plot(z1); set(gca,'fontsize',14); legend('Observation Mic #1');
% axis tight;ylim([-0.06, 0.06]);

%% A1) Frequency domain Frost optimum weights

for k=1:Khalf
    R = zeros(M,M); % R is the spatial covariance of the inputs
    for l=1:L
        Ztmp = squeeze(Z(k,l,:));
        R = R + Ztmp*Ztmp'; % Sum the covariance over l
    end
    R = R + 1e-6*eye(M); % Diagonal loading
    Rinv = inv(R);
    Ak = A(k,:)';
    W(k,:) = Rinv*Ak/(Ak'*Rinv*Ak); % Calculate optimum weights vector 
end
Wopt = W;
 
% Find output using optimum weights
Yopt = zeros(Khalf,L);
for l=1:L
    Ztmp = squeeze(Z(:,l,:));
    Yopt(:,l) = sum((W).*Ztmp,2).';
end
Yopt = [zeros(1,L);Yopt;zeros(2,L);conj(flipud(Yopt))];
yopt = myOverlapAdd(Yopt);

%% A2) Adaptive Frost
% 
for k = 1:Khalf
    Atmp = A(k,:).';
    P(k,:,:) = eye(M) - (Atmp*Atmp')/(norm(Atmp)^2);
    F(k,:) = Atmp/(norm(Atmp)^2);
end
W = F; % Initialize weight vector

% Iterate
mu = 1e-5 ; % mu = step size
for l=1:L
    Ztmp = squeeze(Z(:,l,:));
    Y(:,l) = sum(conj(W).*Ztmp,2); 
    for k = 1:Khalf
        Ptmp = squeeze(P(k,:,:));
        Ztmp = squeeze(Z(k,l,:));
        Ftmp = F(k,:).';
        Wtmp = W(k,:).'; %.'
        Rtmp = Ztmp*Ztmp';
%         W(k,:) = Ptmp*(Wtmp-mu*Ztmp*Y(k,l)')+Ftmp; % Update weights
        W(k,:) = Ptmp*((Wtmp)-mu*Rtmp*(Wtmp))+Ftmp; % Update weights
    end
%     SINR = WDDW/WQW;
%     SINR = 
end

% Create two sided Y for ifft
Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];

% Take ifft
y = myOverlapAdd(Y);
figure; plot(y);

%% Check final weight vector
st = linspace(-pi/2,pi/2,L); % Set up DOA space
dt = dz*sin(st)/c; % Set up time delay space
figure; hold on;
plot([sAng(1), sAng(1)],[-200, 50],'--', [sAng(2), sAng(2)],[-200, 50],'--');
% Psave = zeros(Khalf,1);
for k=1:Khalf
    delays = exp(j*2*pi*fdom(k)*[1:M]'*dt);
    P = 20*log10(abs(W(k,:)*delays).^2);
    plot(st,P);
%     Psave(k)=P(1745);
end

figure; hold on;
plot([sAng(1), sAng(1)],[-200, 50],'--', [sAng(2), sAng(2)],[-200, 50],'--');
for k=1:Khalf
    delays = exp(j*2*pi*fdom(k)*[1:M]'*dt);
    P = 20*log10(abs(Wopt(k,:)*delays).^2);
    plot(st,P);
%     Psave(k)=P(1745);
end
grid on; xlabel('Look direction (radians)');ylabel('Power (dB)');
legend('Target','Interferer'); axis tight; set(gca,'fontsize',14);

% figure; plot(Psave,'*'); legend('Power at interferer direction by frequency bin');

%% Use final weight vector to process mixed z's
Zb = zeros(Khalf,L,M);

for l = 1:L
    for m = 1:M
%         Zb(:,l,m) = A(:,m).*S1half(:,l)+A2(:,m).*S2half(:,l) ; % S1save
        Zbtarget(:,l,m) = A(:,m).*S1half(:,l); % S1save
        Zbinterf(:,l,m) = A2(:,m).*S2half(:,l) ; % S1save
        Zb(:,l,m) = Zbtarget(:,l,m)+ Zbinterf(:,l,m);

    end
end
% z1b = myOverlapAdd(Z(:,:,1));

for l=1:L
    Zbttmp = squeeze(Zbtarget(:,l,:));
    Zbitmp = squeeze(Zbinterf(:,l,:));
    Zbtmp = squeeze(Zb(:,l,:));
    Ybt(:,l) = sum(conj(W).*Zbttmp,2); 
    Ybi(:,l) = sum(conj(W).*Zbitmp,2); 
    Yb(:,l) = sum(conj(W).*Zbtmp,2);
end
% Yb = [zeros(1,L) ; conj(flipud(Yb)) ; zeros(1,L) ; Yb; zeros(1,L)];
Ybt = [zeros(1,L);Ybt;zeros(2,L);conj(flipud(Ybt))];
Ybi = [zeros(1,L);Ybi;zeros(2,L);conj(flipud(Ybi))];
Yb = [zeros(1,L);Yb;zeros(2,L);conj(flipud(Yb))];

ybt = myOverlapAdd(Ybt);
ybi = myOverlapAdd(Ybi);
yb = myOverlapAdd(Yb);

% Find SIR 
SIR = 10*log10(ybt'*ybt / (ybi'*ybi));


