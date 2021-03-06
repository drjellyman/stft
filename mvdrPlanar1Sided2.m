% BF conjugate symmetric
close all; clear all;

% Import audio files
[s1,fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency

s1 = resample(s1,1,3);
s2 = resample(s2,1,3);
fs = fs/3;

% s1 = zeros(length(s1),1);
s2 = zeros(length(s2),1);

% Shorten the signals
NSig = 2^16; % Truncate signals to be a y*0.5 multiple of the window length, where y is an integer.
start = 1;
s1 = s1(start:NSig-1+start); s2 = s2(start:NSig-1+start);
s = [s1,s2]; 

%% STFT 
K = 2^8+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.

% pad the source signals so the 1st half window doesn't distort the data
s1Padded = [zeros((K-1)/2,1);s(:,1);zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s(:,2);zeros((K-1)/2,1)];

% Take stft of both sources
[S1,L] = stft(s1Padded,K);
[S2,L2] = stft(s2Padded,K);

%% Place sensors
NSources = length(s(1,:));
M = 8; % M = number of sensors
dz = 0.06; % ds = distance between sensors
zPos = ones(3,M);
zPos(1,:) = zPos(1,:).*([0:M-1]*dz+1);
sPos = [1,2,1 ; 8,2.5,1]';
% figure; plot3(zPos(1,:),zPos(2,:),zPos(3,:),'o',sPos(1,:),sPos(2,:),sPos(3,:),'*'); grid on; 
% xlabel('x'); ylabel('y'); zlabel('z'); xlim([0, 9]); ylim([0,5]); zlim([0,3]); 

% Find the minimum distance from each source to the array
for m = 1:M
    for ns = 1:NSources
        dmin(m,ns) = norm(sPos(:,ns)-zPos(:,m));
    end
end
[dmin,i] = min(dmin);     % d contains the minimum distance, i contains the index of the closest sensor
zmean = mean(zPos')';

% find the delay (in meters) between sensors for each source
% theta = acos(u'*v/(norm(u)*norm(v)));
u = sPos-repmat(zmean,1,NSources);
v = (zmean+[0,1,0]')-zmean;
for ns = 1:NSources
    theta(ns) = acos((u(:,ns)'*v)/(norm(u(:,ns))*norm(v)));
    d(ns) = dz*sin(theta(ns));
end

% Find the initial delay (in meters) for both sensors
for ns = 1:NSources
    di(ns) = norm(u(:,ns))-((M-1)/2)*d(ns);
end

% Make delay vectors for ATF
D = [[0:M-1]'*d(1)+di(1),[0:M-1]'*d(2)+di(2)];

% Flip delay vectors for sources that hit sensor 8 first
for ns = 1:NSources
    if i(ns) == M
        D(:,ns) = flipud(D(:,ns));
    end
end

%% Create atf for both sources
c = 343; % Speed of sound in m.s^-1
% kdom = [1:K]'*fs/(K-1);
% kdom = [0:K-1]'*fs/(K-1);
kdom = (fs/(K-1)) * [0:(K-1)/2 , -(K-1)/2:-1]';
for m = 1:M
    A(:,m) = exp(-j*2*pi*kdom'*D(m,1)/c) ;%/ D(m,1);
    A2(:,m) = exp(-j*2*pi*kdom'*D(m,2)/c);%/ D(m,2);
end

%% Create observations
for l = 1:L
    for m = 1:M
        Z(:,l,m) = A(:,m).*S1(:,l)+A2(:,m).*S2(:,l); %A(:,m).*
    end
end

% check time domain 
for m=1:M
    z(:,m) = myOverlapAdd(Z(:,:,m));
end
figure; hold on; 
plot(z(:,1)); plot(z(:,8))
legend('1','8');
xc12 = xcorr(z(:,1),z(:,8));
figure; plot(xc12)
[mx, mxi] = max(xc12)
mxi-length(z(:,1))

%% GSC - Griffiths + Jim freq domain
% Dayle version
% ZZ is wiht the look direction re-aligned
% for k=1:K
%     for l=1:L
%         ZZ(k,l,:) = A(k,:)'.*squeeze(Z(k,l,:));
%     end    
% end
% Yfbf = sum(ZZ,3);
% 
% % blocking matrix
% B = zeros(M-1,M);
% for m = 1:M-1
%     B(m,m:m+1) = [1,-1];
% end
% % ZZZ is post blocking matrix
% for k=1:K
%     for l=1:L
%         ZZZ(k,l,:) = B*squeeze(ZZ(k,l,:));
%     end
% end
% % ZZZZ is post blocking matrix and post adaptive weights
% W = ones(K,M-1);
% Y = zeros(K,L);
% mu = 0.21;
% for l=1:L
%     ZZZZ = W(:,:).*squeeze(ZZZ(:,l,:));
% %     muNorm = mu/(squeeze(sum(ZZZ(:,l,:),3))'*squeeze(sum(ZZZ(:,l,:),3)));
%     Yanc = sum(ZZZZ,2);
%     Y(:,l) = Yfbf(:,l)-Yanc;
%     for k=1:K
%         W(k,:) = W(k,:)+mu*Y(k,l)*squeeze(ZZZ(k,l,:)).';%MLS % muNorm
%     end
% end
% y = myOverlapAdd(Y);
% figure; plot(y);
% audiowrite('y.flac',myNormalize(y),fs);
% % 
% % % What can I compare it to? sum(ZZ)?
% ZZsum = sum(ZZ,3);
% zzsum = myOverlapAdd(ZZsum);
% audiowrite('zzsum.flac',myNormalize(zzsum),fs);

%% MVDR using analytical solution i.e. optimal weights
% W = (R^-1*A) * (A*R^-1*A)^-1 % From Lorenz and Boyd, Robust Minimum Variance
% Beamforming, but they are referring to the Capon beamformer.
% W = ones(M,1);
% for l=1:L
%     for k=1:K
% %         Y(k,l) = sum(conj(W).*squeeze(Z(k,l,:)));
%         Rz = squeeze(Z(k,l,:))*squeeze(Z(k,l,:))';
%         Rz = Rz + 0.0000000001*eye(M);
%         Rzi = inv(Rz);
%         AA = squeeze(A(k,:)).';
%         W = Rzi*AA/(AA'*Rzi*AA);
%         Y(k,l) = sum(conj(W).*squeeze(Z(k,l,:)));
%     end
% end
% y=myOverlapAdd(Y);
% audiowrite('y.flac',myNormalize(y),fs);
% figure; plot(real(y));
% sound(myNormalize(real(y)),fs);

%% MVDR using adaptive LMS
% W = ones(M,1);
% Y = zeros(K,L);
% for l=1:L
%     for k=1:K
%         Rz = squeeze(Z(k,l,:))*squeeze(Z(k,l,:))';
%         Rz = Rz + 0.00000001*eye(M);
%         Rzi = inv(Rz);
%         AA = squeeze(A(k,:)).';
%         
%         W = Rzi*AA/(AA'*Rzi*AA);      
%         Y(k,l) = sum(W.*squeeze(Z(k,l,:)));
% %         W = W + mu*Y(k,l)*squeeze(Z(k,l,:));
%     end
% end
% y=myOverlapAdd(Y);
% % audiowrite('y.flac',myNormalize(y),fs);
% figure; plot(real(y));


%% Freq dom frost
% P = eye()-AAH/(norm(A));
% for k=1:K
%     P(k,:,:) = eye(M)-A(k,:)'*A(k,:)/(norm(A(k))^2);
%     F(k,:) = A(k,:)/(norm(A(k))^2);
% end
% W = F(1,:).';
% mu = 0.1;
% Y = zeros(K,L);
% for l=1:L
%     for k=1:K
%         Zt = squeeze(Z(k,l,:));
%         Y(k,l) = W'*Zt;
%         Pt = squeeze(P(k,:,:));        
%         Ft = squeeze(F(k,:)).';
%         W = Pt*(W-mu*Zt*Y(k,l)')+Ft;
%     end
% end
% y=myOverlapAdd(Y);
% figure; plot(real(y));
% sound(myNormalize(real(y)),fs);