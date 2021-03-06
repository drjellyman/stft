% BF conjugate symmetric
close all; clear all;

% Import audio files
[s1,fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % ns = noise stationary, Fs2 = sampling frequency

% s1 = resample(s1,1,3);
% s2 = resample(s2,1,3);
% fs = fs/3;

% s1 = zeros(length(s1),1);
% s2 = zeros(length(s2),1);

% Shorten the signals
NSig = 2^17; % Truncate signals to be a y*0.5 multiple of the window length, where y is an integer.
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
clear L2;

% Discard top half of the spectrum, dc, and fs/2.
% S1 = S1(2:(K-1)/2,:);
% S2 = S2(2:(K-1)/2,:);
% K1 = length(S1(:,1)); % K1 = the truncated number of frequency bins

%% Place sensors
NSources = length(s(1,:));
M = 8; % M = number of sensors
dz = 0.9; % ds = distance between sensors
zPos = ones(3,M);
zPos(1,:) = zPos(1,:).*([0:M-1]*dz+1);
sPos = [1,2,1 ; 8,2.5,1]';
figure; plot3(zPos(1,:),zPos(2,:),zPos(3,:),'o',sPos(1,:),sPos(2,:),sPos(3,:),'*'); grid on; 
xlabel('x'); ylabel('y'); zlabel('z'); xlim([0, 9]); ylim([0,5]); zlim([0,3]); 

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


% for ns = 1:NSources
%     a = sPos(:,ns) - zmean';
%     b = zPos(:,ns) - zmean';
%     phi(ns) = atan2(norm(cross(a,b)),dot(a,b));
%     if phi(ns) > pi/2
%         phi(ns) = phi(ns)-pi/2;
%     end
%     dzphi(:,ns) = [0:M-1]'*dz*cos(phi(ns))+d(ns);
% end
% dzphi(:,2) = flipud(dzphi(:,2));



% for ns = 1:NSources 
%     dib(ns) = atan(norm(cross(sPos(:,ns)-zmean',sPos(:,ns)-zPos(:,i(ns))))); % di = dinitial for planar wave so trades between first and last sensors, i.e. 1st sensor is delayed less, 
% end
% for ns = 1:NSources
%     dzi(ns) = norm(sPos(:,ns)-zPos(:,i(ns)))*cos(dib(ns));
% end
% 
% for ns = 1:NSources
%     a = zmean - sPos(:,ns)';
%     b = zmean - (zmean+[0,2,0]);
%     phi(ns) = atan2(norm(cross(a,b)),dot(a,b));
%     dzphi(:,ns) = [0:M-1]'*dz*cos(phi(ns))+dzi(ns);
% end



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
        Z(:,l,m) = A(:,m).*S1(:,l)+A2(:,m).*S2(:,l);
    end
end

%% ifft of Z1 (sensor 1)
% Z1 = squeeze(Z(:,:,1));
% Z8 = squeeze(Z(:,:,8));
% Z1 = [zeros(1,L) ; Z1 ; zeros(1,L) ; (flipud(Z1)) ; zeros(1,L)];
% mySpectrogram(Z1)
% z1 = myOverlapAdd(Z1);
% z8 = myOverlapAdd(Z8);
% figure; plot(myNormalize(z1)); hold on; plot(myNormalize(z8(24:end))); legend('z1','z8');

%% Filter and sum fixed bf
% W0 = A/||A||^2
% Springer version
% for m = 1:M
%     W0(:,m) =  A(:,m)/(norm(A(:,m))^2); %pinv(A(:,m));
% end
% 
% % Yfbf = W0^H * Z
% Yfbf = zeros(K,L,M);
% for l = 1:L
%     for m = 1:M
%         Yfbf(:,l,m) = squeeze(conj(W0(:,m))).*squeeze(Z(:,l,m));
%     end
% end
% Yfbf = sum(Yfbf,3);

% Dayle version
for k=1:K
    for l=1:L
        ZZ(k,l,:) = A(k,:)'.*squeeze(Z(k,l,:));
    end    
end
Yfbf = sum(ZZ,3);

% figure; hold on;
% for m=1:M
%     zz(:,m) = myOverlapAdd(ZZ(:,:,m));
% end

% blocking matrix
B = zeros(M-1,M);
for m = 1:M-1
    B(m,m:m+1) = [1,-1];
end
for k=1:K
    for l=1:L
        ZZZ(k,l,:) = B*squeeze(ZZ(k,l,:));
    end
end
W = ones(K,M-1);
Y = zeros(K,L);
mu = 0.8;
for l=1:L
    ZZZZ = W(:,:).*squeeze(ZZZ(:,l,:));
    muNorm = mu/(squeeze(sum(ZZZ(:,l,:),3))'*squeeze(sum(ZZZ(:,l,:),3)))
    Yanc = sum(ZZZZ,2);
    Y(:,l) = Yfbf(:,l)-Yanc;
    for k=1:K
        W(k,:) = W(k,:)+muNorm*Y(k,l)*squeeze(ZZZ(k,l,:)).';%MLS 
    end
end
y = myOverlapAdd(Y);
figure; plot(y);

% What can I compare it to? sum(ZZ)?
ZZsum = sum(ZZ,3);
zzsum = myOverlapAdd(ZZsum);

%% ifft of Yfbf 
% yfbf = myOverlapAdd(Yfbf);
% % figure; plot(myNormalize(abs(fft(yfbf)))); hold on;
% % plot(abs(myNormalize(fft(s(:,2)))));
% % 
% % Zsum = sum(Z,3);
% % zsum = myOverlapAdd(Zsum);
% % 
% % audiowrite('yfbf.flac',myNormalize(yfbf),fs);
% % audiowrite('s1.flac',myNormalize(s(:,1)),fs);
% % audiowrite('s2.flac',myNormalize(s(:,2)),fs);
% % audiowrite('zsum.flac',myNormalize(zsum),fs);
% 
% % xc = xcorr(z1,z8);
% % figure; plot(xc);
% % [mxc, mxci] = max(xc)
% 
%% Adaptive noise canceller
% % Hones = repmat([1,-1],K,M/2);
% % 
% % H = Hones.*squeeze(W0(:,:));
% H = zeros(M-1,M);
% 
% 
% G = ones(K,L,M);
% % Find Yanc = GHZ
% Y = zeros(K,L,M);
% U = zeros(K,L,M);
% alpha = 0.5;
% mu = 0.5;
% Pest = zeros(K,L);
% Yanc = zeros(K,L);
% for l = 1:L
%     for k = 1:K
%         U(k,l,:) = squeeze(H(k,:))'.*squeeze(Z(k,l,:));
%         Yanc(k,l) = sum(conj(squeeze(G(k,l,:))).*squeeze(U(k,l,:)));%conj(squeeze(H(k,:)))*squeeze(Z(k,l,:));
%         Y(k,l) = Yfbf(k,l) - Yanc(k,l);
%         Pest(k,l+1) = alpha*Pest(k,l)+(1-alpha)*dot(Z(k,l,:),Z(k,l,:));%(Z(k,l,:)'*Z(k,l,:))
%         G(k,l+1,:) = G(k,l,:)+mu*((U(k,l,:)*Y(k,l)')/(Pest(k,l)));
%         
% %         for m = 1:M
% % %         Yanc(k,l) = squeeze(G(k,l,:))' *conj(squeeze(H(k,:)))*squeeze(Z(k,l,:));
% %             %U(k,l,m) = H(k,m)'*Z(k,l,m);
% %             Gtilde(k,l+1,m) = G(k,l,m)+mu*((U(k,l,m))/());
% %         end
%     end
% %     Gtilde(k,l+1,m) = G(k,l,m)+mu*((H())/());
% end
% 
% % Yanc = sum(Yanc,3);
% % yanc = myOverlapAdd(Yanc);
% 
% % %% Adaptive filtering of Yanc
% % G = ones(K,L,M);
% % for l = 1:L
% %     Yanc = G .* Yanc;
% % end
% 
% % 
% % Y = Yfbf-Yanc;
% y = myOverlapAdd(Y);
% yanc = myOverlapAdd(Yanc);
% % figure; plot(real(yanc));legend('yanc');
% % figure; plot(real(yfbf));legend('yfbf');
% % figure; plot(real(y));legend('real(y)');
% figure; plot(real(yanc(1:128)));hold on; plot(real(yfbf(1:128)));
% plot(real(y(1:128)));legend('yanc','yfbf','y');
% 
% 
