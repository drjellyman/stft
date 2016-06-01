close all; clear all; 

% Import audio files
[s1,fs] = audioread('50_male_speech_english_ch10_orth_2Y.flac'); % s1 = source, Fs1 = sampling frequency
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac'); % s2 = noise stationary, Fs2 = sampling frequency
% Make them longer for longer adaptation
s1 = [s1;s1;s1;s1;s1];
s2 = [s2;s2;s2;s2;s2];

% Option to make s2 noise only
% lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 0.45,...
%     'StopbandFrequency', 0.55, 'PassbandRipple', 1, ...
%     'StopbandAttenuation', 60, 'DesignMethod', 'equiripple');
% s2 = filter(lpFilt,0.5*max(s1)*  randn(length(s2),1))  ;
s2 = 0.5*max(s1)*  randn(length(s2),1)  ;

% Downsample
ds_by = 3;
s1 = resample(s1,1,ds_by);
s2 = resample(s2,1,ds_by);
fs = fs/ds_by;

% s1 = zeros(length(s1),1);
% s2 = zeros(length(s2),1);

% Shorten the signals
K = 2^8+1; % Window length in samples. +1 makes the window length odd and symmetric. stft() will either truncate the window by 1, or overlap the 1st and last elements (which are both zero). This means we could consider the window length to be K-1 also depending on whether we would prefer to consider a symmetric or periodic window.
len_s = min(length(s1), length(s2));
len_s = len_s - mod(len_s,K-1);
s1 = s1(1:len_s);
s2 = s2(1:len_s);
s = [s1,s2]; 

%% STFT 
% pad the source signals so the 1st half window doesn't distort the data
s1Padded = [zeros((K-1)/2,1);s(:,1);zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s(:,2);zeros((K-1)/2,1)];

% Take stft of both sources
[S1,L] = stft(s1Padded,K);
[S2,L2] = stft(s2Padded,K);

%% Place sensors
NSources = length(s(1,:));
M = 8; % M = number of sensors
dz = 0.3; % ds = distance between sensors (m)
zPos = ones(3,M);
zPos(1,:) = zPos(1,:).*([0:M-1]*dz+1); % Set sensor position
sAng = [pi/3, 0]'; % Set source angle of arrival
c = 343; % Speed of sound (m/s)
dt = dz*sin(sAng)/c; % dt = time delay between sensors (s)

%% Create atf for both sources
fdom = (fs/(K-1)) * [0:(K-1)/2 , -(K-1)/2:-1]';
for m = 1:M
        A(:,m) = exp(-j*2*pi*fdom'*m*dt(1)) ;
        A2(:,m) = exp(-j*2*pi*fdom'*m*dt(2));
end

%% Create observations Z
Z = zeros(K,L,M);
for l = 1:L
    for m = 1:M
        Z(:,l,m) = A(:,m).*S1(:,l)+A2(:,m).*S2(:,l) ;
    end
end
z1 = myOverlapAdd(Z(:,:,1));

%% Remove redundant frequency info

%% A1) Frequency domain Frost optimum weights

for k=1:K
    R = zeros(M,M); % R is the spatial covariance of the inputs
    for l=1:L
        Ztmp = squeeze(Z(k,l,:));
        R = R + Ztmp*Ztmp'; % Sum the covariance over l
    end
    R = R + 1e-9*eye(M); % Diagonal loading
    Rinv = inv(R);
    Ak = A(k,:).';
    W(k,:) = Rinv*Ak/(Ak'*Rinv*Ak); % Calculate optimum weights vector 
end
Wopt = W;
 
% % Find output using optimum weights
% Y = zeros(K,L);
% for l=1:L
%     Ztmp = squeeze(Z(:,l,:));
%     Y(:,l) = sum(conj(W).*Ztmp,2).';
% end
% y = myOverlapAdd(Y);
% figure; plot(z1); hold on; plot(y); legend('z1','y optimum Frost');

%% A2) Adaptive Frost

% for k = 1:K
%     Atmp = A(k,:).';
%     P(k,:,:) = eye(M) - (Atmp*Atmp')/(norm(Atmp)^2);
%     F(k,:) = Atmp/(norm(Atmp)^2);
% end
% W = F; % Initialize weight vector
% 
% % Iterate
% mu = 0.0001; % mu = step size
% for l=1:L
%     Ztmp = squeeze(Z(:,l,:));
%     Y(:,l) = sum(W.*Ztmp,2); 
%     for k = 1:K
%         Ptmp = squeeze(P(k,:,:));
%         Ztmp = squeeze(Z(k,l,:));
%         Ftmp = F(k,:).';
%         Wtmp = W(k,:).';
%         W(k,:) = Ptmp*(Wtmp-mu*Ztmp*Y(k,l)')+Ftmp; % Update weights
%     end
%     Wdiff(l) = norm(Wopt-W);
% end
% y = myOverlapAdd(Y);
% figure; plot(y);

%% A2.b) Adaptive Frost using actual covariance
% mu = 0.0001; % mu = step size
% R = zeros(K,M,M);
% for k=1:K
%     for l=1:L
%         Ztmp = squeeze(Z(k,l,:));
%         R(k,:,:) = squeeze(R(k,:,:)) + Ztmp*Ztmp'; % Sum the covariance over l
%     end    
% end
% % Iteration over time
% Wdiff = zeros(l,1);
% for l=1:L
%     Ztmp = squeeze(Z(:,l,:));
%     Y(:,l) = sum(W.*Ztmp,2); 
%     for k = 1:K
%         Ptmp = squeeze(P(k,:,:));
%         Ztmp = squeeze(Z(k,l,:));
%         Ftmp = F(k,:).';
%         Wtmp = W(k,:).';
%         Rtmp = squeeze(R(k,:,:));
%         W(k,:) = Ptmp*(Wtmp-mu*Rtmp*Wtmp)+Ftmp; % Update weights
%     end
%     Wdiff(l) = norm(Wopt-W);
% end
% y = myOverlapAdd(Y);
% figure; plot(y);

%% B) GSC freq domain

% % Zal is with the look direction re-aligned
for k=1:K
    for l=1:L
        Zal(k,l,:) = (A(k,:)/(norm(A(k,:))^2))'.*squeeze(Z(k,l,:));
    end    
end
Yfbf = sum(Zal,3);
yfbf = myOverlapAdd(Yfbf); % For listening

% Create blocking matrix
B = zeros(M-1,M);
for m = 1:M-1
    B(m,m:m+1) = [1,-1];
end

% U is post blocking matrix
for k=1:K
    for l=1:L
        U(k,l,:) = B*squeeze(Zal(k,l,:));
    end
end
u1 = myOverlapAdd(U(:,:,1)); % For listening

% % Multi-channel Weiner
% for l=1:L
%     for k=1:K
%         Utmp = squeeze(U(k,l,:)); %(7x1)
%         Ucov = Utmp*Utmp'+1e-6*eye(M-1); % Covariance with diagonal loading
%         G(k,l,:) = inv(Ucov)*Utmp*Yfbf(k,l)';
%     end
%     Y(:,l) = Yfbf(:,l) - sum(squeeze(G(:,l,:)).*squeeze(U(:,l,:)) , 2);
% end
% y = myOverlapAdd(Y);
% figure; plot(y);

% % % % Optimal filter
% % % G = zeros(K,L,M-1);
% % % UUH = zeros(M-1,M-1);
% % % for k=1:K
% % %     for l=1:L
% % %         UUH = UUH + squeeze(U(k,l,:))*squeeze(U(k,l,:))';        
% % %     end
% % %     UUH = UUH + 1e-9*eye(M-1);
% % %     UUHinv(k,:,:) = inv(UUH); 
% % % end
% % % 
% % % for l=1:L
% % %     for k=1:K
% % %         G(:,l,:) = squeeze(UUHinv(k,:,:))*squeeze(U(k,l,:))*Yfbf(k,l)';        
% % %     end
% % % end
% % % %   G(k,:) = inv(UUH)*squeeze(U(k,l,:)).'*Yfbf(k,l);
  
% NLMS (eqn 23)
mu = 0.001; 
G = ones(K,M-1);
rho = 0.95;
Y = zeros(K,L);
Pest = zeros(K+1,1);
for l=1:L 
    Y(:,l) = Yfbf(:,l) - sum(G.*squeeze(U(:,l,:)) , 2);    
    for k=1:K
        Ztmp = squeeze(Z(k,l,:));
        Pest(k+1) = rho*Pest(k) +(1-rho)*(Ztmp'*Ztmp);
        G(k,:) = G(k,:) + mu*squeeze(U(k,l,:)).'*Y(k,l)'/Pest(k+1);        
    end  
    Pest2(l) = Pest'*Pest;
end
y = myOverlapAdd(Y);

%% Check final weight vector
st = linspace(-pi/2,pi/2,L); % Set up DOA space
dt = dz*sin(st)/c; % Set up time delay space
figure; hold on;
for k=2:127
    delays = exp(j*2*pi*fdom(k)*[1:M]'*dt);
%     P = 20*log10(abs(W(k,:)*delays).^2);
    P = 20*log10(abs(G(k,:)*delays(1:7,:)).^2);
    plot(st,P);
end
plot([sAng(1), sAng(1)],[-200, 50],'--', [sAng(2), sAng(2)],[-200, 50],'--');
grid on;


