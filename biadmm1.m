%% Start again mvdr beamformer, hopefully distributed using biadmm. Needs 
% to have 50 sensors randomly placed within a 100x100x100 m free space. One
% target speech signal, and one noise signal, also randomly placed. Signal
% of interest is a 20 s speech signal chosen randomly from a 60 s file. fs
% = 16 kHz. Window length is 25 ms with a Hann window and 50% overlap.
% Interference is a randomly placed, zero meaqn gaussian point source with
% power equal to -5, 0, 5 dB when compared to the target signal. 
close all; clear;

% % Import target audio
% [s1,fs1] = audioread('273177__xserra__la-vaca-cega-eva.wav'); 
% 
% % Downsample
% fsd = 16e3; % fsd = desired sampling frequency
% fs = fs1/round(fs1/fsd); % fs = actual sampling frequency post resample
% s1 = resample(s1,1,round(fs1/fsd));
% 
% % Truncate to desired length, ensure that the length is a multiple of 
% % the window length, and randomly select a section of the audio file.
% K = 2^9+1; % K = window length in samples, and the number of frequency bins
% tls = 20; % tls = target length in seconds
% tl = tls*fs-mod(tls*fs,K-1); % tl = target length in samples, adjusted for window length and sampling frequency
% start = floor((length(s1)-tl)*rand);
% s1 = s1(start:start+tl-1,1); % Truncate s1 to be one channel, and 20 s long
% 
% % Normalize the target audio file to make it easy to change files
% s1rms = rms(s1);
% s1rmsinv = 1./s1rms;
% s1 = s1 * (0.1*diag(s1rmsinv));% This should probably be scaled down to avoid clipping
% 
% % Set up interferer with equal power, i.e. snr = 0 dB 
% s1Pow = (s1'*s1) / length(s1);
% s2 = sqrt(s1Pow) *  randn(length(s1),1); % Currently set up for equal power
% s2Pow = (s2'*s2) / length(s2);
% SourceSNRdB = 10*log10(s1Pow/s2Pow)
% 
% %% STFT 
% % pad the source signals so the 1st half window doesn't distort the data
% s1Padded = [zeros((K-1)/2,1);s1;zeros((K-1)/2,1)];
% s2Padded = [zeros((K-1)/2,1);s2;zeros((K-1)/2,1)];
% 
% % Take stft of both sources and truncate to exclude negative frequencies
% % as well as dc and fs/2.
% [S1,L] = stft(s1Padded,K);
% S1half = S1(2:(K+1)/2-1,:);
% S1halfreal = real(S1half);
% S1halfimag = imag(S1half);
% [S2,L2] = stft(s2Padded,K);
% S2half = S2(2:(K+1)/2-1,:);
% S2halfreal = real(S2half);
% S2halfimag = imag(S2half);
% 
% save('S1halfreal.txt','S1halfreal','-ASCII');
% save('S2halfreal.txt','S2halfreal','-ASCII');
% save('S1halfimag.txt','S1halfimag','-ASCII');
% save('S2halfimag.txt','S2halfimag','-ASCII');

% import data from file for testing with consistent data
fs = 14700;
K = 513;
L = 1149;
S1halfreal = importdata('S1halfreal.txt');
S1halfimag = importdata('S1halfimag.txt');
S1half = S1halfreal + j*S1halfimag;
S2halfreal = importdata('S2halfreal.txt');
S2halfimag = importdata('S2halfimag.txt');
S2half = S2halfreal + j*S2halfimag;


%% Place sensors
M = 100; % M = number of sensors
Nsrcs = 2; % Nsrcs = number of sources
spSize = 10; % spSize = size of the room (m)
space = [spSize, spSize, spSize]'; % Dimensions of the space
spcDim = length(space);
% Mloc = (rand(M,spcDim)*diag(space)).'; % Mloc = matrix containing 3d sensor locations
% sloc = ((rand(Nsrcs,spcDim)*diag(space))).';%+[0,0,2;0,0,2]).'; % sloc = matrix condaining 3d source locations
Mloc = importdata('Mloc.txt'); % For testing with consistent sensor placement
Mloc = Mloc(:,1:M); % Truncate unrequired sensors
sloc = importdata('sloc.txt'); % For testing with consistent source placement
% save('Mloc.txt', 'Mloc', '-ASCII');
% save('sloc.txt', 'sloc', '-ASCII');

% Calculate distances
ssd = zeros(Nsrcs,M);
for ns=1:Nsrcs
    for m=1:M
        ssd(ns,m) = norm(Mloc(:,m)-sloc(:,ns));
    end
end

% Display layout
figure; plot3(Mloc(1,:), Mloc(2,:), Mloc(3,:), '*'); grid on; hold on; 
plot3(sloc(1,1), sloc(2,1), sloc(3,1), 'o'); 
plot3(sloc(1,2), sloc(2,2), sloc(3,2), '^'); legend('Sensors','Target','Interferer')
set(gca, 'fontsize', 14);

%% Create ATFs
Khalf = (K-1)/2-1;
fdom = (fs/(K-1)) * [1:Khalf]';
c = 343; % c = speed of sound in m.s^-1
At = zeros(Khalf,M);
Ai = zeros(Khalf,M);
Atnogain = zeros(Khalf,M);
for m=1:M
    At(:,m) = exp(-1i*2*pi*fdom'*ssd(1,m)/c) / (4*pi*ssd(1,m)^2);
    Ai(:,m) = exp(-1i*2*pi*fdom'*ssd(2,m)/c) / (4*pi*ssd(2,m)^2);
    Atnogain(:,m) = exp(-1i*2*pi*fdom'*ssd(1,m)/c); 
end

%% Create observations
X = zeros(Khalf,L,M); Xt = zeros(Khalf,L,M); Xi = zeros(Khalf,L,M);
for l = 1:L
    for m = 1:M
        Xt(:,l,m) = At(:,m).*S1half(:,l); % These are used for calculating SNR 
        Xi(:,l,m) = Ai(:,m).*S2half(:,l); % These are used for calculating SNR 
        X(:,l,m) = Xt(:,l,m)+Xi(:,l,m);
    end
end

%% Delay and sum
W = At; % Atnogain
[yds,ydsSNRdb] = myBfOp(X,Xt,Xi,W);
ydsSNRdb

% dsSNRdb = mySnr();

% Covariance of noise and interference, should be I = identity
% Ri = zeros(Khalf,M,M);
% for k = 1:Khalf
%     for l = 1:L
%         Xtmp = squeeze(Xi(k,l,:));
%         Ri(k,:,:) = (1/L)*squeeze(Ri(k,:,:)) + Xtmp*Xtmp';
%     end
% end

% White noise gain
% G = zeros(Khalf,1);
% for k = 1:Khalf
%     Atmp = Atnogain(k,:).';
%     Wtmp = W(k,:).';
%     G(k) = (Atmp'*Wtmp*Wtmp'*Atmp)/(Wtmp'*Wtmp);
% 
% end
% Gsum = sum(real(G));
% 10*log10(Gsum)

%% MVDR optimum weights
dl = 1e-9; % dl = diagonal loading factor - ensures that the covariance is invertible
Wopt = myMvdrOpt(At,X,dl);
[yopt,yFrSNRdb] = myBfOp(X,Xt,Xi,Wopt);
yFrSNRdb

output = [M, ydsSNRdb, yFrSNRdb];
fid = fopen('SnrVsSensors.txt', 'at');
fprintf(fid, '%f %f %f\n', output);
fclose(fid);

% White noise gain
% G = zeros(Khalf,1);
% for k = 1:Khalf
%     Atmp = Atnogain(k,:).';
%     Wtmp = Wopt(k,:).';
%     G(k) = (Atmp'*Wtmp*Wtmp'*Atmp)/(Wtmp'*Wtmp);
% 
% end
% Gsum = sum(real(G));
% 10*log10(Gsum)


% Vorobyov principles of min...
% SINR = E{|w^H s|^2} / E{|w^H (i+n)|^2} = sigma_s^2 |w^H a(theta_s)|^2 /
% w^H R_i+n w where sigma_s^2 = E{|s(k)|^2} and R_i+n=E{(i(k)+n(k))(i(k)+n(k)^H)}
% sigma_S_sq = zeros(Khalf,1);
% Rin = zeros(Khalf,M,M);
% SINR = zeros(Khalf,1);
% for k=1:Khalf
%     for l=1:L
%         sigma_S_sq(k) = sigma_S_sq(k) + (1/L)*S1(k,l)'*S1(k,l);
%         Rin(k,:,:) = squeeze(Rin(k,:,:)) + (1/L)*(squeeze(Xi(k,l,:))*squeeze(Xi(k,l,:))');
%     end
%     SINR(k) = (sigma_S_sq(k)*(abs(conj(Wopt(k,:))*At(k,:).')^2))...
%         / (conj(Wopt(k,:))*squeeze(Rin(k,:,:))*Wopt(k,:).');
%     
%     SINRnosig(k) = (conj(Wopt(k,:))*At(k,:).')'*(conj(Wopt(k,:))*At(k,:).')...
%         / ((1/M)*(conj(Wopt(k,:))*Wopt(k,:).'));
% end

% Kelson
% WHA = zeros(Khalf,2,2);
% for k=1:Khalf
%     WHA(k,:,:) = [conj(Wopt(k,:));conj(Wopt(k,:))] * [At(k,:).', Ai(k,:).'] ;
%     
% end
% for l=1:L
%     Ytc(:,l) = WHA(:,1,1) .* S1half(:,l);
%     Yti(:,l) = WHA(:,1,2) .* S2half(:,l);
% end
% Ytc = [zeros(1,L);Ytc;zeros(2,L);conj(flipud(Ytc))];
% Yti = [zeros(1,L);Yti;zeros(2,L);conj(flipud(Yti))];
% ytc = myOverlapAdd(Ytc);
% yti = myOverlapAdd(Yti);
% SNRKelson = 10*log10((ytc'*ytc) / (yti'*yti))




%% Adaptive Frost MVDR
% mu = 200; % mu = step size 
% Iter = 10; % Iter = number of iterations per window
% [Y,W,Wmse] = myFrostAdapt(At,X,mu,Iter,Wopt);
% % [yFrA,yFrASNRdb] = myBfOp(X,Xt,Xi,W);
% 
% % Create two sided Y and take istft
% Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
% yFrO = myOverlapAdd(Y); % yFrO = y Frost Online
% figure; plot(yFrO);% 
% figure; semilogy(Wmse);
% figure; plot(Wmse);
% 
% %% Adaptive Frost using actual covariance
% mu = 200; % mu = step size 
% Iter = 10; % Iter = number of iterations per window
% [Y,W,Wmse] = myFrostAdaptTrueCov(At,X,mu,Iter,Wopt);
% 
% % Create two sided Y and take istft
% Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
% yFrO = myOverlapAdd(Y); % yFrO = y Frost Online
% figure; plot(yFrO);
% figure; semilogy(Wmse);
% figure; plot(Wmse);

%% Sparse distributed BiADMM MVDR
% Calculate the adjacency matrix
% Find sensor to sensor distances
% sensd = zeros(M,M);
% for m=1:M
%     for mm=m:M
%         sensd(m,mm) = norm(Mloc(:,m)-Mloc(:,mm));
%     end
% end
% sensd = sensd + triu(sensd,1).'; % Convert from upper triangular to symmetric
% 
% sensdtmp = sensd + diag(ones(1,M)*2*spSize); % 2* spSize is guaranteed to be larger than any distance between sensors
% Nnghbrs = 2; % Nnghbrs sets the minimum number of neighbors per node
% NNghbrsCheck = 10;
% N = -ones(NNghbrsCheck,M); % N = matrix containing indices of the neighbours to each node
% while min(min(N(1:Nnghbrs,:))) == -1
%     sensdmin = find(sensdtmp(:)==min(sensdtmp(:)));
%     [sensdminRow, sensdminCol] = ind2sub(size(sensdtmp),sensdmin(1)); % Find the row and column of the smallest distance left, note that there will be two values as the matrix is symmetric
%     if N(1,sensdminRow) == -1 || N(2,sensdminRow) == -1 || ... % You only get to add a neighbor if you have less than two already
%             N(1,sensdminCol) == -1 || N(2,sensdminCol) == -1 
%         for aa = 1:NNghbrsCheck
%             if N(aa,sensdminRow) == -1
%                 N(aa,sensdminRow) = sensdminCol; % Add the neighbor to the first empty slot
%                 break
%             end
%         end
%         for aa = 1:NNghbrsCheck
%             if N(aa,sensdminCol) == -1
%             N(aa,sensdminCol) = sensdminRow; % Add the neighbor to the first empty slot
%                 break
%             end
%         end
%     end
%     
%     % Set both copies of the distance just used higher than the possible
%     % distance, so they cannot be the minimum distance again.
%     sensdtmp(sensdminRow,sensdminCol) = 2*spSize; 
%     sensdtmp(sensdminCol,sensdminRow) = 2*spSize;
% end
% % Truncate unused rows of N
% maxNoNghbrs = min(find(sum(N,2)==-M))-1;
% N = N(1:maxNoNghbrs,:);



