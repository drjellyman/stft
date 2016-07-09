%% Start again mvdr beamformer, hopefully distributed using biadmm. Needs 
% to have 50 sensors randomly placed within a 100x100x100 m free space. One
% target speech signal, and one noise signal, also randomly placed. Signal
% of interest is a 20 s speech signal chosen randomly from a 60 s file. fs
% = 16 kHz. Window length is 25 ms with a Hann window and 50% overlap.
% Interference is a randomly placed, zero meaqn gaussian point source with
% power equal to -5, 0, 5 dB when compared to the target signal. 
close all; clear;

% Import target audio
[s1,fs1] = audioread('273177__xserra__la-vaca-cega-eva.wav'); 

% Downsample
fsd = 16e3; % fsd = desired sampling frequency
fs = fs1/round(fs1/fsd); % fs = actual sampling frequency post resample
s1 = resample(s1,1,round(fs1/fsd));

% Truncate to desired length, ensure that the length is a multiple of 
% the window length, and randomly select a section of the audio file.
K = 2^9+1; % K = window length in samples, and the number of frequency bins
tls = 5; % tls = target length in seconds
tl = tls*fs-mod(tls*fs,K-1); % tl = target length in samples, adjusted for window length and sampling frequency
start = floor((length(s1)-tl)*rand);
s1 = s1(start:start+tl-1,1); % Truncate s1 to be one channel, and 20 s long

% Normalize the target audio file to make it easy to change files
s1rms = rms(s1);
s1rmsinv = 1./s1rms;
s1 = s1 * (0.1*diag(s1rmsinv));% This should probably be scaled down to avoid clipping

% Set up interferer with equal power, i.e. snr = 0 dB 
s1Pow = (s1'*s1) / length(s1);
s2 = sqrt(s1Pow) *  randn(length(s1),1); % Currently set up for equal power
s2Pow = (s2'*s2) / length(s2);
SourceSNRdB = 10*log10(s1Pow/s2Pow);

%% STFT 

% Zero the target for testing
% s1 = zeros(length(s1),1);
% Zero the interferer for testing
% s2 = zeros(length(s2),1);

% pad the source signals so the 1st half window doesn't distort the data
s1Padded = [zeros((K-1)/2,1);s1;zeros((K-1)/2,1)];
s2Padded = [zeros((K-1)/2,1);s2;zeros((K-1)/2,1)];

% Take stft of both sources and truncate to exclude negative frequencies
% as well as dc and fs/2.
[S1,L] = stft(s1Padded,K);
S1half = S1(2:(K+1)/2-1,:);
S1halfreal = real(S1half);
S1halfimag = imag(S1half);
[S2,L2] = stft(s2Padded,K);
S2half = S2(2:(K+1)/2-1,:);
S2halfreal = real(S2half);
S2halfimag = imag(S2half);
% 
% save('S1halfreal.txt','S1halfreal','-ASCII');
% save('S2halfreal.txt','S2halfreal','-ASCII');
% save('S1halfimag.txt','S1halfimag','-ASCII');
% save('S2halfimag.txt','S2halfimag','-ASCII');

% % import data from file for testing with consistent data
% fs = 14700;
% K = 513;
% L = 1149;
% S1halfreal = importdata('S1halfreal.txt');
% S1halfimag = importdata('S1halfimag.txt');
% S1half = S1halfreal + j*S1halfimag;
% S2halfreal = importdata('S2halfreal.txt');
% S2halfimag = importdata('S2halfimag.txt');
% S2half = S2halfreal + j*S2halfimag;



%% Place sensors
M = 50; % M = number of sensors

% Create nodes
node = cell(2*M,1);
for m=1:2*M
    node{m} = myNode;
end

Nsrcs = 2; % Nsrcs = number of sources
spSize = 10; % spSize = size of the room (m)
space = [spSize, spSize, spSize]'; % Dimensions of the space
spcDim = length(space);
Mloc = (rand(M,spcDim)*diag(space)).'; % Mloc = matrix containing 3d sensor locations
sloc = ((rand(Nsrcs,spcDim)*diag(space))).';%+[0,0,2;0,0,2]).'; % sloc = matrix containing 3d source locations
% Mloc = importdata('Mloc.txt'); % For testing with consistent sensor placement
% Mloc = Mloc(:,1:M); % Truncate unrequired sensors
% sloc = importdata('sloc.txt'); % For testing with consistent source placement
% save('Mloc.txt', 'Mloc', '-ASCII');
% save('sloc.txt', 'sloc', '-ASCII');

% Set location for each node
for m=1:2*M
    if m<M+1
        node{m}.loc = Mloc(:,m);
    elseif m>M
        node{m}.loc = 0;
    end
end

% Calculate distances
ssd = zeros(Nsrcs,M); % ssd = source to sensor distances
for ns=1:Nsrcs
    for m=1:M
        ssd(ns,m) = norm(Mloc(:,m)-sloc(:,ns));
    end
end

%% Display layout
% figure; plot3(Mloc(1,:), Mloc(2,:), Mloc(3,:), '*'); grid on; hold on; 
% plot3(sloc(1,1), sloc(2,1), sloc(3,1), 'o'); 
% plot3(sloc(1,2), sloc(2,2), sloc(3,2), '^'); legend('Sensors','Target','Interferer')
% set(gca, 'fontsize', 14);

%% Create ATFs
Khalf = (K-1)/2-1;
fdom = (fs/(K-1)) * (1:Khalf)';
c = 343; % c = speed of sound in m.s^-1
At = zeros(Khalf,M);
Ai = zeros(Khalf,M);
Atnogain = zeros(Khalf,M);
for m=1:M
    At(:,m) = exp(-1i*2*pi*fdom'*ssd(1,m)/c) / (4*pi*ssd(1,m)^2);
    Ai(:,m) = exp(-1i*2*pi*fdom'*ssd(2,m)/c) / (4*pi*ssd(2,m)^2);
    Atnogain(:,m) = exp(-1i*2*pi*fdom'*ssd(1,m)/c); 
    
    node{m}.d = At(:,m); % Store each nodes Atf
    node{m+M}.d = 0; % The virtual nodes have no need for an Atf
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
whichNode = 24;
Xtcheck = [zeros(1,L);Xt(:,:,whichNode);zeros(2,L);conj(flipud(Xt(:,:,whichNode)))];
mySpectrogram(Xtcheck);
xtcheck = myOverlapAdd(Xtcheck);
figure; plot(xtcheck);

Xcheck = [zeros(1,L);X(:,:,whichNode);zeros(2,L);conj(flipud(X(:,:,whichNode)))];
mySpectrogram(Xcheck);
xcheck = myOverlapAdd(Xcheck);
figure; plot(xcheck);

return

%% Create true covariance matrix
% Rtrue = cell(Khalf,1);
% for l=1:L
%     for k=1:Khalf
%         % Initialize Rtrue on first time
%         if l==1
%             Rtrue{k} = zeros(M);
%         end
%         Rtrue{k} = Rtrue{k} + squeeze(X(k,l,:))*squeeze(X(k,l,:)).';
%     end
% end
% rcond(Rtrue{17})
% 
% % Alternative calculation
% for k=1:Khalf
%     XTmp = squeeze(X(k,:,:));
%     Rtrue2(k,:,:) = XTmp' * XTmp;
% end
% Rtrue2mean = squeeze(mean(Rtrue2,1));
% rcond(Rtrue2mean)

%% Calculate noise only, target only, and combined covariances over all time
% Ri = cell(Khalf,1);
% Rt = cell(Khalf,1);
R = cell(Khalf,1);
for k=1:Khalf
    for l=1:L
        if l==1
%             Ri{k} = zeros(M);
%             Rt{k} = zeros(M);
            R{k} = zeros(M);
        end
%         Ri{k} = Ri{k} + (1/L)*squeeze(Xi(k,l,:))*squeeze(Xi(k,l,:))';
%         Rt{k} = Rt{k} + (1/L)*squeeze(Xt(k,l,:))*squeeze(Xt(k,l,:))';
        R{k} = R{k} + (1/L)*squeeze(X(k,l,:))*squeeze(X(k,l,:))';
    end
end

%% Sparse distributed BiADMM MVDR Matt
% Calculate the adjacency matrix
% Find sensor to sensor distances
sensd = zeros(M,M); % sensd = sensor to sensor distance
for m=1:M
    for mm=m:M
        sensd(m,mm) = norm(Mloc(:,m)-Mloc(:,mm));
    end
end
sensd = sensd + triu(sensd,1).'; % Convert from upper triangular to symmetric

%% Find neighbors, closest first
% sensdtmp = sensd + diag(ones(1,M)*2*spSize); % 2* spSize is guaranteed to be larger than any distance between sensors
% Nnghbrs = 5; % Nnghbrs sets the minimum number of neighbors per node
% NNghbrsCheck = 20;
% N = -ones(NNghbrsCheck,M); % N = matrix containing indices of the neighbours to each node
% C = zeros(M,M);
% while min(min(N(1:Nnghbrs,:))) == -1
%     sensdmin = find(sensdtmp(:)==min(sensdtmp(:)));
%     [sensdminRow, sensdminCol] = ind2sub(size(sensdtmp),sensdmin(1)); % Find the row and column of the smallest distance left, note that there will be two values as the matrix is symmetric
%     % You only get to add a neighbor if you have less than Nnghbrs already
%     if N(1,sensdminRow) == -1 || N(2,sensdminRow) == -1 || N(3,sensdminRow) == -1 || N(4,sensdminRow) == -1 || N(5,sensdminRow) == -1 ||... 
%             N(1,sensdminCol) == -1 || N(2,sensdminCol) == -1 || N(3,sensdminCol) == -1 || N(4,sensdminCol) == -1 || N(4,sensdminCol) == -1
%         for aa = 1:NNghbrsCheck
%             if N(aa,sensdminRow) == -1
%                 N(aa,sensdminRow) = sensdminCol; % Add the neighbor to the first empty slot
%                 break
%             end
%         end
%         for aa = 1:NNghbrsCheck
%             if N(aa,sensdminCol) == -1
%                 N(aa,sensdminCol) = sensdminRow; % Add the neighbor to the first empty slot
%                 break
%             end
%         end
%         C(sensdminRow,sensdminCol) = 1;
%         C(sensdminCol,sensdminRow) = 1;
%     end    
%     % Set both copies of the distance just used higher than the possible
%     % distance, so they cannot be the minimum distance again.
%     sensdtmp(sensdminRow,sensdminCol) = 2*spSize; 
%     sensdtmp(sensdminCol,sensdminRow) = 2*spSize;
% end
% 
% % Truncate unused rows of N
% maxNoNghbrs = find(sum(N,2)==-M,1)-1;
% N = N(1:maxNoNghbrs,:);

% Save neighbor lists to nodes
% for m=1:M*2
%     if m<M+1
%         node{m}.N = [m;N(1:find(N(:,m)==-1,1)-1,m);m+M];        
%     elseif m>M
%         node{m}.N = [m;m-M]; % A virtual node has only two neighbors - itself and its real node
%     end
%     node{m}.Nlen = length(node{m}.N);
%     if node{m}.Nlen < Nnghbrs && m <= M
%         fprintf('Warning: Node %d has only %d neighbors. \n',m,node{m}.Nlen);
%     end
% end

%% Find neighbors, everyone within 0.6*spSize (Aiming for >5 neighbors each)
Nneighbors = zeros(M,1);
for m=1:M
    node{m}.N = [find(sensd(:,m)<0.6*spSize) ; m+M];
    node{m}.Nlen = length(node{m}.N);
    Nneighbors(m) = node{m}.Nlen;
    
    node{m+M}.N = [m+M;m];
    node{m+M}.Nlen = 2;
end
fprintf('The minimum number of neighbors was %d. \nThe mean number of neighbors was %d. \n\n',min(Nneighbors),mean(Nneighbors));




%% Updates
% Initialize Lambdas and weights
for m=1:M*2
    node{m}.L = zeros(Khalf,2,node{m}.Nlen-1);
    node{m}.W = ones(Khalf,node{m}.Nlen);
end

% Initialize Amn, Anm, and d
for m=1:2*M
    Amn = cell(node{m}.Nlen,1); % zeros(2,node{m}.Nlen);
%     Anm = cell(node{m}.Nlen,1); % zeros(2,node{node{m.N(m)}}.Nlen);
    for n = 1:node{m}.Nlen
        Amn{n} = myAconsistency(node{m}.N,node{node{m}.N(n)}.N);            
    end
    node{m}.Amn = Amn;
%     node{m}.Anm = Anm;
%     node{m}.Nlen = length(node{m}.N);

    % Save look direction for node m
    if m<M+1
        node{m}.d = At(:,m);
    end
end

% Step through window indexes (i.e. blocks of time)

% Initialize regularization parameter alpha
alpha = 1;

% Initialize output Y
Y = zeros(Khalf,L);

for l=1:10
   
    % Step through real nodes for primal/dual update
    for m=1:M 
        
        % iter sets the number of iterative updates per window
        for iter=1
        
            % Get observations
            for n=1:node{m}.Nlen-1
                node{m}.X(:,n) = X(:,l,node{m}.N(n));       
            end

%             % Calculate covariance using local observations
%             for k=1:Khalf
%                 node{m}.R{k} = node{m}.X(k,:)'*node{m}.X(k,:);
%             end

            % Assign covariance from actual covariance calculated above
            if l==1
                for k=1:Khalf
                    RTemp = zeros(node{m}.Nlen-1);
                    for n=1:node{m}.Nlen-1
                        RTemp(:,n) = R{k}(m,node{m}.N(1:end-1)).';
                        RTemp(n,:) = R{k}(m,node{m}.N(1:end-1));
                    end                
                    node{m}.R{k} = RTemp;
                end
            end

            % Calculate primal update for real node
            AARsum = zeros(node{m}.Nlen-1,node{m}.Nlen-1);
            ALAWdsum = zeros(node{m}.Nlen-1,1);
            AARIsum = zeros(node{m}.Nlen-1,node{m}.Nlen-1);
            ALAWARdsum = zeros(node{m}.Nlen-1,1);
            WRealNew = ones(Khalf,node{m}.Nlen);
            WVirtNew = ones(Khalf);
            LambdaRealNew = ones(Khalf,2,node{m}.Nlen);        
            LambdaVirtNew = ones(Khalf,2);
            LambdaNew = ones(Khalf,2,M);

            bk = zeros(Khalf,1);
            for k=1:Khalf
                dTmp =  [node{m}.d(k);zeros(node{m}.Nlen-2,1)];
                for n=1:node{m}.Nlen-1 % -1 because I want to exclude the virtual node for now

                    % for W update
                    AmnTmp = node{m}.Amn{n}(:,1:node{m}.Nlen-1);
                    AARsum = AARsum+(AmnTmp.'*AmnTmp+node{m}.R{k});
                    LambdaTmp = node{node{m}.N(n)}.L(k,:,find(node{node{1}.N(1)}.N==1,1)).';  
                    AnmTmp = -node{node{m}.N(n)}.Amn{(node{node{m}.N(n)}.N==m)};% node{m}.Anm{n}(:,1:end-1);
                    WTmp = node{node{m}.N(n)}.W(k,1:end).';
                    ALAWdsum = ALAWdsum + (AmnTmp.'*(LambdaTmp-AnmTmp*WTmp)+node{m}.d(k));

                    % For Lambda update, zk
                    AARIsum = AARIsum + (AmnTmp.'*AmnTmp/node{m}.R{k}+eye(node{m}.Nlen-1));   
                    ALAWARdsum = ALAWARdsum + (AmnTmp.'*(LambdaTmp-AnmTmp*WTmp-AmnTmp*node{m}.R{k}*dTmp));                
                end
    %             node{m}.W(k,1:node{m}.Nlen-1) = (AARsum+node{m}.R{k})\(ALAWdsum);
                WRealNew(k,1:node{m}.Nlen-1) = (AARsum+node{m}.R{k})\(ALAWdsum);
                zk = AARIsum\ALAWARdsum;

                % Lambda update requires twice around the neighboring nodes:
                % once for finding the summations in zk, and once for Lambda
                % itself.
                
                for n=1:node{m}.Nlen-1
                    LambdaTmp = node{node{m}.N(n)}.L(k,:,find(node{node{1}.N(1)}.N==1,1)).';  
                    AnmTmp = -node{node{m}.N(n)}.Amn{(node{node{m}.N(n)}.N==m)};  %AnmTmp = node{m}.Anm{n}(:,1:end-1);
                    WTmp = node{node{m}.N(n)}.W(k,1:end).';
                    AmnTmp = node{m}.Amn{n}(:,1:node{m}.Nlen-1);                
    %                 node{m}.L(k,:,n) = LambdaTmp-AnmTmp*WTmp-AmnTmp*node{m}.R{k}*(dTmp+zk);
                    LambdaNew(k,:,n) = LambdaTmp-AnmTmp*WTmp-AmnTmp*node{m}.R{k}*(dTmp+zk);
                end 

                % Primal update for virtual node m+M
                WTmp = node{m}.W(k,:);
                AnmTmp = -node{m}.Amn{node{m}.Nlen};
                AmnTmp = node{m+M}.Amn{1}(:,2);
                LambdaTmp = node{m}.L(k,:,end);
                bk = 2*WTmp*AnmTmp.'*AmnTmp-AmnTmp.'*LambdaTmp.'; 
                WVirtNew(k) = (-bk+sign(bk)*min(abs(bk),alpha))/(2);

                % Dual update for virtual node m+M
                LambdamnTmp = node{m+M}.L(k,:);
                LambdaVirtNew(k,:) = 2*LambdaTmp*LambdamnTmp.'-WTmp*AnmTmp.';
                if LambdaVirtNew(k,1) > alpha
                    LambdaVirtNew(k,1) = alpha; end
                if LambdaVirtNew(k,2) > alpha
                    LambdaVirtNew(k,2) = alpha; end
                if LambdaVirtNew(k,1) < -alpha
                    LambdaVirtNew(k,1) = -alpha; end
                if LambdaVirtNew(k,1) < -alpha
                    LambdaVirtNew(k,1) = -alpha; end                
            end   

            % Assign new W and L values to node m and m+M
            node{m}.W = WRealNew;
            node{m+M}.W = WVirtNew;
            node{m}.L = LambdaRealNew;
            node{m+M}.L = LambdaVirtNew;
        end
    end  
    
    % Generate output using updated primal and dual weights
    for m=1:M
        Y(:,l) = Y(:,l) + sum(node{m}.W(:,1:end-1).*node{m}.X,2);
    end    
end


%% Calculate BF output
Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
mySpectrogram(Y);
y = myOverlapAdd(Y);
figure; plot(y);

















