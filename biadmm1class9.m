%% Start again mvdr beamformer, hopefully distributed using biadmm. 

% % % History % % %

% biadmm1class6.m is a working realnodeonly version. It has issues though,
% for example, Wopt is worse than listening to the very close microphone.
% This suggests there is a problem with the covariance.

% biadmm1class7.m uses the spatial covariance only for the first time,
% hopefully that clears up a few issues with W and Wopt.

% biadmm1class8.m uses the setup from 'on simplifying pdmm...' zhang,
% specifically in the setup of the Aij matrices and the use of uij for
% setting the sign. 

% biadmm1class9.m is just a tidy up for showing Aryan

close all; clear;

% Import target audio
Nsrcs = 1; % Nsrcs = number of sources
s = cell(Nsrcs,1);
AudioFileNames = {'422-122949-0013.flac';'2078-142845-0005.flac'};
for ns=1:Nsrcs
    s{ns} = audioread(strcat('/audio/',AudioFileNames{ns})); 
end
fs = 16e3;

% Truncate to desired length, ensuring that the length is a multiple of 
% the window length.
K = 2^12+1; % K = window length in samples, and the number of frequency bins
Khalf = (K-1)/2-1;
tls = 5; % tls = target length in seconds
tl = tls*fs-mod(tls*fs,K-1)+1; % tl = target length in samples, adjusted for window length and sampling frequency
for ns=1:Nsrcs
    s{ns} = s{ns}(1:tl);
end

%% FFT
S = cell(Nsrcs,1);
for ns=1:Nsrcs
    S{ns} = fft(s{ns});
    S{ns} = S{ns}(2:(tl-1)/2); % Truncate to half spectrum
end

%% Place sensors
M = 3; % M = number of sensors

% Create nodes
node = cell(M,1); % 2*M
for m=1:M%2*M
    node{m} = myNode;
end

spSize = 1; % spSize = size of the room (m)
space = [spSize, spSize, spSize]'; % Dimensions of the space
spcDim = length(space);
% xloc = (rand(M,spcDim)*diag(space)).'; % xloc = matrix containing 3d sensor locations
% sloc = ((rand(Nsrcs,spcDim)*diag(space))).';%+[0,0,2;0,0,2]).'; % sloc = matrix containing 3d source locations
% sloc =   spSize*[0.1,0.92;
%                 0.1,0.92;
%                 0.1,0.92];
sloc =   spSize*[0.1;
                0.1;
                0.1];
xloc(:,1:3) = spSize*[0.11,0.3,0.91;
                      0.11,0.4,0.91;
                      0.11,0.5,0.91];
    
% Set location for each node
for m=1:M
    node{m}.loc = xloc(:,m);
end

% Calculate distances
ssd = myGetDist(xloc,sloc);

%% Display layout
myDisplayLayout(xloc,sloc);

%% Create ATFs and observations for full fft version
fdom = (fs/(tl-1)) * (1:(tl-1)/2-1);
c = 343;
L = (length(s{1}(1:end-1))/(K-1))*2+1;
X = zeros(Khalf,L,M);
xsave = zeros(length(s{1}),M);
for m=1:M
    Xfft = zeros((tl-1)/2-1,1);
    for ns=1:Nsrcs
        A = exp(-1i*2*pi*fdom.'*ssd(ns,m)/c) / (4*pi*ssd(ns,m));
        Xfft = Xfft + (A .* S{ns});
    end
    Xfft = [0;Xfft;0;0;conj(flipud(Xfft))];
    x = ifft(Xfft) + 0.0001*randn(tl,1);
    xsave(:,m) = x;
    xPadded = [zeros((K-1)/2,1);x(1:end-1);zeros((K-1)/2,1)];
    XTmp = stft(xPadded,K);
    X(:,:,m) = XTmp(2:(K-1)/2,:);
end

%% Find sensor to sensor distances
sensd = myFindSensDist(xloc);

%% Find neighbors, everyone within 0.5*spSize (Aiming for >5 neighbors each, but depends on the number of nodes as well)
Nneighbors = zeros(M,1);
for m=1:M
    node{m}.N = [find(sensd(:,m)<2*spSize) ];
    node{m}.Nlen = length(node{m}.N);
    Nneighbors(m) = node{m}.Nlen;    
end
fprintf('The minimum number of neighbors was %d. \nThe mean number of neighbors was %d. \n\n',min(Nneighbors),mean(Nneighbors));

%% Covariance based only on sensor and source location
fdomShort = (fs/(K-1)) * (1:((K-1)/2)-1);
R = zeros(Khalf,M,M);
for k=1:Khalf
    for ii=1:M
        for jj=1:M
            temp1 = norm(xloc(:,ii)-sloc(:,1));
            temp2 = norm(xloc(:,jj)-sloc(:,1));
            R(k,ii,jj) = exp(-1i*2*pi*fdomShort(k)*(temp1-temp2)/c)/((4*pi)^2 * temp1 * temp2);
        end
    end
    R(k,:,:) = squeeze(R(k,:,:)) + 1e-3*eye(M);
end

%% Initialization
for m=1:M
    % Initialize Lambdas and weights
    node{m}.L = ones(Khalf,2,node{m}.Nlen); % These are for node m's real neighbors including itself
    node{m}.W = ones(Khalf,node{m}.Nlen); % Note that this includes node m itself, node m's real neighbors, and node m's virtual neighbor (m+M)
 
    % Initialize Amn for all nodes
    Amn = cell(node{m}.Nlen,1); 
    for n = 1:node{m}.Nlen
        if m==node{m}.N(n)
            Amn{n} = zeros(2,node{m}.Nlen);
        else
            Amn{n} = double([(node{m}.N==m).';-(node{m}.N==node{m}.N(n)).']);
        end
    end
    node{m}.Amn = Amn;
   
    % Save look direction d for node m
    node{m}.d = exp(-1i*2*pi*fdomShort.'*ssd(1,m)/c) / (4*pi*ssd(1,m));
end

% Initialize output Y
Y = zeros(Khalf,L);

%% A1) Frequency domain Frost optimum weights
d = zeros(Khalf,M);
for m=1:M
    d(:,m) = node{m}.d;
end
Wopt = zeros(Khalf,M);
dk = zeros(Khalf,M);
for k=1:Khalf
    dk = d(k,:).';
    Wopt(k,:) = (squeeze(R(k,:,:))\dk)/(dk'/squeeze(R(k,:,:))*dk); 
end

% Find output using optimum weights
Yopt = zeros(Khalf,L);
for l=1:L
    Xtmp = squeeze(X(:,l,:));
    Yopt(:,l) = sum(conj(Wopt).*Xtmp,2);
end
Yopt = [zeros(1,L);Yopt;zeros(2,L);conj(flipud(Yopt))];
yopt = myOverlapAdd(Yopt);

%% Adaptive algorithm (new based on biadmm_1bin2.m)
Ltmp = 20; % For shorter run times
bin = 53;
ITER1 = 1;
ITER2 = 1;
Wsave = zeros(Khalf,Ltmp,M);

% Initialize W to Wopt
for m=1:M
    node{m}.W = Wopt;     
end

% % Store Wopt in Wsave (for single freq bin opt)
% for l=1:Ltmp
%     Wsave(:,l,:) = Wopt;
% end

for l=1:Ltmp
    for iter1=1:ITER1
        for k=1:Khalf
            for m=1:M
                for iter2=1:ITER2
                    [l,k,m]
                    Nlen = node{m}.Nlen;
                    AA = zeros(Nlen);
                    ALAW = zeros(Nlen,1);
                    dm = zeros(Nlen,1);
                    dm(m) = node{m}.d(k);

                    % W update
                    for n=1:Nlen
                        Amn = node{m}.Amn{n};
                        AA = AA + (Amn.'*Amn);
                        Lnm = node{node{m}.N(n)}.L(k,:,node{node{m}.N(n)}.N==m).';
                        Anm = flipud(node{node{m}.N(n)}.Amn{node{node{m}.N(n)}.N==m});
                        Wn = node{node{m}.N(n)}.W(k,:).';
                        ALAW = ALAW + (Amn.'*(Lnm-Anm*Wn));
                    end
                    node{m}.W(k,:) = (AA + squeeze(R(k,:,:)))\(ALAW + dm);
                    
                    % Lambda update
                    for n=1:Nlen
                        Amn = node{m}.Amn{n};
                        Anm = flipud(node{node{m}.N(n)}.Amn{node{node{m}.N(n)}.N==m});
                        Wn = node{node{m}.N(n)}.W(k,:).';
                        Wm = node{m}.W(k,:).';
                        node{m}.L(k,:,n) = node{n}.L(k,:,m).' - 0.7*(Anm*Wn + Amn*Wm);
                    end
                end
            end
        end
        
        % Save the current weights 
        Wtmp = zeros(Khalf,M);
        for m=1:M
            Wtmp(:,m) = node{m}.W(:,node{m}.N==m);
        end          
        Wsave(:,l,:) = Wtmp;
        
    end
    
    % Calculate output Y
    Y(:,l) = (1/M)*sum(squeeze(conj(Wsave(:,l,:))).*squeeze(X(:,l,:)),2);
    
end


%% Calculate BF output
Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
y = myOverlapAdd(Y);
figure; plot(y);

%% MSE between W and Wopt
WWoptMSE = zeros(Ltmp,1);
for l=1:Ltmp
    Wtmp = squeeze(Wsave(:,l,:));
    WWoptMSE(l) = mean(mean((Wtmp-Wopt).*conj(Wtmp-Wopt)));
end
figure; semilogy(WWoptMSE); grid on; title('WWoptMSE');

%% W vs Wopt full spectrum
figure; imagesc(abs(Wopt)); title('Wopt');
figure; imagesc(abs(squeeze(Wsave(:,Ltmp,:))));title('Wsave');

%% W vs Wopt single bin
% figure; imagesc(abs(squeeze(Wopt(bin,:)))); title('Wopt');
% figure; imagesc(abs(squeeze(Wsave(bin,Ltmp,:)).'));title('Wsave');

%% MSE between W and Wopt single bin
WSaveTmp = squeeze(Wsave(bin,:,:));
WOptTmp = squeeze(Wopt(bin,:));
LL = length(WSaveTmp);
WWoptMSE = zeros(LL,1);
for l=1:LL
    WWoptMSE(l) = mean((WOptTmp-WSaveTmp(l,:))*(WOptTmp-WSaveTmp(l,:))');
end
figure ; semilogy(WWoptMSE); grid on; title('WWoptMSE single bin');

%% Variance of the sensor weights
VarWsave = zeros(M,1);
VarWopt = zeros(M,1);
for m=1:M
    VarWsave(m) = Wsave(:,2,m)'*Wsave(:,2,m);
    VarWopt(m) = Wopt(:,m)'*Wopt(:,m);
end

figure; plot(VarWsave,'*--'); grid on; hold on; plot(VarWopt,'o--'); title('Variance in sensors'); legend('VarWsave','VarWopt');
% ratio = max(VarWsave)/max(VarWopt);
% figure; plot((1/ratio)*VarWsave,'*--'); grid on; hold on; plot(VarWopt,'o--'); title('Variance in sensors scaled'); legend('VarWsave','VarWopt');


%% Print setup for records
fprintf('Nsrcs = %d, K = %d, tls = %d, M = %d, spSize = %d, bin = %d, Ltmp = %d\n\n',Nsrcs,K,tls,M,spSize,bin,Ltmp);

%%
% Wopt(54,1)'*node{1}.d(54)+Wopt(54,2)'*node{2}.d(54)
% Wopt(123,1)'*node{1}.d(123)+Wopt(123,2)'*node{2}.d(123)

