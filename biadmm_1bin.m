%% Start again mvdr beamformer, hopefully distributed using biadmm. Needs 
% to have 50 sensors randomly placed within a 100x100x100 m free space. One
% target speech signal, and one noise signal, also randomly placed. Signal
% of interest is a 20 s speech signal chosen randomly from a 60 s file. fs
% = 16 kHz. Window length is 25 ms with a Hann window and 50% overlap.
% Interference is a randomly placed, zero meaqn gaussian point source with
% power equal to -5, 0, 5 dB when compared to the target signal. 
close all; clear;

% Import target audio
Nsrcs = 2; % Nsrcs = number of sources
s = cell(Nsrcs,1);

AudioFileNames = {'422-122949-0013.flac';'2078-142845-0005.flac'};
for ns=1:Nsrcs
    s{ns} = audioread(strcat('/audio/',AudioFileNames{ns})); 
end
fs = 16e3;
% s{1} = randn(length(s{2}),1);
% s{2} = randn(length(s{2}),1);

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
for ns=1:Nsrcs
    S{ns} = fft(s{ns});
    S{ns} = S{ns}(2:(tl-1)/2); % Truncate to half spectrum
end

%% Place sensors
M = 3; % M = number of sensors

% Create nodes
node = cell(2*M,1);
for m=1:M%2*M
    node{m} = myNode;
end

spSize = 10; % spSize = size of the room (m)
space = [spSize, spSize, spSize]'; % Dimensions of the space
spcDim = length(space);
% Mloc = (rand(M,spcDim)*diag(space)).'; % Mloc = matrix containing 3d sensor locations
% sloc = ((rand(Nsrcs,spcDim)*diag(space))).';%+[0,0,2;0,0,2]).'; % sloc = matrix containing 3d source locations
sloc =   [1,9.2;
          1,9.2;
          1,9.2];
Mloc = [1.1,3,9.1;
        1.1,4,9.1;
        1.1,5,9.1];


% Set location for each node
for m=1:M
    node{m}.loc = Mloc(:,m);
%     node{m+M}.loc = 0; % Virtual nodes don't have a real location
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
% plot3(sloc(1,2:end), sloc(2,2:end), sloc(3,2:end), '^'); legend('Sensors','Target','Interferer')
% set(gca, 'fontsize', 14);

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
    x = ifft(Xfft) + 0.001*randn(tl,1);
    xsave(:,m) = x;
    xPadded = [zeros((K-1)/2,1);x(1:end-1);zeros((K-1)/2,1)];
    XTmp = stft(xPadded,K);
    X(:,:,m) = XTmp(2:(K-1)/2,:);
end

%% Find sensor to sensor distances
sensd = zeros(M,M); % sensd = sensor to sensor distance
for m=1:M
    for mm=m:M
        sensd(m,mm) = norm(Mloc(:,m)-Mloc(:,mm));
    end
end
sensd = sensd + triu(sensd,1).'; % Convert from upper triangular to symmetric

%% Find neighbors, everyone within 0.5*spSize (Aiming for >5 neighbors each, but depends on the number of nodes as well)
Nneighbors = zeros(M,1);
for m=1:M
    node{m}.N = [find(sensd(:,m)<2*spSize) ];
    node{m}.Nlen = length(node{m}.N);
    Nneighbors(m) = node{m}.Nlen;    
end
fprintf('The minimum number of neighbors was %d. \nThe mean number of neighbors was %d. \n\n',min(Nneighbors),mean(Nneighbors));

%% Calculate covariances over all time
% Rtmp = cell(Khalf,1);
% R = zeros(node{m}.Nlen);
% for k=1:Khalf
%     for l=1:L
%         if l==1
%             R{k} = zeros(M);
%         end
%         XTmp = squeeze(X(k,l,:));
%         Rtmp{k} = Rtmp{k} + (1/L)*(XTmp*XTmp');
%     end
% end
% for k=1:Khalf
%     R = R + Rtmp{k};
% end
R = zeros(M);
for k=1:Khalf
    for l=1:L
        Xtmp = squeeze(X(k,l,:));
        R = R + (1/(L*Khalf))*(Xtmp*Xtmp');
    end
end
R = R + 1*eye(3);
rcond(R)

%% Initialize
bin=49; % roughly 250 Hz
f=bin*(fs/(K-1));

% % Covariance (over all time)
% Xtmp = sum(squeeze(X(bin,:,:)),1).';
% R = Xtmp*Xtmp' + 0*eye(node{m}.Nlen);
% Rrcond = rcond(R)

for m=1:M
    node{m}.L = ones(2,node{m}.Nlen);
    node{m}.W = ones(node{m}.Nlen,1);
        
    % Initialize Amn for all nodes
    Amn = cell(node{m}.Nlen,1); 
    for n = 1:node{m}.Nlen
        if m==node{m}.N(n)
            Amn{n} = zeros(2,node{m}.Nlen);
        else
            Amn{n} = [(node{m}.N==m).';-(node{m}.N==node{m}.N(n)).'];
        end
    end
    node{m}.Amn = Amn;
   
    % Save look direction d for node m
    node{m}.d = exp(-1i*2*pi*f.'*ssd(1,m)/c) ;%/ (4*pi*ssd(1,m));
end

Y = zeros(L,1);

%% Adaptive algorithm
Ltmp = 25; % For shorter run times
ITER1 = 1;
ITER2 = 1;
Wsave = zeros(Ltmp,M);
for l=1:Ltmp
    for iter1=1:ITER1
        for k=bin
            for m=1:M
                for iter2=1:ITER2
                    Nlen = node{m}.Nlen;
                    AApR = zeros(Nlen);
                    ALAWD = zeros(Nlen,1);
                    AARpI = zeros(Nlen);
                    ALAWARD = zeros(Nlen,1);
                    for n=1:Nlen
                        Amn = node{m}.Amn{n};
                        AApR = AApR + (Amn.'*Amn+R);
                        Lnm = node{node{m}.N(n)}.L(:,node{node{m}.N(n)}.N==m);
                        Anm = flipud(node{node{m}.N(n)}.Amn{node{node{m}.N(n)}.N==m});
                        Wn = node{node{m}.N(n)}.W;
                        dm = zeros(Nlen,1);
                        dm(m) = node{m}.d;
                        ALAWD = ALAWD + (Amn.'*(Lnm-Anm*Wn)+dm);
                        AARpI = AARpI + (Amn.'*Amn/R+eye(Nlen));
                        ALAWARD = ALAWARD + (Amn.'*(Lnm-Anm*Wn-Amn/R*dm));
                        
                    end
                    zm = AARpI\ALAWARD;
                    node{m}.W = AApR\ALAWD;
                    %             node{m}.W
                    for n=1:Nlen
                        Lnm = node{node{m}.N(n)}.L(:,node{node{m}.N(n)}.N==m);
                        Anm = flipud(node{node{m}.N(n)}.Amn{node{node{m}.N(n)}.N==m});
                        Wn = node{node{m}.N(n)}.W;
                        Amn = node{m}.Amn{n};
                        dm = zeros(Nlen,1);
                        dm(m) = node{m}.d;
                        node{m}.L(:,n) = Lnm-Anm*Wn-Amn/R*(dm+zm);
                    end
                end
            end
            Wsave(l,:) = [node{1}.W(1),node{2}.W(2),node{3}.W(3)];
        end
    end
end



%% Optimal weights
d = [node{1}.d;node{2}.d;node{3}.d];
Wopt = (R\d)/((d'/R)*d);

%% Calculate BF output
% Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
% mySpectrogram(Y);
% y = myOverlapAdd(Y);
% figure; plot(y);

%% MSE(W-Wopt)
for l=1:Ltmp
    W = squeeze(Wsave(l,:)).';
    MSEWWopt(l) = mean((W-Wopt).*conj(W-Wopt));
    varW(l) = var(W);
end
figure; semilogy(MSEWWopt); title('MSE(W-Wopt)'); grid on; 
varWopt = var(Wopt)
varWopt_varW = varWopt/var(W)
SNR = varWopt./var(W-Wopt)
SNRdb = 10*log10(SNR)
% figure; semilogy(varW); title('varW'); grid on; 
myNormalize(abs(W))
myNormalize(abs(Wopt))



%% MSE between W and Wopt
% WWoptMSE = zeros(Ltmp,1);
% for l=1:Ltmp
%     Wtmp = squeeze(Wsave(:,l,:));
%     WWoptMSE(l) = mean(mean((Wtmp-Wopt).*conj(Wtmp-Wopt)));
% %     WWoptMSE(l) = mean(mean((squeeze(Wsave(:,l,:))-Wopt).^2));
% end
% figure; plot((WWoptMSE)); grid on;

%% Let's have a look at W and Wopt
% figure; imagesc(abs(Wopt))
% figure; imagesc(abs(squeeze(Wsave(:,Ltmp,:))))

%% Print setup for records
% fprintf('Nsrcs = %d, K = %d, tls = %d, M = %d, spSize = %d, bin = %d, Ltmp = %d\n\n',Nsrcs,K,tls,M,spSize,bin,Ltmp);

