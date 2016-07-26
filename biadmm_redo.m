%% Start again again biadmm mvdr beamformer. 25/7/16

close all; clear;

% Import target audio
Nsrcs = 2; % Nsrcs = number of sources
s = cell(Nsrcs,1);

AudioFileNames = {'422-122949-0013.flac';'2078-142845-0005.flac';'2902-9006-0015.flac';'1272-128104-0004.flac';'422-122949-0014.flac';'2078-142845-0025.flac';'1919-142785-0007.flac';'422-122949-0009.flac';'174-168635-0018.flac';...
                  '2902-9006-0017.flac';'251-136532-0004.flac';'2078-142845-0004.flac';'2035-152373-0005.flac';'2902-9006-0001.flac';'1993-147964-0010.flac';'1673-143396-0004.flac';'2902-9008-0002.flac';'422-122949-0000.flac';...
                  '2078-142845-0039.flac';'1919-142785-0008.flac';'2412-153948-0004.flac';'2078-142845-0007.flac';'2078-142845-0043.flac';'1988-148538-0006.flac';'1919-142785-0005.flac';'174-84280-0004.flac'; '422-122949-0019.flac';'1993-147149-0006.flac'};
for ns=1:Nsrcs
    s{ns} = audioread(strcat('/audio/',AudioFileNames{ns})); 
end
fs = 16e3;
s{2} = 0.001*randn(length(s{2}),1); % For replacing the source with noise

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
Mloc = (rand(M,spcDim)*diag(space)).'; % Mloc = matrix containing 3d sensor locations
sloc = ((rand(Nsrcs,spcDim)*diag(space))).';%+[0,0,2;0,0,2]).'; % sloc = matrix containing 3d source locations

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

%% Create ATFs and observations using full length fft
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
    x = ifft(Xfft) + 0*randn(tl,1);
    xsave(:,m) = x;
    xPadded = [zeros((K-1)/2,1);x(1:end-1);zeros((K-1)/2,1)];
    XTmp = stft(xPadded,K);
    X(:,:,m) = XTmp(2:(K-1)/2,:);
end

%% Calculate covariances over all time
R = cell(Khalf,1);
RSum = zeros(M);
for k=1:Khalf
    for l=1:L
        if l==1
            R{k} = zeros(M);
        end
        XTmp = squeeze(X(k,l,:));
        R{k} = R{k} + (1/L)*(XTmp*XTmp');
    end
    RSum = RSum + R{k}; % RSum is the covariance summed over time and frequency bins
end
figure; imagesc(abs(RSum)); title('RSum');

%% Find sensor to sensor distances
sensd = zeros(M,M); % sensd = sensor to sensor distance
for m=1:M
    for mm=m:M
        sensd(m,mm) = norm(Mloc(:,m)-Mloc(:,mm));
    end
end
sensd = sensd + triu(sensd,1).'; % Convert from upper triangular to symmetric

%% Find neighbors
Nneighbors = zeros(M,1);
for m=1:M
    node{m}.N = [find(sensd(:,m)<2*spSize) ];
    node{m}.Nlen = length(node{m}.N);
    Nneighbors(m) = node{m}.Nlen;    
end
fprintf('The minimum number of neighbors was %d. \nThe mean number of neighbors was %d. \n\n',min(Nneighbors),mean(Nneighbors));

%% Initialization
fdomShort = (fs/(K-1)) * (1:((K-1)/2)-1);
for m=1:M
    % Initialize Lambdas and weights
    node{m}.L = ones(Khalf,2,node{m}.Nlen); % These are for node m's real neighbors including itself
    node{m}.W = ones(Khalf,node{m}.Nlen); % Note that this includes node m itself, node m's real neighbors, and node m's virtual neighbor (m+M)
    
%     % Initialize Amn for all nodes
%     Amn = cell(node{m}.Nlen,1); 
%     for n = 1:node{m}.Nlen
%         if node{m}.N(n)==m
%             Amn{n} = zeros(2,node{m}.Nlen);
%         else
%             Amn{n} = [(node{m}.N==m).';-(node{m}.N==node{m}.N(n)).'];
%         end
%     end
%     node{m}.Amn = Amn;
   
    % Save look direction d for node m
    node{m}.d = exp(-1i*2*pi*fdomShort.'*ssd(1,m)/c) / (4*pi*ssd(1,m));
    
%     % Initialize R
%     for k=1:Khalf
%         node{m}.R{k} = 1e-6*eye(node{m}.Nlen); % Note that R now includes a dimension for the virtual node, which seems meaningless at first glance
%     end
end

% Initialize regularization parameter alpha
alpha = 1;

% Initialize output Y
Y = zeros(Khalf,L);



%% Iterative update
% Step through windows
wcnt=1;
Ltmp = 10   ;
ITER2 = 2;
Wsave = zeros(Khalf,Ltmp,M);
for l=1:Ltmp
    
    % Step through real (and virtual) nodes for primal/dual update
    for iter2 = 1:ITER2
        for m=1:M
            [l,m]
            Nlen = node{m}.Nlen;
%             Rm = zeros(Khalf,Nlen,Nlen);
            Xm = zeros(Khalf,Nlen);
            Nm = node{m}.N;
            dm = node{m}.d;
            Wm = node{m}.W;
            
            % iter sets the number of iterative updates per window
            for iter=1
                
                % Get observations
                %             node{m}.X = squeeze(X(:,l,node{m}.N));
                %             Xm = squeeze(X(:,l,node{m}.N)); % virtual node has real node's observations
                
                %             % Assign covariance from actual covariance calculated above
                %             if l==1 && iter==1
                %                 for k=1:Khalf
                %                     RTemp = zeros(node{m}.Nlen-1);
                %                     for n=1:node{m}.Nlen-1
                %                         RTemp(:,n) = R{k}(m,node{m}.N(1:end-1)).';
                %                         RTemp(n,:) = R{k}(m,node{m}.N(1:end-1));
                %                     end
                %                     node{m}.R{k} = RTemp;
                %                 end
                %             end
                
                % Calculate local covariance
                %             % Calculation only needs to happen when the node changes
                %             if iter==1
                %                 for k=1:Khalf
                %                     Rm(k,:,:) = (1-gamma)*node{m}.R{k} + gamma*(Xm(k,:).'*Xm(k,:));
                %                 end
                %             end
                
                %             if iter==1
                %                 for k=1:Khalf
                % %                     Rm(k,:,:) = Xm(k,:).'*Xm(k,:) + 1e-6*eye(Nlen);
                %                     Rm(k,:,:) = RSum;
                %                 end
                %             end
                
                
                %             % Initialize temp vars
                %             AAsum = zeros(Nlen,Nlen);
                %             ALAWsum = zeros(Nlen,1);
                %             AARsum = zeros(Nlen,Nlen);
                %             ALAWARdsum = zeros(Nlen,1);
                %             WRealNew = zeros(Khalf,Nlen);
                %             WVirtNew = zeros(Khalf,1);
                %             LambdaRealNew = zeros(Khalf,2,Nlen);  % WHY DO THESE START AS ONES? because they are similar to weights? which means to leave the system unchanged they need to be ones, not zeros ?
                %             bk = zeros(Khalf,1);
                %
                %             % Calculate primal update for real node
                %             for k=1:Khalf
                %                 dTmp = (Nm(1:end)==m) * dm(k);
                %
                %                 for n=1:Nlen
                %                     % for W update
                %                     AmnTmp = +node{m}.Amn{n};
                %                     AAsum = AAsum+(AmnTmp.'*AmnTmp);
                %                     LambdanmTmp = node{Nm(n)}.L(k,:,(node{Nm(n)}.N==m)).'; % This is just the real-real lambdas
                %                     AnmTmp = flipud(node{Nm(n)}.Amn{node{Nm(n)}.N==m});
                %                     WnTmp = node{Nm(n)}.W(k,:).';
                %                     ALAWsum = ALAWsum + (AmnTmp.'*(LambdanmTmp-AnmTmp*WnTmp));
                %
                %                     % For Lambda,zk update
                %                     AmnRinv = AmnTmp/squeeze(Rm(k,:,:));
                %                     AARsum = AARsum + (AmnTmp.'*AmnRinv);
                %                     ALAWARdsum = ALAWARdsum + (AmnTmp.'*(LambdanmTmp-AnmTmp*WnTmp-AmnRinv*dTmp));
                %                 end
                %                 WRealNew(k,1:Nlen) = (AAsum+squeeze(Rm(k,:,:)))\(ALAWsum+dTmp);
                %                 zk = (AARsum+eye(Nlen))\ALAWARdsum;
                %
                %                 % Lambda update requires twice around the neighboring nodes:
                %                 % once for finding the summations in zk, and once for Lambda
                %                 % itself.
                %                 for n=1:node{m}.Nlen
                %                     LambdanmTmp = node{Nm(n)}.L(k,:,(node{Nm(n)}.N==m)).';
                %                     AnmTmp = flipud(node{Nm(n)}.Amn{node{Nm(n)}.N==m});
                %                     WnTmp = node{Nm(n)}.W(k,:).';
                %                     AmnTmp = node{m}.Amn{n};
                %                     LambdaRealNew(k,:,n) = LambdanmTmp-AnmTmp*WnTmp-AmnTmp/squeeze(Rm(k,:,:))*(dTmp+zk);
                %                 end
                %
                %                 % Primal update for virtual node m+M
                %                 WnTmp = Wm(k,:); % Comes from m because m is n for virtual node m+M
                %                 AnmTmp = node{m}.Amn{end};
                %                 LambdanmTmp = node{m}.L(k);
                %             end
                %
                %             % Assign new W and L values to node m and m+M
                %             node{m}.W = WRealNew;
                %             node{m}.L = LambdaRealNew;
                
                % Calculate primal update for real node
                Wnew = zeros(Khalf,Nlen);
                Lnew = zeros(Khalf,2,Nlen);
                for k=1:Khalf
                    AAplusR = zeros(Nlen);
                    ALAWD = zeros(Nlen,1);
                    AARI = zeros(Nlen);
                    ALAWARD = zeros(Nlen,1);
                    for n=1:Nlen
                        if m==Nm(n)
                            Amn = zeros(2,Nlen);
                            Anm = zeros(2,node{node{m}.N(n)}.Nlen);
                        else
                            Amn = [Nm==m,-(Nm==n)].';
                            Anm = [-(node{node{m}.N(n)}.N==m),(node{node{m}.N(n)}.N==n)].';
                        end
                        % For primal
                        AAplusR = AAplusR + (Amn.'*Amn+R{k});
                        Lnm = node{node{m}.N(n)}.L(k,:,node{node{m}.N(n)}.N==m).';
                        Wn = node{node{m}.N(n)}.W(k,:).';
                        d = zeros(Nlen,1);
                        d(Nm==m) = dm(k);
                        ALAWD = ALAWD + (Amn.'*(Lnm-Anm*Wn)+d);
                        
                        % For dual
                        AmnRinv = Amn/R{k};
                        AARI = AARI + (Amn.'*AmnRinv+eye(Nlen)) ;
                        ALAWARD = ALAWARD + (Amn.'*(Lnm-Anm*Wn-AmnRinv*d));
                    end
                    Wnew(k,:) = AAplusR\ALAWD;
                    z = AARI\ALAWARD;
                    for n=1:Nlen
                        if m==Nm(n)
                            Amn = zeros(2,Nlen);
                            Anm = zeros(2,node{node{m}.N(n)}.Nlen);
                        else
                            Amn = [Nm==m,-(Nm==n)].';
                            Anm = [-(node{node{m}.N(n)}.N==m),(node{node{m}.N(n)}.N==n)].';
                        end
                        Lnm = node{node{m}.N(n)}.L(k,:,node{node{m}.N(n)}.N==m).';
                        Wn = node{node{m}.N(n)}.W(k,:).';
                        AmnRinv = Amn/RSum;
                        d = zeros(Nlen,1);
                        d(Nm==m) = dm(k);
                        Lnew(k,:,n) = Lnm-Anm*Wn-AmnRinv*(d+z);
                    end
                end
                node{m}.W = Wnew;
                node{m}.L = Lnew;
            end
            W(:,m) = node{m}.W(:,node{m}.N==m);
        end
        Wsave(:,wcnt,:) = W;
        wcnt=wcnt+1;
    end
    for k=1:Khalf
%         node{m}.R{k} = squeeze(Rm(k,:,:));
        
        Y(k,l) = (W(k,:))*squeeze(X(k,l,:));
    end
%     Wsave(:,l,:) = W;
    
    % Generate output using updated primal and dual weights
%     for m=1:M
% %         Y(:,l) = Y(:,l) + (1/M)*conj(node{m}.W(:,node{m}.N==m)).*node{m}.X(:,node{m}.N==m);
%         
%         Wsave(:,l,m) = (1/M)*node{m}.W(:,node{m}.N==m);
%         Wnorm(l,m) = norm(Wsave);
%         Y(:,l) = Y(:,l) + (1/M)*sum(node{m}.W(:,1:end-1).*node{m}.X(:,1:end-1),2); % -1 because the virtual node's weight and observation should not contribute to the output
%     end     
    
%     figure; imagesc(Wsave);

%     Wnorm(l) = norm(Wsave);
end


%% Calculate BF output
Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
mySpectrogram(Y);
y = myOverlapAdd(Y);
figure; plot(y); title('bf output');



%%
% figure; plot(Wnorm)

%% A1) Frequency domain Frost optimum weights
d = zeros(Khalf,M);
for m=1:M
    d(:,m) = node{m}.d;
end

Wopt = zeros(Khalf,M);
dk = zeros(Khalf,M);
for k=1:Khalf
    dk = d(k,:).';
%     Wopt(k,:) = RSum\dk/(dk'/RSum*dk); 
    Wopt(k,:) = R{k}\dk/(dk'/R{k}*dk);
end

% Find output using optimum weights
Yopt = zeros(Khalf,L);
for l=1:L
    Xtmp = squeeze(X(:,l,:));
    Yopt(:,l) = sum(conj(Wopt).*Xtmp,2);
end
Yopt = [zeros(1,L);Yopt;zeros(2,L);conj(flipud(Yopt))];
yopt = myOverlapAdd(Yopt);


%% 
WWoptMSE = zeros(wcnt,1);
wcnt=wcnt-1;
for l=1:wcnt
    WWoptMSE(l) = mean(mean((squeeze(Wsave(:,l,:))-Wopt).^2));
end
figure; plot(abs(WWoptMSE)); grid on; title('WWoptMSE')


%% 
WWoptNORM = zeros(wcnt,1);
for l=1:wcnt
    WWoptNORM(l) = norm(squeeze(Wsave(:,l,:))-Wopt);
end
figure; plot(abs(WWoptNORM)); grid on; title('WWoptNORM');

%% Visualize weights
figure; subplot(2,1,1); imagesc(abs(Wopt)); subplot(2,1,2); imagesc(abs(squeeze(Wsave(:,wcnt,:))));
% figure; subplot(2,1,1); imagesc(angle(Wopt)); subplot(2,1,2); imagesc(angle(squeeze(Wsave(:,l,:))));

