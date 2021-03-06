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

AudioFileNames = {'422-122949-0013.flac';
                  '2078-142845-0005.flac';
                  '2902-9006-0015.flac';
                  '1272-128104-0004.flac';
                  '422-122949-0014.flac';
                  '2078-142845-0025.flac';
                  '1919-142785-0007.flac';
                  '422-122949-0009.flac';
                  '174-168635-0018.flac';
                  '2902-9006-0017.flac';
                  '251-136532-0004.flac';
                  '2078-142845-0004.flac';
                  '2035-152373-0005.flac';
                  '2902-9006-0001.flac';
                  '1993-147964-0010.flac';
                  '1673-143396-0004.flac';
                  '2902-9008-0002.flac';
                  '422-122949-0000.flac';
                  '2078-142845-0039.flac';
                  '1919-142785-0008.flac';
                  '2412-153948-0004.flac';
                  '2078-142845-0007.flac';
                  '2078-142845-0043.flac';
                  '1988-148538-0006.flac';
                  '1919-142785-0005.flac';
                  '174-84280-0004.flac';
                  '422-122949-0019.flac';
                  '1993-147149-0006.flac'};

for ns=1:Nsrcs
    s{ns} = audioread(strcat('/audio/',AudioFileNames{ns})); 
end
fs = 16e3;
s{2} = 0.001*randn(length(s{2}),1);

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
for m=1:2*M
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
    node{m+M}.loc = 0; % Virtual nodes don't have a real location
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
    x = ifft(Xfft) + 0.000001*randn(tl,1);
    xsave(:,m) = x;
    xPadded = [zeros((K-1)/2,1);x(1:end-1);zeros((K-1)/2,1)];
    XTmp = stft(xPadded,K);
    X(:,:,m) = XTmp(2:(K-1)/2,:);
end

%% Calculate covariances over all time
% R = cell(Khalf,1);
% for k=1:Khalf
%     for l=1:L
%         if l==1
%             R{k} = zeros(M);
%         end
%         R{k} = R{k} + (1/L)*squeeze(X(k,l,:))*squeeze(X(k,l,:))';
%     end
% end

%% Calculate covariances over all time
R = cell(Khalf,1);
for k=1:Khalf
    for l=1:L
        if l==1
            R{k} = zeros(M);
        end
        XTmp = squeeze(X(k,l,:));
        R{k} = R{k} + (1/L)*(XTmp*XTmp');
    end
end

%% 
RR = xsave'*xsave;
RRcondition = rcond(RR)
figure; imagesc(abs(RR));

RSum = zeros(M);
for k=1:Khalf
    RSum = RSum + R{k}; % RSum is the covariance summed over time and frequency bins
end
figure; imagesc(abs(RSum));

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
    node{m}.N = [find(sensd(:,m)<2*spSize) ; m+M];
    node{m}.Nlen = length(node{m}.N);
    Nneighbors(m) = node{m}.Nlen;    
    node{m+M}.N = [m+M;m];
    node{m+M}.Nlen = 2;
end
fprintf('The minimum number of neighbors was %d. \nThe mean number of neighbors was %d. \n\n',min(Nneighbors),mean(Nneighbors));

%% Initialization
fdomShort = (fs/(K-1)) * (1:((K-1)/2)-1);
for m=1:M
    % Initialize Lambdas and weights
%     node{m}.L = zeros(Khalf,2,node{m}.Nlen-1); % Note that the virtual node has L = 1x1, so the second element is redundant. I have kept it in for ease of coding
    node{m}.L = cell(2,1); % I need two different sizes for a real node's lambdas: one for real-real, one for real-virtual
    node{m}.L{1} = ones(Khalf,2,node{m}.Nlen-1); % These are for node m's real neighbors including itself
    node{m}.L{2} = ones(Khalf,1); % This is for node m's virtual node
%     node{m+M}.L = zeros(Khalf,1,1); 
    node{m+M}.L{1} = ones(Khalf,1,1); % This is for virtual node (m+M)'s 1x1 lambda
    node{m}.W = ones(Khalf,node{m}.Nlen); % Note that this includes node m itself, node m's real neighbors, and node m's virtual neighbor (m+M)
    node{m+M}.W = ones(Khalf,1); % There is only one weight stored on a virtual node
    
    % Initialize Amn for all nodes
    Amn = cell(node{m}.Nlen,1); 
    for n = 1:node{m}.Nlen-1 % -1 because I have to set up the virtual node differently
        Amn{n} = [(node{m}.N==m).';-(node{m}.N==node{m}.N(n)).'];
    end
%     Amn{n+1} = [zeros(1,node{m}.Nlen);-(node{m}.N==node{m}.N(n+1)).'];
    Amn{n+1} = (node{m}.N==m).';
    node{m}.Amn = Amn;
%     node{m+M}.Amn{1} = [1;0]; % This is a special consistency matrix A
%     specifically for virtual nodes. % Don't need this now as all virtual
%     nodes have Amn = 1 i.e. (1x1).
    
    % Save look direction d for node m
    node{m}.d = exp(-1i*2*pi*fdomShort.'*ssd(1,m)/c) / (4*pi*ssd(1,m));
    
    % Initialize R
    for k=1:Khalf
%         node{m}.R{k} = 1e-12*eye(node{m}.Nlen-1); 
        node{m}.R{k} = 1e-6*eye(node{m}.Nlen); % Note that R now includes a dimension for the virtual node, which seems meaningless at first glance

    end
    
end

% Initialize regularization parameter alpha
alpha = 1;

% Initialize output Y
Y = zeros(Khalf,L);

%% Iterative update
% Step through windows
gamma = 0.4; % No reason for 0.4, this is just what seems to work
Ltmp = 10;
for l=1:Ltmp
    
    % Step through real (and virtual) nodes for primal/dual update
    for m=1:M
        [l,m]
        Nlen = node{m}.Nlen;
%         Rm = zeros(Khalf,Nlen-1,Nlen-1);
        Rm = zeros(Khalf,Nlen,Nlen);
%         Xm = zeros(Khalf,Nlen-1);
        Xm = zeros(Khalf,Nlen);
        Nm = node{m}.N;
        dm = node{m}.d;
        Wm = node{m}.W;
        
        % iter sets the number of iterative updates per window
        for iter=1
            
            % Get observations
%             node{m}.X = squeeze(X(:,l,node{m}.N(1:end-1))); % -1 to exclude the virtual node
            node{m}.X = squeeze(X(:,l,[node{m}.N(1:end-1);m])); % DODGINESS ALERT --> Here I have concatenated the real node observations into the virtual node's dimension, as the virtual node has no observations of its own.
%             Xm = squeeze(X(:,l,node{m}.N(1:end-1))); % -1 to exclude the virtual node
            Xm = squeeze(X(:,l,[node{m}.N(1:end-1);m])); % virtual node has real node's observations
            
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
            % Initialize
%             if l==1
%                 for k=1:Khalf
% %                     node{m}.R{k} = 1e-12*eye(Nlen-1);
%                     Rm(k,:,:) = 1e-12*eye(Nlen-1);
%                 end
%             end
            % Calculation only needs to happen when the node changes
            if iter==1
                for k=1:Khalf
%                     node{m}.R{k} = (1-gamma)*node{m}.R{k} + gamma*(node{m}.X(k,:).'*node{m}.X(k,:));
%                     Rm(k,:,:) = (1-gamma)*node{m}.R{k} + gamma*(Xm(k,:).'*Xm(k,:));
                    Rm(k,:,:) = RSum;
                end
            end
            
            % Initialize temp vars
%             AAsum = zeros(Nlen-1,Nlen-1);
            AAsum = zeros(Nlen,Nlen);
%             ALAWsum = zeros(Nlen-1,1);
            ALAWsum = zeros(Nlen,1);
%             AARsum = zeros(Nlen-1,Nlen-1);
            AARsum = zeros(Nlen,Nlen);
%             ALAWARdsum = zeros(Nlen-1,1);
            ALAWARdsum = zeros(Nlen,1);
            WRealNew = ones(Khalf,Nlen);
            WVirtNew = ones(Khalf,1);
            LambdaRealNew{1} = ones(Khalf,2,Nlen);  % WHY DO THESE START AS ONES? because they are similar to weights? which means to leave the system unchanged they need to be ones, not zeros ? 
            LambdaRealNew{2} = ones(Khalf,Nlen);    % For the virtual node
            LambdaVirtNew = ones(Khalf,1);
%             LambdaNew = ones(Khalf,2,M); 
            bk = zeros(Khalf,1);
            
            % Calculate primal update for real node
            for k=1:Khalf
                
%                 dTmp = (Nm(1:end-1)==m) * dm(k);
                dTmp = (Nm(1:end)==m) * dm(k);
                
%                 for n=1:node{m}.Nlen-1 % -1 because I want to exclude the virtual node for now
                for n=1:Nlen
                    
                    % for W update
%                     AmnTmp = node{m}.Amn{n}(:,1:node{m}.Nlen-1);
                    AmnTmp = +node{m}.Amn{n};
                    AAsum = AAsum+(AmnTmp.'*AmnTmp);
                    if n<Nlen
                        LambdanmTmp = node{Nm(n)}.L{1}(k,:,(node{Nm(n)}.N==m)).'; % This is just the real-real lambdas
                        AnmTmp = flipud(node{Nm(n)}.Amn{node{Nm(n)}.N==m});
                        WnTmp = node{Nm(n)}.W(k,:).';
                        ALAWsum = ALAWsum + (AmnTmp.'*(LambdanmTmp-AnmTmp*WnTmp));
                    else 
                        LambdanmTmp = node{Nm(n)}.L{1}(k);
%                         AnmTmp = node{Nm(n)}.Amn{node{Nm(n)}.N==m};
                        AnmTmp = 1;
                        WnTmp = node{Nm(n)}.W(k);
                        ALAWsum = ALAWsum + (AmnTmp.'*(LambdanmTmp-AnmTmp*WnTmp));
                    end
                    


                    % For Lambda,zk update
                    AmnRinv = AmnTmp/squeeze(Rm(k,:,:));
                    AARsum = AARsum + (AmnTmp.'*AmnRinv); 
                    ALAWARdsum = ALAWARdsum + (AmnTmp.'*(LambdanmTmp-AnmTmp*WnTmp-AmnRinv*dTmp));  
                end

                
%                 WRealNew(k,1:Nlen-1) = (AAsum+squeeze(Rm(k,:,:)))\(ALAWsum+dTmp);
                WRealNew(k,1:Nlen) = (AAsum+squeeze(Rm(k,:,:)))\(ALAWsum+dTmp);
%                 zk = (AARsum+eye(Nlen-1))\ALAWARdsum;
                zk = (AARsum+eye(Nlen))\ALAWARdsum;

                % Lambda update requires twice around the neighboring nodes:
                % once for finding the summations in zk, and once for Lambda
                % itself.
                for n=1:node{m}.Nlen
                    
                    if n<Nlen
%                     LambdanmTmp = node{Nm(n)}.L(k,:,(node{Nm(n)}.N==m)).';
                        LambdanmTmp = node{Nm(n)}.L{1}(k,:,(node{Nm(n)}.N==m)).';
                        AnmTmp = flipud(node{Nm(n)}.Amn{node{Nm(n)}.N==m});
                        WnTmp = node{Nm(n)}.W(k,:).';
%                       AmnTmp = node{m}.Amn{n}(:,1:end-1);      
                        AmnTmp = node{m}.Amn{n};  
%                       LambdaNew(k,:,n) = LambdanmTmp-AnmTmp*WnTmp-AmnTmp/squeeze(Rm(k,:,:))*(dTmp+zk);
                        LambdaRealNew{1}(k,:,n) = LambdanmTmp-AnmTmp*WnTmp-AmnTmp/squeeze(Rm(k,:,:))*(dTmp+zk);
                    else
                        LambdanmTmp = node{Nm(n)}.L{1}(k);
                        AnmTmp = 1;
                        WnTmp = node{Nm(n)}.W(k,:).';
%                       AmnTmp = node{m}.Amn{n}(:,1:end-1);      
                        AmnTmp = node{m}.Amn{n};  
%                       LambdaNew(k,:,n) = LambdanmTmp-AnmTmp*WnTmp-AmnTmp/squeeze(Rm(k,:,:))*(dTmp+zk);
                        LambdaRealNew{2}(k,n) = LambdanmTmp-AnmTmp*WnTmp-AmnTmp/squeeze(Rm(k,:,:))*(dTmp+zk);
                    end
                end 
                
                % Primal update for virtual node m+M
                WnTmp = Wm(k,:); % Comes from m because m is n for virtual node m+M
%                 AnmTmp = flipud(node{m}.Amn{end}); % Its ok to use end because m+M is always last 
                AnmTmp = node{m}.Amn{end};
%                 AmnTmp = node{m+M}.Amn{1};
%                 LambdanmTmp = node{m}.L(k,:,end);
                LambdanmTmp = node{m}.L{2}(k);
                bk = 2*WnTmp*AnmTmp.'*1-1*LambdanmTmp.'; % 1 because Amn is always 1 when m is a virtual node
                WVirtNew(k) = (-bk+sign(bk)*min(abs(bk),alpha))/(1);

                % Dual update for virtual node m+M
                LambdamnTmp = node{m+M}.L{1}(k);
                LambdaVirtNew(k) = 2*LambdanmTmp*LambdamnTmp.'-WnTmp*AnmTmp.';
                if LambdaVirtNew(k) > alpha
                    LambdaVirtNew(k) = alpha; 
                elseif LambdaVirtNew(k) < -alpha
                    LambdaVirtNew(k) = -alpha; 
                end
%                 if LambdaVirtNew(k,2) > alpha
%                     LambdaVirtNew(k,2) = alpha; 
%                 elseif LambdaVirtNew(k,2) < -alpha
%                     LambdaVirtNew(k,2) = -alpha; 
%                 end                  
            end   

            % Assign new W and L values to node m and m+M
            node{m}.W = WRealNew;
            node{m+M}.W = WVirtNew;
            node{m}.L{1} = LambdaRealNew{1};
            node{m}.L{2} = LambdaRealNew{2};
            node{m+M}.L{1} = LambdaVirtNew;
        end
    end  
    for k=1:Khalf
        node{m}.R{k} = squeeze(Rm(k,:,:));
    end
    
    % Generate output using updated primal and dual weights
    for m=1:M
        Wsave(:,l,:) = node{m}.W(node{m}.N==m);
        Y(:,l) = Y(:,l) + (1/M)*sum(conj(node{m}.W(:,1:end-1)).*node{m}.X(:,1:end-1),2); % -1 because the virtual node's weight and observation should not contribute to the output
    end       
end


%% Calculate BF output
Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
mySpectrogram(Y);
y = myOverlapAdd(Y);
figure; plot(y);

%% A1) Frequency domain Frost optimum weights
d = zeros(Khalf,M);
for m=1:M
    d(:,m) = node{m}.d;
end

Wopt = zeros(Khalf,M);
dk = zeros(Khalf,M);
for k=1:Khalf
    dk = d(k,:).';
    Wopt(k,:) = RSum\dk/(dk'/RSum*dk); 
end

% Find output using optimum weights
Yopt = zeros(Khalf,L);
for l=1:L
    Xtmp = squeeze(X(:,l,:));
    Yopt(:,l) = sum(conj(Wopt).*Xtmp,2);
end
Yopt = [zeros(1,L);Yopt;zeros(2,L);conj(flipud(Yopt))];
yopt = myOverlapAdd(Yopt);

%% MSE between W and Wopt
for l=1:Ltmp
    WWoptMSE(l) = mean(mean(conj(Wopt-squeeze(Wsave(:,l,:))).*(Wopt-squeeze(Wsave(:,l,:)))));
end
figure; grid on; plot(WWoptMSE);

%% Let's have a look a W and Wopt
figure; imagesc(abs(Wopt))
figure; imagesc(abs(squeeze(Wsave(:,Ltmp,:))))













