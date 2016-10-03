
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

% pdmm_ls.m is a working pdmm for least squares with guaranteed convergence
% to the optimum.

% biadmm1class10.m is working pdmm - note that it cannot guarantee
% convergece to Wopt, but it sounds great and the weights are only a factor
% out

% biadmm1class11.m has the virtual nodes added, and it runs. It looks like
% it could be working, though it is still fully connected to can be hard to
% say that the sparsity is really working. It looks like it pushes the
% smaller weights down in comparison to the largest weight

% biadmm1class12.m, make the network consist of local neighbors only. 12 is
% working and it looks like sparsity is working. I can run 100 sensors (with
% min=4,mean=18 neighbors) for one frequency bin very quickly. 200 sensors
% freezes up the computer, so uses too much memory(?).

% biadmm1class12a.m: look at forcing sparsity to use separated mics using
% two sources and four sensors, two sensors close to each source at the
% same position, or very close, or the same distance, and use sparsity to
% cause only two of the sensors to be active. Ideally the regularization
% would choose one sensor close to each source.

% biadmm1class13.m: Diffuse noise source added

% biadmm1class14.m - can I start from one node and grow? Using a pragmatic
% approach i.e. a condition on inclusion, not part of the objective.

% biadmm1class15.m - copy of biadmm1class10.m because it's pre-sparsity.
% Implement a simple case of starting with a single node and growing to
% meet SNR target.

% biadmm1class16.m - similar to 15, I took 15 home, changed it, saved it as
% 16. I'm scared of syncing my home computer to git.

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

% Truncate to desired length, ensuring that the length is a multiple of
% the window length.
K = 2^9+1; % K = window length in samples, and the number of frequency bins
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
M = 10^2; % M = number of sensors

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
% sloc =   spSize*[0.1;
%                 0.1;
%                 0.1];
% xloc(:,1:3) = spSize*[0.11,0.3,0.91;
%                       0.11,0.4,0.91;
%                       0.11,0.5,0.91];
% xloc = spSize*[0.2,0.1,0.1;
%                 0.1,0.2,0.1;
%                 0.1,0.1,0.2];
sloc = (rand(Nsrcs,spcDim-1)*diag(space(1:end-1))).';
sloc = [sloc;0.5*ones(1,Nsrcs)];
xloc = spSize*[combvec(linspace(0,1,sqrt(M)), linspace(0,1,sqrt(M)));repmat(1,1,M)]+(0.03*randn(3,M));

% Set location for each node
for m=1:M
    node{m}.loc = xloc(:,m);
end

% Calculate distances
ssd = myGetDist(xloc,sloc);
[dontcare,nearestSensor] = min(ssd(1,:));

%% Display layout
% myDisplayLayout(xloc,sloc);

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
sensd = myFindSensDist(xloc);

%% Find neighbors, everyone within 0.5*spSize (Aiming for >5 neighbors each, but depends on the number of nodes as well)
% Nneighbors = zeros(M,1);
% for m=1:M
%     node{m}.N = [find(sensd(:,m)<2*spSize) ];
%     node{m}.Nlen = length(node{m}.N);
%     Nneighbors(m) = node{m}.Nlen;    
% end
% fprintf('The minimum number of neighbors was %d. \nThe mean number of neighbors was %d. \n\n',min(Nneighbors),mean(Nneighbors));

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
    R(k,:,:) = squeeze(R(k,:,:)) + 1e-9*eye(M);
end

%% Initialization
for m=1:M
%     % Initialize Lambdas and weights
%     node{m}.L = zeros(Khalf,2,node{m}.Nlen); % These are for node m's real neighbors including itself
%     node{m}.W = zeros(Khalf,node{m}.Nlen); % Note that this includes node m itself, node m's real neighbors, and node m's virtual neighbor (m+M)
%     node{m}.Lnew = zeros(Khalf,2,node{m}.Nlen);
%     node{m}.Wnew = zeros(Khalf,node{m}.Nlen);
    
%     % Initialize Amn for all nodes
%     Amn = cell(node{m}.Nlen,1);
%     for n = m:node{m}.Nlen % Note that this loop starts from m
%         if m==node{m}.N(n)
%             node{m}.Amn{n} = zeros(2,node{m}.Nlen);
%         else
%             node{m}.Amn{n} = double([(node{m}.N==m).';(node{m}.N==node{m}.N(n)).']);
%             node{node{m}.N(n)}.Amn{m} = -node{m}.Amn{n};
%         end
%     end


 
    % Save look direction d for node m
    node{m}.d = exp(-1i*2*pi*fdomShort.'*ssd(1,m)/c) / (4*pi*ssd(1,m));
    
    % Save covariance for nearestSensor
%     node{m}.R = R(:,[node{m}.N(1:end-1)],[node{m}.N(1:end-1)]);
    

end

% Initialize Lambdas and weights
node{nearestSensor}.L = zeros(Khalf,2); % These are for node m's real neighbors including itself
node{nearestSensor}.W = zeros(Khalf,1); % Note that this includes node m itself, node m's real neighbors, and node m's virtual neighbor (m+M)
node{nearestSensor}.Lnew = zeros(Khalf,2);
node{nearestSensor}.Wnew = zeros(Khalf,1);

node{nearestSensor}.N = [nearestSensor];
node{nearestSensor}.Amn{1} = [1;0];
node{nearestSensor}.Nlen = 1;
node{nearestSensor}.R = R(:,nearestSensor,nearestSensor);

% Initialize output Y
Y = zeros(Khalf,L);

%% Find input SIR

bin = 53;

[mindist,mindisti] = min(ssd(1,:));                 

% % Get stft of the sources
% Target = stft(s{1}(1:tl),K);
% Target = Target(2:(K-1)/2,:);
% Target = Target(bin,:);
% Interferer = stft(s{2}(1:tl),K);
% Interferer = Interferer(2:(K-1)/2,:);
% Interferer = Interferer(bin,:);
dT = exp(-1i*2*pi*fdomShort(bin)*ssd(1,mindisti)/c)/(4*pi*ssd(1,mindisti));
dI = exp(-1i*2*pi*fdomShort(bin)*ssd(2,mindisti)/c)/(4*pi*ssd(2,mindisti));

% iSIR = 10*log10((Target*dT*dT'*Target')/(Interferer*dI*dI'*Interferer'))
oSIRtf(1) = 10*log10((dT'*dT)/(dI'*dI));
% iSIRnotf = 10*log10((Target*Target')/(Interferer*Interferer'))

% % Increase size of network based on initial SIR
% % Choose size
% networkSize = 2;
% 
% % Find neighborhood of closest sensors
% % list sensors in terms of distance from nearestSensor
% sortsensd = sort(sensd(:,nearestSensor));
% [dontcare,dontcare,activeNodes] = intersect(sortsensd(1:networkSize),sensd(:,nearestSensor));
% 
% % take the closest networkSize sensors and set up network between them
% for m=activeNodes.'
%     node{m}.N = activeNodes;
%     node{m}.Nlen = length(activeNodes);
% end




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
% Ltmp = L; % For shorter run times

ITER1 = 20;
ITER2 = 1;

% % Initialize W to Wopt
% for m=1:M
%     node{m}.W = Wopt;     
% end

ftmp = zeros(ITER1,M);
rho = 1.3; % scaling for consensus
alpha = 1; % scaling for lambda consensus

% Check pdmm using ||Bx-b||^2
% B = randn(M);
% b = randn(M,1);
% Wopt = inv(B)*b;




for ani=2:5 % ani = active node iteration, i.e. each time through this loop has a different network size. If the network hasn't changed then there is no reason to recalculate the weights
    ani
    
    % list sensors in terms of distance from nearestSensor
    sortsensd = sort(sensd(:,nearestSensor));
    
    % Find the active nodes
    [dontcare,dontcare,activeNodes] = intersect(sortsensd(1:ani),sensd(:,nearestSensor));
    
    % Divide active nodes into neighborhoods
    Nneighbors = zeros(length(activeNodes),1);
    a=1;
    activeNodesLogical = zeros(M,1);
    for m=activeNodes.' 
        activeNodesLogical(m) = 1;
    end
    for m=activeNodes.'
        
        node{m}.N = [find((sensd(:,m)<0.5*spSize) & (activeNodesLogical)) ];
        node{m}.Nlen = length(node{m}.N);
        Nneighbors(a) = node{m}.Nlen;
        a=a+1;
    end
    fprintf('The minimum number of neighbors was %d. \nThe mean number of neighbors was %d. \n\n',min(Nneighbors),mean(Nneighbors));

    
%     % take the closest networkSize sensors and set up network between them
%     for m=activeNodes.'
%         node{m}.N = activeNodes;
%         node{m}.Nlen = length(activeNodes);
%     end
    % Initialize network
    for m=activeNodes.'
        mi=find(activeNodes==m); % mi = m index
        % Initialize Lambdas and weights
        node{m}.L = zeros(Khalf,2,node{m}.Nlen); % These are for node m's real neighbors including itself
        node{m}.W = zeros(Khalf,node{m}.Nlen); % Note that this includes node m itself, node m's real neighbors, and node m's virtual neighbor (m+M)
        node{m}.Lnew = zeros(Khalf,2,node{m}.Nlen);
        node{m}.Wnew = zeros(Khalf,node{m}.Nlen);
        
        % Initialize Amn for all nodes
        Amn = cell(node{m}.Nlen,1);
        for ni = mi:node{m}.Nlen % Note that this loop starts from m
            if m==node{m}.N(ni)
                node{m}.Amn{ni} = zeros(2,node{m}.Nlen);
            else
                node{m}.Amn{ni} = double([(node{m}.N==m).';(node{m}.N==node{m}.N(ni)).']);
                node{node{m}.N(ni)}.Amn{mi} = -node{m}.Amn{ni};
            end
        end
        
        
        
        % Save look direction d for node m
        node{m}.d = exp(-1i*2*pi*fdomShort.'*ssd(1,m)/c) / (4*pi*ssd(1,m));
        
        % Save covariance for each active node
        node{m}.R = R(:,[node{m}.N],[node{m}.N]);
        
        
    end
    
    for l=1
        for iter1=1:ITER1
            for k=1:Khalf
                %             for m=1:M
                for m=activeNodes.'
                    
                    mi=find(activeNodes==m); % mi = m index
                    for iter2=1:ITER2
                        Nlen = node{m}.Nlen;
                        AA = zeros(Nlen);
                        ALAW = zeros(Nlen,1);
                        dm = zeros(Nlen,1);
                        % dm = node{m}.d(k);
                        dm(node{m}.N(1:end-1)==m) = node{m}.d(k);
                        
                        % W update
                        for n=1:Nlen
                            Amn = node{m}.Amn{n};
                            AA = AA + (Amn.'*Amn);
%                             if Nlen == 1
%                                 Lnm = node{node{m}.N(n)}.L(k,:).';
%                             else
                                Lnm = node{node{m}.N(n)}.L(k,:,node{node{m}.N(n)}.N==m).';
%                             end
                            Anm = node{node{m}.N(n)}.Amn{node{node{m}.N(n)}.N==m};
                            Wn = node{node{m}.N(n)}.W(k,:).';
                            ALAW = ALAW + (Amn.'*(Lnm-Anm*Wn));
                        end
                        %                     node{m}.Wnew(k,:) = (rho*AA + squeeze(R(k,:,:)))\(ALAW + dm);
                        Rtmp = squeeze(node{m}.R(k,:,:));
                        %                     Rtmp = [Rtmp,Rtmp(:,m)];
                        %                     Rtmp = [Rtmp;Rtmp(m,1:M),Rtmp(m,m)];
                        %                     node{m}.Wnew(k,:) = (rho*AA + squeeze(R(k,:,:)))\(ALAW + dm);
                        node{m}.Wnew(k,:) = (rho*AA + Rtmp)\(ALAW + dm);
                        
                        % Check pdmm using ||Bx-b||^2
                        %                     node{m}.Wnew(k,:) = (AA+B(m,:).'*B(m,:))\(ALAW+B(m,:).'*b(m));
                        
                        % Lambda update
                        for n=1:Nlen
                            Amn = node{m}.Amn{n};
                            Anm = node{node{m}.N(n)}.Amn{node{node{m}.N(n)}.N==m};
                            Wn = node{node{m}.N(n)}.W(k,:).';
                            Wm = node{m}.Wnew(k,:).';
%                             if Nlen ==1
%                                 node{m}.Lnew(k,:) = node{node{m}.N(n)}.L(k,:).' - alpha*(Anm*Wn + Amn*Wm);
%                             else
                                node{m}.Lnew(k,:,n) = node{node{m}.N(n)}.L(k,:,mi).' - alpha*(Anm*Wn + Amn*Wm);
%                             end
                        end
                    end
                end
            end
            
            % Full Spectrum - Save the new weights to the nodes
            for m=activeNodes.'%1:M
                dtmp(1,m) = node{m}.d(k);
            end
            for m=activeNodes.'%1:M
                node{m}.W = node{m}.Wnew;
                node{m}.L = node{m}.Lnew;
                %             ftmp(iter1,m) = 0.5*((node{m}.W(k,:)))*squeeze(R(k,:,:))*(node{m}.W(k,:).')-(dtmp)*(node{m}.W(k,:).');
%                 WsaveAll(iter1,m,:) = node{m}.W(k,:);
%                 LsaveAll(iter1,m,:,:) = node{m}.L(k,:,:);
            end
            W = zeros(Khalf,M);
            for m=activeNodes.'
                Wtmp = [];
                for n=1:node{m}.Nlen
                    Wtmp = cat(2,Wtmp,node{node{m}.N(n)}.W(:,node{node{m}.N(n)}.N(1:end)==m));
                end
                
                W(:,m) = mean(Wtmp,2);
            end
            for k=bin
                f(iter1,k) = 0.5*W(k,:)*squeeze(R(k,:,:))*W(k,:)'-d(k,:)*W(k,:)'+norm(W(k,:),1);
            end
        end
    end
        
%     % Get stft of the sources
%     Target = stft(s{1}(1:tl),K);
%     Target = Target(2:(K-1)/2,:);
%     Target = Target(bin,:);
%     Interferer = stft(s{2}(1:tl),K);
%     Interferer = Interferer(2:(K-1)/2,:);
%     Interferer = Interferer(bin,:);
%     dT = exp(-1i*2*pi*fdomShort(bin)*ssd(1,mindisti)/c)/(4*pi*ssd(1,mindisti));
%     dI = exp(-1i*2*pi*fdomShort(bin)*ssd(2,mindisti)/c)/(4*pi*ssd(2,mindisti));
    
%     iSIR = 10*log10((Target*dT*dT'*Target')/(Interferer*dI*dI'*Interferer'));
%     iSIRtf = 10*log10((dT'*dT)/(dI'*dI));
%     iSIRnotf = 10*log10((Target*Target')/(Interferer*Interferer'));
    
    % Find output SIR
    weights = W(bin,activeNodes);
    dT = exp(-1i*2*pi*fdomShort(bin)*ssd(1,activeNodes)/c)./(4*pi*ssd(1,activeNodes));
    dI = exp(-1i*2*pi*fdomShort(bin)*ssd(2,activeNodes)/c)./(4*pi*ssd(2,activeNodes));
%     xT = (dT.'*Target).';
%     xI = (dI.'*Interferer).';
%     num = sum(repmat(weights',1,L-2).*xT.');
%     den = sum(repmat(weights',1,L-2).*xI.');
    
%     oSIR = 10*log10((num*num')/(den*den'));
    oSIRtf(ani) = 10*log10((conj(dT)*weights.'*conj(weights)*dT.')/(conj(dI)*weights.'*conj(weights)*dI.'));
%     SIRGain = oSIR-iSIR;
%     SIRGaintf = oSIRtf-iSIRtf;
    
%     table([iSIR;oSIR;SIRGain],[iSIRtf;oSIRtf;SIRGaintf],[iSIRnotf;0;0],'RowNames',{'Input';'Output';'Gain'},'VariableNames',{'SIR';'SIRtf';'SIRnotf'})

end

figure; plot(abs(oSIRtf)); grid on; 

%% Calculate BF output
% W = mean(cat(3,node{1}.W,node{2}.W,node{3}.W),3);
% for l=1:L
%     Y(:,l) = (1/M)*sum(squeeze(conj(W)).*squeeze(X(:,l,:)),2);
% end
% Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
% y = myOverlapAdd(Y);
% figure; plot(y); grid on; title('BF output y');


%% MSE between W and Wopt
% a = length(Wsave(1,:,1))
% WWoptMSE = zeros(Ltmp,1);
% for b=1:a
%     Wtmp = squeeze(Wsave(:,b,:));
%     WWoptMSE(b) = mean(mean((Wtmp-Wopt).*conj(Wtmp-Wopt)));
% end
% figure; semilogy(WWoptMSE); grid on; title('WWoptMSE');

%% W vs Wopt full spectrum
figure; imagesc(abs(Wopt)); title('Wopt');
figure; imagesc(abs(squeeze(W))); title('W');

%% MSE between W and Wopt single bin
% WSaveTmp = squeeze(Wsave(bin,:,:));
% WOptTmp = squeeze(Wopt(bin,:));
% LL = length(WSaveTmp);
% WWoptMSE = zeros(a,1);
% for b=1:a
%     WWoptMSE(b) = mean((WOptTmp-WSaveTmp(b,:))*(WOptTmp-WSaveTmp(b,:))');
% end
% figure ; semilogy(WWoptMSE); grid on; title('WWoptMSE single bin');

%% Variance of the sensor weights
VarWsave = zeros(M,1);
VarWopt = zeros(M,1);
for m=1:M
    VarWsave(m) = W(:,m)'*W(:,m);
    VarWopt(m) = Wopt(:,m)'*Wopt(:,m);
end

figure; plot(VarWsave,'*--'); grid on; hold on; plot(VarWopt,'o--'); title('Variance in sensors'); legend('VarWsave','VarWopt');
ratio = max(VarWsave)/max(VarWopt)
figure; plot((1/ratio)*VarWsave,'*--'); grid on; hold on; plot(VarWopt,'o--'); title('Variance in sensors with Wopt scaled'); legend('VarWsave','VarWopt');


%% Print setup for records
% fprintf('Nsrcs = %d, K = %d, tls = %d, M = %d, spSize = %d, bin = %d, Ltmp = %d\n\n',Nsrcs,K,tls,M,spSize,bin,Ltmp);

%%
% Wopt(54,1)'*node{1}.d(54)+Wopt(54,2)'*node{2}.d(54)
% Wopt(123,1)'*node{1}.d(123)+Wopt(123,2)'*node{2}.d(123)

%%
% fopt = 0.5*squeeze((Wopt(bin,:)))*squeeze(R(bin,:,:))*(squeeze(Wopt(bin,:)).')-(dtmp)*(squeeze(Wopt(bin,:)).');
% figure; plot(abs(ftmp(:,1)),'o--'); hold on; plot(abs(ftmp(:,2)),'.--'); plot(abs(ftmp(:,1)),'^--');
% plot(repmat(abs(fopt),ITER1,1));
% grid on; legend('node 1','node 2','node 3','fopt'); title('objective fn at each node');


%%
% % Find ||xi-xiopt|| for all i, note that xiopt is the same for all i
% xi_xiopt_norm = zeros(ITER1,M);
% for iter1=1:ITER1
%     for m=1:M
%         xi_xiopt_norm(iter1,m) = norm(squeeze(WsaveAll(iter1,m,:)) - squeeze(Wopt(bin,:)).');
% %         xi_xiopt_norm(iter1,m) = norm(squeeze(WsaveAll(iter1,m,:)) - (Wopt));
%     end
% end
% % figure; semilogy(xi_xiopt_norm(:,1), '.--'); hold on; semilogy(xi_xiopt_norm(:,2),'.--');
% % semilogy(xi_xiopt_norm(:,3),'.--'); grid on; legend('node 1','node 2','node 3'); title('norm(xi-xopt) for i=1,2,3, single bin');
% figure; plot(xi_xiopt_norm(:,1), '.--'); hold on; plot(xi_xiopt_norm(:,2),'.--');
% plot(xi_xiopt_norm(:,3),'.--'); grid on; legend('node 1','node 2','node 3'); title('norm(xi-xopt) for i=1,2,3, single bin');
% 
% 
% % Find ||x1-x2|| + ||x1-x3|| + ||x2-x3|| = 0
% xNormSum = zeros(ITER1,1);
% x1 = squeeze(WsaveAll(:,1,:));
% x2 = squeeze(WsaveAll(:,2,:));
% x3 = squeeze(WsaveAll(:,3,:));
% for iter1=1:ITER1
%     xNormSum(iter1) = norm(x1(iter1,:)-x2(iter1,:)) + norm(x1(iter1,:)-x3(iter1,:)) + norm(x2(iter1,:)-x3(iter1,:));
% end
% % figure; semilogy(xNormSum); grid on; title('norm(x1-x2)+norm(x1-x3)+norm(x2-x3), single bin');
% figure; plot(xNormSum); grid on; title('norm(x1-x2)+norm(x1-x3)+norm(x2-x3), single bin');
% 
% fprintf('\nfinal norm(x1-x2)+norm(x1-x3)+norm(x2-x3) = %d\n',xNormSum(end));
% 
% %%
% % Look at the L's
% figure; hold on;
% for m=1:M
%     for n=1:M
%         plot((abs(LsaveAll(:,m,1,n))));
%     end
% end
% grid on; legend('11','12','13','21','22','23','31','32','33'); title('Lambda values');

%% Find input SIR
% [mindist,mindisti] = min(ssd(1,:));
% 
% % Get stft of the sources
% Target = stft(s{1}(1:tl),K);
% Target = Target(2:(K-1)/2,:);
% Target = Target(bin,:);
% Interferer = stft(s{2}(1:tl),K);
% Interferer = Interferer(2:(K-1)/2,:);
% Interferer = Interferer(bin,:);
% dT = exp(-1i*2*pi*fdomShort(bin)*ssd(1,mindisti)/c)/(4*pi*ssd(1,mindisti));
% dI = exp(-1i*2*pi*fdomShort(bin)*ssd(2,mindisti)/c)/(4*pi*ssd(2,mindisti));
% 
% iSIR = 10*log10((Target*dT*dT'*Target')/(Interferer*dI*dI'*Interferer'));
% iSIRtf = 10*log10((dT'*dT)/(dI'*dI));
% iSIRnotf = 10*log10((Target*Target')/(Interferer*Interferer'));

%% Find output SIR
% weights = W(bin,:);
% dT = exp(-1i*2*pi*fdomShort(bin)*ssd(1,:)/c)./(4*pi*ssd(1,:));
% dI = exp(-1i*2*pi*fdomShort(bin)*ssd(2,:)/c)./(4*pi*ssd(2,:));
% xT = (dT.'*Target).';
% xI = (dI.'*Interferer).';
% num = sum(repmat(weights',1,L-2).*xT.');
% den = sum(repmat(weights',1,L-2).*xI.');
% 
% oSIR = 10*log10((num*num')/(den*den'));
% oSIRtf = 10*log10((conj(dT)*weights.'*conj(weights)*dT.')/(conj(dI)*weights.'*conj(weights)*dI.'));
% SIRGain = oSIR-iSIR;
% SIRGaintf = oSIRtf-iSIRtf;
% 
% table([iSIR;oSIR;SIRGain],[iSIRtf;oSIRtf;SIRGaintf],[iSIRnotf;0;0],'RowNames',{'Input';'Output';'Gain'},'VariableNames',{'SIR';'SIRtf';'SIRnotf'})



