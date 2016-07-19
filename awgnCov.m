close all; clear

N = 10;
fs = 16e3;
K = 513;
Khalf = (K-1)/2-1;
T = 1; % length in time
S = T*fs-mod(T*fs,(K-1)/2);; % length in samples
x = randn(S,N);
[X,L] = stft(x,K);
X = X(2:(K-1)/2,:,:);
Rf = zeros(Khalf,L,N,N);
RfRunSum = zeros(Khalf,N,N);
alpha = 0.5;
runRcond = zeros(L,3);
for l=1:L
    for k=1:Khalf
        Rf(k,l,:,:) = squeeze(X(k,l,:))*squeeze(X(k,l,:)).';
        RfRunSum(k,:,:) = alpha*squeeze(RfRunSum(k,:,:)) + (1-alpha)*(squeeze(Rf(k,l,:,:))); 
        
    end
    runRcond(l,:) = [rcond(squeeze(RfRunSum(34,:,:))),rcond(squeeze(RfRunSum(12,:,:))),rcond(squeeze(RfRunSum(145,:,:)))];
end
figure; plot(runRcond)

rc1 = rcond(squeeze(Rf(23,41,:,:)))
RfSum = squeeze(sum(Rf,2));
rc2 = rcond(squeeze(RfSum(23,:,:)))
RfklSum = squeeze(sum(RfSum,1));
rc3 = rcond(RfklSum)

%% What about time domain spatial covariance? 
Rt = zeros(T,N,N);
for s=1:S
    Rt(s,:,:) = x(s,:).'*x(s,:);
end
rc4 = rcond(squeeze(Rt(14,:,:)))
RtSum = squeeze(sum(Rt,1));
rc5 = rcond(RtSum)

