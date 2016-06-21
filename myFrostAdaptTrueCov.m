function [ Y,W,Wmse ] = myFrostAdaptTrueCov( A,X,mu,Iter,Wopt )
    % computes the output of a Frost MVDR beamformer using the true signal
    % covariance. 
    %
    % Inputs: 
    % A = target source ATF (K x M)
    % X = observation signals (K x L x M)
    % mu = step size
    % Iter = number of iterations per window
    % Wopt = mvdr optimal weights (K x M)
    %
    % Outputs
    % Y = Beamformer output (K x L)
    % W = Adapted weights matrix (K x M)
    % Wmse = mean squared error between W and Wopt
    %
    % Indices
    % k,Khalf = frequency bins (half = one sided frequency domain)
    % l,L = window index (time)
    % m,M = sensor index
    
    % Initialize vars
    Khalf = length(X(:,1,1));
    M = length(X(1,1,:));
    L = length(X(1,:,1));
    
    % Initialize matrices
    Wmse = zeros(L,1);
    Y = zeros(Khalf,L);
    R = zeros(Khalf,M,M);
    W = ones(Khalf,M);
    
    % Setup precomputable matrices
    for k = 1:Khalf
        Ak = A(k,:).';
        P(k,:,:) = eye(M) - (Ak*Ak')/(norm(Ak)^2);
        F(k,:) = Ak/(norm(Ak)^2);
    end
    
    % Initialize weight vector
    W = F; 
    
    % Calculate the true covariance
    for l=1:L
        for k=1:Khalf
            R(k,:,:) = squeeze(R(k,:,:)) + (1/L)*squeeze(X(k,l,:))*squeeze((X(k,l,:)))';
        end
    end
    
    % Iteratively update the weights W and output Y
    for l=1:L
        Xtmp = squeeze(X(:,l,:));
        Y(:,l) = sum(conj(W).*Xtmp,2);     
        for k = 1:Khalf
            Xtmp = squeeze(X(k,l,:));
            Ptmp = squeeze(P(k,:,:));
            Ftmp = F(k,:).';
            Rtmp = squeeze(R(k,:,:));
            for iter = 1:Iter
                Wtmp = W(k,:).'; 
                W(k,:) = Ptmp*(Wtmp-mu*Rtmp*Wtmp)+Ftmp; 
            end
        end
        Wmse(l) = myMse(Wopt,W);
    end

end

