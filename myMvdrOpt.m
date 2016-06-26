function [ W ] = myMvdrOpt( A,X,dl )
% Calculates the optimum MVDR weights. A is the target signal ATF, X is the
% observations, dl is the diagonal loading factor, and W is the weights. 
    
    % Initialize vars
    Khalf = length(X(:,1,1));
    M = length(X(1,1,:));
    L = length(X(1,:,1));

    % Initialize matrices
    W = ones(Khalf,M);
    for k=1:Khalf
        R = zeros(M,M); % R is the spatial covariance of the inputs
        for l=1:L
            Xtmp = squeeze(X(k,l,:));
            R = R + Xtmp*Xtmp'; % Sum the covariance over l
        end
        R = R + dl*eye(M); % Diagonal loading
        Ak = A(k,:).';
        W(k,:) = (R\Ak)/(Ak'*(R\Ak)); % Calculate optimum weights vector 
    end
end

