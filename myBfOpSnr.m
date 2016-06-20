function [ y,SNRdB ] = myBfOp( X,Xt,Xi,W )
% This function calculates the output of a beamformer and the SNR. It finds
% the SNR by creating outputs for independently for the two sources and
% comparing their powers.

    % Initialize output matrices
    Y = zeros(Khalf,L);
    Yt= zeros(Khalf,L);
    Yi= zeros(Khalf,L);
    
    % Calculate outputs
    for l=1:L
        Xtmp = squeeze(X(:,l,:));
        Xttmp = squeeze(Xt(:,l,:));
        Xitmp = squeeze(Xi(:,l,:));
        Y(:,l) = sum(conj(W).*Xtmp,2);
        Yt(:,l) = sum(conj(W).*Xttmp,2);
        Yi(:,l) = sum(conj(W).*Xitmp,2);
    end
    
    % Make 2 sided spectrums
    Y = [zeros(1,L);Y;zeros(2,L);conj(flipud(Y))];
    Yt = [zeros(1,L);Yt;zeros(2,L);conj(flipud(Yt))];
    Yi = [zeros(1,L);Yi;zeros(2,L);conj(flipud(Yi))];
    
    % Calculate time domain signals
    y = myOverlapAdd(Y);
    yt = myOverlapAdd(Yt);
    yi = myOverlapAdd(Yi);
    
    % Calculate SNRdB
    SNRdB = 10*log10(((yt'*yt)/length(yt)) / ((yi'*yi)/length(yi)));
end

