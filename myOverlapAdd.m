function [x ] = myOverlapAdd( X )
    % Reconstruct original signal x from time-frequency domain X
    
    K = length(X(:,1)) + 1;
    L = length(X(1,:));

    % Sqrt Hann window
    w = sqrt(0.5*(1-cos((2*pi*[0:K-1]')/(K-1))));
        
    % Take the ifft of each time block
    xBlocks = ifft(X);
    
    % Overlap add the blocks
    x = zeros((K-1)*(L+1)/2,1);
    sp = 1 + floor(K/2) * [0:(L-1)]'; % sp = start points of each window
    for l = 1:L
        x(sp(l):sp(l)+K-2) = x(sp(l):sp(l)+K-2) + (xBlocks(1:end,l).*w(1:end-1));
    end
end

