function [ y ] = myNormalize( x )
    % Takes vector x and scales it to have amplitude -1 <= A <= 1
    
    m = max(abs(x));
    y = (1/m)*x;
end

