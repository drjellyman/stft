function [ S ,L] = stft( s, K )
    % Take the short time fourier transform of x and return the
    % time-frequency domain version X. K is the length of the window and
    % should be odd. L is the number of windows used. 
    
    % Check if K is odd, add one if not
    if mod(K,2) ~= 1
        K = K+1; 
    end
    
    
        
    % Hann window for stft
    w = sqrt(0.5*(1-cos((2*pi*[0:K-1]')/(K-1)))); 
    
    if isa(s,'cell') 
        % Find number of sensors
        M = length(s);
        
        % Calculate number of windows
        L = (length(s{1})/(K-1))*2-1;

        % Calculate STFT
        sp = 1 + floor(K/2) * [0:L-1]'; % sp = start points of each window
        for m = 1:M
            s{m} = [s{m};0];
            for l = 1:L
                S{m}(:,l) = fft(s{m}(sp(l):sp(l)+K-1) .* w(1:end));
            end
        end
    else 
        % Find number of sensors
        M = length(s(1,:));
        
        % Calculate number of windows
        L = (length(s)/(K-1))*2-1;

        % Calculate STFT
        sp = 1 + floor(K/2) * [0:L-1]'; % sp = start points of each window
        for m = 1:M
            s = [s;zeros(1,M)];
            for l = 1:L
                S(:,l,m) = fft(s(sp(l):sp(l)+K-1,m) .* w(1:end));
            end
        end
    end
end

