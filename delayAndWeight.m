function [ y ] = delayAndWeight( s,ssd,Fs )

    velSnd = 334; % The velocity of sound in air (ms^-1)
    S = fft(s);
    sd = (ssd/velSnd)/(1/Fs) % sd = sample delay
    w = 1/((ssd+1)^2); % I've put the plus one in here to avoid amplifying the signal, but it seems a bit weird
    N = length(S);
    
    % This shifts the signal by a non integer amount. 
    for k = 1:N
        if k < N-k
            SShifted(k) = S(k) * exp(2*pi*i * sd*k/N);
        elseif k > N-k
            SShifted(k) = S(k) * exp(2*pi*i * sd*(k-N)/N);
        elseif k == N/2
            SShifted(k) = S(k) * cos(2*pi * sd/2);
        end
    end
    y = ifft(SShifted); % This actually results in a complex signal, which leads to the output being attenuated when only the real part is considered. Best I've got right now though. 

end

