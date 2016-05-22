function [ ] = mySpectrogram( X )
    % Plots a spectrogram given the time-frequency domain data X 
    
    N = length(X(:,1));
    fs = 1/(N-1);
    ts = 1/(length(X(1,:))-1);
    t = 0:ts:1; % Normalized
    f = 0:fs:1; % Normalized

    
    % Plot the spectrogram
    X = X(1:N/2+1,:); % Truncate to one sided spectrum
    X(2:end-1) = 2*X(2:end-1); % Truncating halved the power, now double it to get back to equal. Exclude the DC and Nyquist components as there were only ever one of each
    X = (abs(X).^2); % Square to make power. There could be scaling done here too i.e. 1/N but iamgesc() will just undo it anyway
    figure;  imagesc(t,f,10*log10(X)); 
    xlabel('Normalized time'); ylabel('Normalized Frequency'); 
    axis xy; colorbar; set(gca,'fontsize',14);
end

