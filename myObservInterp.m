function [ z,d] = myObservInterp( zPos,sPos,s,Fs1,Fs,nsWt )
    % Takes the position arrays zPos(sensors) and sPos(sources) and creates
    % the observation signals z based on the source signals(s1,s1,s2) and
    % the distance between the source and the sensor. s1 is the source of
    % interest, s2 is an interferer, and s3 is random
    % noise (non-directional). The signals should be the same length.
    
    NSamples = length(s(:,1));
    NSources = 2;
    NSensors = length(zPos(1,:));
    z = zeros(NSamples, NSensors);
    NOrder = 3;
    c = 343; % c = velocity of sound in air (m.s^-1)
    d = zeros(NSources, NSensors);
    for a = 1:NSources
        for b = 1:NSensors
            ssd = norm(zPos(:,b) - sPos(:,a)); % ssd = source sensor distance (m)
            ssds = (ssd/c) / (1/Fs1); % Delay between the sensor and source in samples
            ssdsInt = floor(ssds) - 1; % Split up the delay into an integer and a rational between 1 and 2
            ssdsSmall = ssds - ssdsInt; 
            weight = 1/(ssd+1)^2; % Calculate the weight based on the source to sensor distance
            sTemp = [s(ssdsInt+1:end,a) ; zeros(ssdsInt,1)]; % Delay the source by ssdsInt
            z(:,b) = z(:,b) + weight * resample(filter(lagrange(NOrder,ssdsSmall),1,sTemp),1,(Fs1/Fs)); % Delay the source by ssdsSmall, then weight and sum into the obervation signal
            
            d(a,b) = ssd; % keep track of the distances to be returned (required for computing A, the ATF)
        end
    end
    z = z + nsWt*randn(NSamples,NSensors);
end

