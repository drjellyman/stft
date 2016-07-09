close all; clear
% distance vs sample delay

fs = 16e3;
distance = [0.01:0.001:100]';
distanceLen = length(distance);
c = 343; 

closeToInt = zeros(distanceLen,1);
for d = 1:distanceLen
    delayInSamples = (distance(d)/c)/(1/fs);
    closeToInt(d) = mod(delayInSamples,1);
end

figure; plot(distance,closeToInt,'.'); grid on;
