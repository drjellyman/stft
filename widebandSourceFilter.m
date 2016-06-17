close all; 
clear all; 


% lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 0.45,...
%     'StopbandFrequency', 0.55, 'PassbandRipple', 1, ...
%     'StopbandAttenuation', 60, 'DesignMethod', 'equiripple');

lpFilt = designfilt('lowpassfir', 'SampleRate', 16000, ...
                'PassbandFrequency', 6e3,...
      'StopbandFrequency', 8e3,...
         'PassbandRipple', 1,...
    'StopbandAttenuation', 60,...
           'DesignMethod', 'equiripple');

s = filter(lpFilt, randn(16e3,1))  ;
figure; plot(s);
figure; plot(abs(fft(s)));
d = 7; % d = delay
M = 8; % M = number of sensors
x = zeros(length(s),M); % x = observation signals
for m=1:M
    x(d*(m-1)+1:end,m) = s(1:end-d*(m-1)); 
end
figure; plot(abs(fft(s))); legend('s'); set(gca,'fontsize', 14);
xsum = sum(x,2);
figure; plot(abs(fft(xsum))); legend('xsum'); set(gca,'fontsize', 14);

