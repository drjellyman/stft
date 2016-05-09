% Randomly place an arbitrary number of sensors within a 10 x 5 x 3 m room.
% Randomly place two audio sources within the same room. 

close all; clear all;

% Variables
NSensors = 4; % The number of sensors
NSources = 2; % The number of sources
Fs = 16e3; % The sample rate at the sensors
roomDim = [10, 5, 3]';
velSnd = 334; % The velocity of sound in air (ms^-1)

% Import audio files
[s1,Fs1] = audioread('50_male_speech_english_ch10_orth_2Y.flac');
[s2,Fs2] = audioread('44_soprano_ch17_orth_1Z.flac');

% Check sample rates of the audio files
if Fs1 ~= Fs2 fprintf('The audio files have different sample rates\n'); end

% Truncate one of the files so they have equivalent lengths
if (length(s1) > length(s2)) s1 = s1(1:length(s1));
elseif (length(s1) < length(s2)) s2 = s2(1:length(s1));
end
sDur = length(s1)/Fs1; % sDur = source duration in seconds

% Randomly place sensors and sources
sensorLocation = diag(roomDim)*rand(3,NSensors); % Each column of sensor_locations gives x,y,z coordinates of a sensor.
sourceLocation = diag(roomDim)*rand(3,NSources);

% Calculate observation signals x_m, as a mixture of the two sources based
% on the distance between the source and the sensor
NSamples = sDur/(1/Fs);
x = zeros(NSamples, NSensors);
for a = 1:NSources
    for b = 1:NSensors
        ssd = norm(sensorLocation(:,b) - sourceLocation(:,a)) / velSnd; % ssd = source sensor delay
%         x(:,b) = x(:,b) + delayed and weighted source a
    end
end



