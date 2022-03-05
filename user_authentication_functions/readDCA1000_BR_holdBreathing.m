function [ retVal] = readDCA1000_BR_holdBreathing(fileName, adcsamples)


%% global variables
% change based on sensor config
numADCSamples = adcsamples;     % number of ADC samples per chirp
numADCBits = 16;                % number of ADC bits per sample
numRX = 1;                      % number of receivers
numLanes = 2;                   % do not change. number of lanes is always 2
isReal = 0;                     % set to 1 if real only data, 0 if complex data0

%% read file
% read .bin file
fid = fopen(fileName,'r');
adcData = fread(fid, 'int16');

fclose(fid);
fileSize = size(adcData, 1);



numChirps = fileSize/2/numADCSamples/numRX;

LVDS = zeros(1, fileSize/2);
%combine real and imaginary part into complex data
%read in file: 2I is followed by 2Q

counter = 1;

for i=1:4:fileSize-1
% LVDS(1,counter) = adcData(i+2) + sqrt(-1)*adcData(i); 
% LVDS(1,counter+1) = adcData(i+3)+sqrt(-1)*adcData(i+1); 
LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
counter = counter + 2;
end


% create column for each chirp
LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
%each row is data from one chirp
LVDS = LVDS.';

%organize data per RX
adcData = zeros(numRX,numChirps*numADCSamples);

for row = 1:numRX
for i = 1: numChirps
adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
end
end


% return receiver data
retVal = LVDS;%adcData;%;LVDS;%adcData;
end