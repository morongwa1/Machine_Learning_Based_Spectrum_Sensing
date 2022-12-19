function [LTETransmData, powers] = rescaleLTEData(LTETransmData, resourceAlloc, FFTLength, numCPSamples, ofdmSymbPerTimeSlot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

resources = sum(sign(resourceAlloc));
resources = resources(1:end-1);

powers = zeros(1,length(resources));

numTimeSlots = length(LTETransmData)/((FFTLength+numCPSamples)*ofdmSymbPerTimeSlot);

if length(resources) ~= numTimeSlots
    error('wrong length of LTETransmData');
end

tInd1 = 1;
tInd2 = tInd1 + ((FFTLength+numCPSamples)*ofdmSymbPerTimeSlot) - 1;
numTSamps = tInd2 - tInd1 + 1;

for k = 1:numTimeSlots         %normalizing power of all the resource blocks
    if resources(k) > 1
        LTETransmData(tInd1:tInd2) = LTETransmData(tInd1:tInd2)/sqrt(resources(k));
    end
    powers(k) = sum(abs(LTETransmData(tInd1:tInd2)).^2)/numTSamps;
    tInd1 = tInd2 + 1;
    tInd2 = tInd1 + ((FFTLength+numCPSamples)*ofdmSymbPerTimeSlot) - 1;
end

