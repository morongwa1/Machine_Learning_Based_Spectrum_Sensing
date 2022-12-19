function [LTETransmData] = transmissionLTE(ofdmVect, firstTimeSlot, numSymbols, LTETransmData, FFTLength, numCPSamples, ofdmSymbPerTimeSlot)
%transmissionLTE adds ofdm data to LTETransmData vector
%   Parameters:
%       ofdmVect - generated ofdm data
%       firstTimeSlot - first time slot, when ofdmVect is send
%       numSymbols - number of ofdm symbols
%       LTETransmData - vector of transmitted signals
%       FFTLength - number of fft samples
%       numCPSamples - number of cyclic prefix samples
%   Output:
%       LTETransmData - modified vector of transmitted signals

firstInd = (firstTimeSlot-1) * ofdmSymbPerTimeSlot * (FFTLength+numCPSamples) + 1;
lastInd = firstInd + numSymbols * ofdmSymbPerTimeSlot * (FFTLength+numCPSamples) - 1;

LTETransmData(firstInd:lastInd) = LTETransmData(firstInd:lastInd) + ofdmVect;


end

