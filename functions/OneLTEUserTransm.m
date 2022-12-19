function [ofdmData, dataPower] = OneLTEUserTransm(numBits, polar, modType, numSymbols, FFTLength, numCPSamples, firstSubcarrier, ofdmSymbPerTimeSlot)
%OneLTEUserTransm generaters OFDM data according to input parameters
%   Parameters:
%       numBits - number of data bits
%       polar - polarity of data bits
%       modType - modulation type
%       numSymbols - number of OFDM symbols
%       FFTLength - number of FFT samples
%       numCPSamples - number of cyclic prefix samples
%       firstSubcarrier - number of first subcarrier of OFDM signal
%   Output:
%       ofdmData - generated ofdm signal

[bits] = BitStreamGen(numBits, polar);
[modData] = Modulator(modType, bits);

modData = vec2mat(modData, numSymbols*ofdmSymbPerTimeSlot);           %generating a resource regions of size allocated to the user

[ofdmData, dataPower] = OFDM(modData, FFTLength, numCPSamples, firstSubcarrier);

end

