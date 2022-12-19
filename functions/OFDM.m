function [ofdmData, dataPower] = OFDM(modData, FFTLength, numCPSamples, firstSubcarrier)
%OFDM - calculates OFDM symbols
%   Parameters:
%       modData - data symbols (complex) - data bits modulated with QPSK 
%                       or 16QAM or 64QAM modulation
%       FFTLength - number of FFT samples
%       numCPSamples - number of samples that will be used as cyclic prefix
%       firstSubcarrier - index of first frequency
%       
%   Output:
%       ofdmSymbols - stream of OFDM symbols
%

[rowModD,colModD] = size(modData);
if FFTLength > rowModD
    modData = [zeros(firstSubcarrier-1, colModD); modData; zeros(FFTLength-rowModD-firstSubcarrier+1, colModD)];
end

dataPower = sum(abs(modData).^2)/length(modData);

ofdmData = ifft(modData, FFTLength);
ofdmData = CyclicPrefixInsertion(ofdmData, numCPSamples);
end

