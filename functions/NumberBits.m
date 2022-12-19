function [numBitsPerSymb] = NumberBits(modType)
%NUMBERBITS - returns number of bits used per symbol in modulation of 
%               type modType
%   Parameters:
%       modType - modulation type (QPSK or 16QAM or 64QAM)
%   Output:
%       numBitsPerSymb - number of bits per modulated symbol
%

if strcmp(modType, 'QPSK')
    numBitsPerSymb = 2;
elseif strcmp(modType, '16QAM')
    numBitsPerSymb = 4;
elseif strcmp(modType, '64QAM')
    numBitsPerSymb = 6;
elseif strcmp(modType, 'BPSK')
    numBitsPerSymb = 1;
else
    error('Incorrect modType argument value')
end
end

