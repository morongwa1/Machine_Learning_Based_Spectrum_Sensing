function [modData] = Modulator(a_ModType, a_Bits)
%MODULATOR applies modulation to input bit stream
%   Parameters:
%       a_ModType - modulation type (Possible modulation types: QPSK, 
%                                           16QAM, 64QAM)
%       a_Bits - sequence of bits to modulate
%       a_Coding - 0 for binary coding, 1 for Grey coding
%   Output:
%       modData - modulated data symbols
%
    
if ~iscolumn(a_Bits)
    a_Bits = a_Bits';
end

if strcmp(a_ModType, 'QPSK')
    k = 2;                                                                 % Number of bits per symbol
    remainder = mod(length(a_Bits),k);
    if mod(length(a_Bits),k) ~= 0
        disp(['Number of input bits is not a multiple of ', num2str(k), '(number of bits per symbol)'])
        disp(['Number of last bits ignored in further processing: ', num2str(remainder)])
        a_Bits = a_Bits(1:length(a_Bits)-remainder);
    end
    
    qpskMod = comm.QPSKModulator('BitInput', true);
    modData = qpskMod(a_Bits);


elseif strcmp(a_ModType, 'BPSK')
    k = 1;                                                                 % Number of bits per symbol
    remainder = mod(length(a_Bits),k);
    if mod(length(a_Bits),k) ~= 0
        disp(['Number of input bits is not a multiple of ', num2str(k), '(number of bits per symbol)'])
        disp(['Number of last bits ignored in further processing: ', num2str(remainder)])
        a_Bits = a_Bits(1:length(a_Bits)-remainder);
    end
    
    bpskMod = comm.BPSKModulator();
    modData = bpskMod(a_Bits);
    
    
elseif strcmp(a_ModType, '16QAM')
    k = 4;                                                                 % Number of bits per symbol
    remainder = mod(length(a_Bits),k);
    if mod(length(a_Bits),k) ~= 0
        disp(['Number of input bits is not a multiple of ', num2str(k), '(number of bits per symbol)'])
        disp(['Number of last bits ignored in further processing: ', num2str(remainder)])
        a_Bits = a_Bits(1:length(a_Bits)-remainder);
    end
    
    bitsMatrix = reshape(a_Bits, length(a_Bits)/k, k);                     % Reshape data into binary k-tuples
    bitsSymbols = bi2de(bitsMatrix);                                       % Convert to integers
    lteSymMap = [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
    modData = qammod(bitsSymbols, 16, lteSymMap, 'UnitAveragePower', true);
    
elseif strcmp(a_ModType, '64QAM')
    k = 6;                                                                 % Number of bits per symbol
    remainder = mod(length(a_Bits),k);
    if mod(length(a_Bits),k) ~= 0
        disp(['Number of input bits is not a multiple of ', num2str(k), '(number of bits per symbol)'])
        disp(['Number of last bits ignored in further processing: ', num2str(remainder)])
        a_Bits = a_Bits(1:length(a_Bits)-remainder);
    end
    
    bitsMatrix = reshape(a_Bits, length(a_Bits)/k, k);                     % Reshape data into binary k-tuples
    bitsSymbols = bi2de(bitsMatrix);                                       % Convert to integers
    modData = qammod(bitsSymbols, 64);
    
else
    error('Incorrect a_ModType argument value')
end

end

