function [bits] = BitStreamGen(a_NumBits, a_Polar)
%BITSTREAMGEN generates sequence of bits
%   Parameters:
%       a_NumBits - number of bits to generate
%       a_Polar - polarization of bits (select 0 for unipolar or 1 for 
%                       bipolar)
%   Output:
%       bits - stream of bits
%

if and(a_Polar ~= 0, a_Polar ~= 1)
    error('Incorrect a_Polar argument value')
end

bits = randi([0 1],[1,a_NumBits]);                                         % if a_Polar == 0

if a_Polar == 1
    bits = 2*bits - 1;
end
end

