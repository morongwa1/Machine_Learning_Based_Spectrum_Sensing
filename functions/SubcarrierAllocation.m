function [ind1, ind2] = SubcarrierAllocation(blocks, FFTLength, numSubcarriers, resBlocksNum)
%SubcarrierAllocation - allocates resource blocks (frequency domain) for LTE user
%   Parameters:
%       blocks - resource blosks that are to be assigned (this should be a
%           vector [a b] - assigning block from number a to b. When only one
%           block is to be assign, use a=b. Maximum b value = resBlocksNum.
%           Minimum a value = 1).
%       FFTLength - 
%       numSubcarriers - 
%       resBlocksNum - 
%
%   Output:
%       ind1 - first sample index (frequency domain)
%       ind2 - last sample index (frequency domain)
%

if length(blocks) ~= 2
    error('wrong size of input "block" vector')
end

if (blocks(2)-blocks(1)+1) > resBlocksNum
    error(['In this band it is impossible to use more than ' num2str(resBlocksNum) ' resource block.'])
end

if blocks(2) > resBlocksNum
    error(['The second value in input vector blocks must be equal or smaller to ' num2str(resBlocksNum)])
end

if blocks(1) > blocks(2)
    error('Second blocks vector value must be bigger or equal to first value')
end


usableSamples = numSubcarriers*resBlocksNum;               
marginSamples = FFTLength - usableSamples;                 

firstFreqIdx = 0.5*marginSamples+1;                       
lastFreqIdx = 0.5*marginSamples + usableSamples;          

ind1 = firstFreqIdx + (blocks(1)-1)*numSubcarriers;      
ind2 = firstFreqIdx + blocks(2)*numSubcarriers - 1;       

if ind2 > lastFreqIdx
    error('exceeded frequency range')
end

end

