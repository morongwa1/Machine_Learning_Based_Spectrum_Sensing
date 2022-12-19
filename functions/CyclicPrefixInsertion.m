function [ofdmData] = CyclicPrefixInsertion(ofdmSymb, numCPSamples)
%CyclicPrefixInsertion - adds cyclic prefix to ofdm data
%   Parameters:
%       ofdmSymb - time ofdm samples (complex values). Every column is an
%           ofdm symbol
%       numCPSamples - length of cyclic prefix
%
%   Output:
%       ofdmData - ofdm data with cyclic prefix

ind1 = length(ofdmSymb)-numCPSamples+1;
ind2 = length(ofdmSymb);

CP = ofdmSymb(ind1:ind2,:);
ofdmData = [CP; ofdmSymb];

end

