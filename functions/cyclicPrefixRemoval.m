function [data] = cyclicPrefixRemoval(data, FFTLength, numCPSamples)

rows = length(data)/(FFTLength+numCPSamples);
data_temp = reshape(data, [FFTLength+numCPSamples, rows])';
data_temp = data_temp(:,numCPSamples+1:FFTLength+numCPSamples);
data = reshape(data_temp', [1, rows*FFTLength]);
end