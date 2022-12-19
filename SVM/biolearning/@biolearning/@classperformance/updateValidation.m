function updateValidation(CP,gps,idx)
%UPDATEVALIDATION updates the CP object after one validation of the classifier

% Copyright 2004-2006 The MathWorks, Inc.

numClasses = numel(CP.ClassLabels);
gps(gps == 0) = numClasses + 1; % indexes for Inconclusive results (pass them to the last)
CP.ValidationCounter = CP.ValidationCounter + 1;
errors = CP.GroundTruth(idx)~=gps;
CP.ErrorDistribution(idx) = CP.ErrorDistribution(idx) + errors;
CP.SampleDistribution(idx) = CP.SampleDistribution(idx) + 1;
CP.ErrorDistributionByClass = CP.ErrorDistributionByClass + ...
                  accumarray(CP.GroundTruth(idx),double(errors),[numClasses,1]);
CP.SampleDistributionByClass = CP.SampleDistributionByClass + ...
                  accumarray(CP.GroundTruth(idx),1,[numClasses,1]);
for k = 1:numClasses
   if ~isempty(gps(CP.GroundTruth(idx)==k))
    CP.CountingMatrix(:,k) = CP.CountingMatrix(:,k) + ...
       accumarray(gps(CP.GroundTruth(idx)==k),1,[numClasses+1,1]);
   end
end
n = sum(gps<=numClasses);
if n == 0
    CP.LastErrorRate = NaN;
    CP.LastCorrectRate = NaN;
else
    CP.LastCorrectRate = sum(~errors)/n;
    CP.LastErrorRate = 1 - CP.LastCorrectRate;
end