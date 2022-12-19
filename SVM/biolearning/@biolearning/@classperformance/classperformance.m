function h = classperformance(CL,GT)
%CLASSPERFORMANCE class constructor for CP objects

% $Revision: 1.1.10.2 $    $Date: 2006/06/16 20:07:16 $
% Copyright 2003-2006 The MathWorks, Inc.

% construct a CP object

h.ListenerEnable = false;
h = biolearning.classperformance;
h.ClassLabels = CL;
h.GroundTruth = GT;
h.IsClassLabelTypeNumeric = isnumeric(CL);
h.NumberOfObservations = numel(GT);
h.ValidationCounter = 0;
numClasses = numel(CL);
h.TargetClasses = 1;
h.ControlClasses = (2:numClasses)';
h.SampleDistribution = zeros(h.NumberOfObservations,1);
h.ErrorDistribution = zeros(h.NumberOfObservations,1);
numClasses = numel(CL);
h.CountingMatrix = zeros(numClasses+1,numClasses);
h.SampleDistributionByClass = zeros(numClasses,1);
h.ErrorDistributionByClass = zeros(numClasses,1);
h.ListenerEnable = true;
