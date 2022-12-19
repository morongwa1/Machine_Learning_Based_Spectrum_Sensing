function schema
%
% $Revision: 1.1.10.2 $    $Date: 2006/06/16 20:07:17 $
% Copyright 2003-2006 The MathWorks, Inc.

%% add classperf class to package
blpk = findpackage('biolearning');
cls = schema.class(blpk,'classperformance');

% user properties
p = schema.prop(cls,'Label','string');  %#ok
p = schema.prop(cls,'Description','string'); %#ok

% properties needed at initialization
p = schema.prop(cls,'IsClassLabelTypeNumeric','bool');   p.Visible = 'off';
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'ListenerEnable','bool');   p.Visible = 'off';
                                  p.AccessFlags.PublicSet = 'off';
                                  p.FactoryValue = false;
p = schema.prop(cls,'ClassLabels','MATLAB array');
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'GroundTruth','MATLAB array'); 
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'NumberOfObservations','double');
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'ControlClasses','MATLAB array'); %#ok
p = schema.prop(cls,'TargetClasses','MATLAB array'); %#ok


% properties for updating
p = schema.prop(cls,'ValidationCounter','double');  
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'SampleDistribution','MATLAB array');
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'ErrorDistribution','MATLAB array');
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'SampleDistributionByClass','MATLAB array');
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'ErrorDistributionByClass','MATLAB array');
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'CountingMatrix','MATLAB array');
                                  p.AccessFlags.PublicSet = 'off';

% properties for reporting 
p = schema.prop(cls,'CorrectRate','double'); 
                                  p.GetFunction = @calculateCorrectRate;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'ErrorRate','double');
                                  p.GetFunction = @calculateErrorRate;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'LastCorrectRate','double'); 
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'LastErrorRate','double');
                                  p.AccessFlags.PublicSet = 'off';                                  
p = schema.prop(cls,'InconclusiveRate','double'); 
                                  p.GetFunction = @calculateInconclusiveRate;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'ClassifiedRate','double');
                                  p.GetFunction = @calculateClassifiedRate;
                                  p.AccessFlags.PublicSet = 'off';                                  
p = schema.prop(cls,'Sensitivity','double');
                                  p.GetFunction = @calculateSensitivity;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'Specificity','double');
                                  p.GetFunction = @calculateSpecificity;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'PositivePredictiveValue','double');
                                  p.GetFunction = @calculatePositivePredictiveValue;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'NegativePredictiveValue','double');
                                  p.GetFunction = @calculateNegativePredictiveValue;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'PositiveLikelihood','double');
                                  p.GetFunction = @calculatePositiveLikelihood;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'NegativeLikelihood','double');
                                  p.GetFunction = @calculateNegativeLikelihood;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'Prevalence','double');
                                  p.GetFunction = @calculatePrevalence;
                                  p.AccessFlags.PublicSet = 'off';
p = schema.prop(cls,'DiagnosticTable','MATLAB array');
                                  p.GetFunction = @calculateDiagnosticTable;
                                  p.AccessFlags.PublicSet = 'off';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateErrorRate(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix(1:end-1,:);
M = size(CM,1);
N = sum(sum(CM));
if N ~= 0
    v = sum(sum(CM.*(1-eye(M))))/N;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateCorrectRate(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix(1:end-1,:);
M = size(CM,1);
N = sum(sum(CM));
if  N ~= 0
    v = 1-sum(sum(CM.*(1-eye(M))))/N;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateInconclusiveRate(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
N = sum(sum(CM));
if N ~= 0
    v = sum(CM(end,:))/N;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateClassifiedRate(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
N = sum(sum(CM));
if N ~= 0
    v = 1-sum(CM(end,:))/N;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateSensitivity(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
%negInd = h.ControlClasses;
TP = sum(sum(CM(posInd,posInd)));
FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
%FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
%TN = sum(sum(CM(negInd,negInd)));
if (TP+FN) ~= 0
    v = TP/(TP+FN);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateSpecificity(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
negInd = h.ControlClasses;
TP = sum(sum(CM(posInd,posInd)));
FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
TN = sum(sum(CM(negInd,negInd)));
if (TP+FN) ~= 0
    v = TN/(FP+TN);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculatePositivePredictiveValue(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
negInd = h.ControlClasses;
TP = sum(sum(CM(posInd,posInd)));
%FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
%TN = sum(sum(CM(negInd,negInd)));
if (TP+FP) ~= 0
    v = TP/(TP+FP);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateNegativePredictiveValue(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
negInd = h.ControlClasses;
%TP = sum(sum(CM(posInd,posInd)));
FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
%FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
TN = sum(sum(CM(negInd,negInd)));
if (FN+TN) ~= 0
    v = TN/(FN+TN);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculatePositiveLikelihood(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
negInd = h.ControlClasses;
TP = sum(sum(CM(posInd,posInd)));
FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
TN = sum(sum(CM(negInd,negInd)));
if (TP+FN)*FP ~= 0
    v = TP*(FP+TN)/(TP+FN)/FP;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateNegativeLikelihood(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
negInd = h.ControlClasses;
TP = sum(sum(CM(posInd,posInd)));
FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
TN = sum(sum(CM(negInd,negInd)));
if (TP+FN)*TN ~= 0
    v = FN*(FP+TN)/(TP+FN)/TN;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculatePrevalence(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
negInd = h.ControlClasses;
TP = sum(sum(CM(posInd,posInd)));
FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
TN = sum(sum(CM(negInd,negInd)));
if (TP+FN+FP+TN) ~= 0
    v = (TP+FN)/(TP+FN+FP+TN);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = calculateDiagnosticTable(h,w) %#ok
v = NaN;
if ~h.ListenerEnable
    return
end
CM = h.CountingMatrix;
numClasses = size(CM,2);
posInd = h.TargetClasses;
negInd = h.ControlClasses;
TP = sum(sum(CM(posInd,posInd)));
FN = sum(sum(CM(setdiff(1:numClasses+1,posInd),posInd)));
FP = sum(sum(CM(setdiff(1:numClasses+1,negInd),negInd)));
TN = sum(sum(CM(negInd,negInd)));
v = [TP FP;FN TN];
