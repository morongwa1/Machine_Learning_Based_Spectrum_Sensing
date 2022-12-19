function [detectProc, falseAlarmProc] = detectionCheck(resourceAlloc, detectedRes)
%DETECTIONCHECK chcecks what percent of transmitted data was correctly
%detected
%   Parameters:
%       resourceAlloc - resource allocation table for comparing
%       detectedRes - detected resource by sensing algorithm
%   Output:
%       detectProc - percent of correctly detected resources
%       falseAlarmProc - percent of resoruces that were detected but not transmitted 

[row,col] = size(resourceAlloc);
falseAlarm = 0;
detect = 0;

RBocc = sum(sum(resourceAlloc));

for m = 1:row
    for n = 1:col
        if resourceAlloc(m,n) == 1
            if detectedRes(m,n) == 1
                detect = detect + 1;
            end
        else
            % resourceAlloc(m,n) == 0
            if detectedRes(m,n) == 1
                falseAlarm = falseAlarm + 1;
            end
        end
    end
end
%disp(['Total Detected block ', num2str(detect)]);
%disp(['Total resource blocks ', num2str(RBocc)]);
detectProc = detect*100/RBocc;
falseAlarmProc = falseAlarm*100/(row*col - RBocc);
end

