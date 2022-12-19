function [resTab, allocPossible] = CheckResourcesAllocation(firstRB, lastRB, firstTS, lastTS, resTab, userIdx, timeSlotsNum, resBlocksNum)
%CheckResourcesAllocation - checks if resources can be allocated in
%specific time and frequency
%   Parameters:
%       firstRB - first index of frequency interval
%       lastRB - last index of frequency interval
%       firstTS - first index of time interval
%       lastTS - last index of time interval
%       resTab - table of frequency and time allocated per user
%   Output:
%       resTab - updated table of frequency and time allocated per user
%       allocPossible - determines if allocation is possible or not
%       (boolean value)

% disp('    resources allocation check')
allocPossible = true;

if lastRB > resBlocksNum
    allocPossible = false;
    return;
end

if lastTS > timeSlotsNum
    allocPossible = false;
    return;
end

[rows, cols] = size(resTab);

if cols ~= 4
    error('wrong resTab dimensions')
end

for r = 1:rows
    if (firstRB >= resTab(r, 1)) && (firstRB <= resTab(r, 2))
        if (firstTS >= resTab(r, 3)) && (firstTS <= resTab(r, 4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
        if (lastTS >= resTab(r, 3)) && (lastTS <= resTab(r, 4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
        if(firstTS <= resTab(r,3)) && (lastTS >= resTab(r,4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
    end
    if (lastRB >= resTab(r, 1)) && (lastRB <= resTab(r, 2))
        if (firstTS >= resTab(r, 3)) && (firstTS <= resTab(r, 4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
        if (lastTS >= resTab(r, 3)) && (lastTS <= resTab(r, 4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
        if(firstTS <= resTab(r,3)) && (lastTS >= resTab(r,4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
    end
    if (firstRB <= resTab(r, 1)) && (lastRB >= resTab(r, 2))
        if (firstTS >= resTab(r, 3)) && (firstTS <= resTab(r, 4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
        if (lastTS >= resTab(r, 3)) && (lastTS <= resTab(r, 4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
        if(firstTS <= resTab(r,3)) && (lastTS >= resTab(r,4))
            allocPossible = false;
%             disp(['    resources allocation impossible' newline])
            return;
        end
    end
end

if allocPossible == true
%     disp(['    allocating resources' newline])
    resTab(userIdx, 1) = firstRB;
    resTab(userIdx, 2) = lastRB;
    resTab(userIdx, 3) = firstTS;
    resTab(userIdx, 4) = lastTS;
end

end

