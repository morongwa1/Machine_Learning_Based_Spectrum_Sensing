function [vector] = meanRepeatingValuesInVector(vector)

[~,idx] = sort(vector(:,1));
vector = vector(idx,:);
[idu, ia, ic]=unique(vector(:,1));
vector2 = [vector(ia,:), ones(length(idu),1)];

if length(idu) == length(vector(:,1))
    disp('vector values are not repeated')
else
    disp('some vector values are repeated')
    disp(['number of repeated values = ' num2str(length(vector(:,1))-length(idu))])
    for i = 2:length(vector(:,1))
        if ic(i) == ic(i-1)
            vector2(ic(i),2) = vector2(ic(i),2) + vector(i,2);
            vector2(ic(i),3) = vector2(ic(i),3) + vector(i,3);
            vector2(ic(i),4) = vector2(ic(i),4)+1;
        end
    end
end

vector2(:,2) = vector2(:,2)./vector2(:,4);
vector2(:,3) = vector2(:,3)./vector2(:,4);
vector = vector2;
end