function hSV = svmplotsvs(hAxis,hLines,groupString,svm_struct)
%SVMPLOTSVS helper function for SVM plotting

%   Copyright 2004-2008 The MathWorks, Inc.
%   $Revision: 1.1.12.6 $  $Date: 2008/06/16 16:32:53 $

% plot the support vectors and the separating line
hold on;
sv = svm_struct.SupportVectors;
% see if we need to unscale the data
scaleData = svm_struct.ScaleData;
if ~isempty(scaleData)
    for c = 1:size(sv, 2)
        sv(:,c) = (sv(:,c)./scaleData.scaleFactor(c)) - scaleData.shift(c);
    end
end
% plot support vectors
hSV = plot(sv(:,1),sv(:,2),'ko');

lims = axis(hAxis);
[X,Y] = meshgrid(linspace(lims(1),lims(2)),linspace(lims(3),lims(4)));
Xorig = X; Yorig = Y;

% need to scale the mesh
if ~isempty(scaleData)
    X = scaleData.scaleFactor(1) * (X + scaleData.shift(1));
    Y = scaleData.scaleFactor(2) * (Y + scaleData.shift(2));
end

[dummy, Z] = svmdecision([X(:),Y(:)],svm_struct); 

contour(Xorig,Yorig,reshape(Z,size(X)),[0 0],'k');
hold off;
labelString = cellstr(groupString);
labelString{end+1} = 'Support Vectors';
legend([hLines(1),hLines(2),hSV],...
    labelString);