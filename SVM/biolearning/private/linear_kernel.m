function K = linear_kernel(u,v,varargin) 
%LINEAR_KERNEL linear kernel for SVM functions

% Copyright 2004-2008 The MathWorks, Inc.

K = (u*v');
