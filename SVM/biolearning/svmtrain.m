function [svm_struct, svIndex] = svmtrain(training, groupnames, varargin)
%SVMTRAIN trains a support vector machine classifier
%
%   SVMStruct = SVMTRAIN(TRAINING,GROUP) trains a support vector machine
%   classifier using data TRAINING taken from two groups given by GROUP.
%   SVMStruct contains information about the trained classifier, including
%   the support vectors, that is used by SVMCLASSIFY for classification.
%   GROUP is a column vector of values of the same length as TRAINING that
%   defines two groups. Each element of GROUP specifies the group the
%   corresponding row of TRAINING belongs to. GROUP can be a numeric
%   vector, a string array, or a cell array of strings. SVMTRAIN treats
%   NaNs or empty strings in GROUP as missing values and ignores the
%   corresponding rows of TRAINING.
%
%   SVMTRAIN(...,'KERNEL_FUNCTION',KFUN) allows you to specify the kernel
%   function KFUN used to map the training data into kernel space. The
%   default kernel function is the dot product. KFUN can be one of the
%   following strings or a function handle:
%
%       'linear'      Linear kernel or dot product
%       'quadratic'   Quadratic kernel
%       'polynomial'  Polynomial kernel (default order 3)
%       'rbf'         Gaussian Radial Basis Function kernel
%       'mlp'         Multilayer Perceptron kernel (default scale 1)
%       function      A kernel function specified using @,
%                     for example @KFUN, or an anonymous function
%
%   A kernel function must be of the form
%
%         function K = KFUN(U, V)
%
%   The returned value, K, is a matrix of size M-by-N, where U and V have M
%   and N rows respectively.  If KFUN is parameterized, you can use
%   anonymous functions to capture the problem-dependent parameters. For
%   example, suppose that your kernel function is
%
%       function k = kfun(u,v,p1,p2)
%       k = tanh(p1*(u*v')+p2);
%
%   You can set values for p1 and p2 and then use an anonymous function:
%       @(u,v) kfun(u,v,p1,p2).
%
%   SVMTRAIN(...,'RBF_SIGMA',SIGMA) allows you to specify the scaling
%   factor, sigma, in the radial basis function kernel.
%
%   SVMTRAIN(...,'POLYORDER',ORDER) allows you to specify the order of a
%   polynomial kernel. The default order is 3.
%
%   SVMTRAIN(...,'MLP_PARAMS',[P1 P2]) allows you to specify the
%   parameters of the Multilayer Perceptron (mlp) kernel. The mlp kernel
%   requires two parameters, P1 and P2, where K = tanh(P1*U*V' + P2) and P1
%   > 0 and P2 < 0. Default values are P1 = 1 and P2 = -1.
%
%   SVMTRAIN(...,'METHOD',METHOD) allows you to specify the method used
%   to find the separating hyperplane. Options are
%
%       'QP'  Use quadratic programming (requires the Optimization Toolbox)
%       'SMO' Use Sequential Minimal Optimization method
%       'LS'  Use least-squares method
%
%   If you have the Optimization Toolbox, then the QP method is the default
%   method. If not, the default method is SMO. When using the QP method,
%   the classifier is a 2-norm soft-margin support vector machine.
%
%   SVMTRAIN(...,'QUADPROG_OPTS',OPTIONS) allows you to pass an OPTIONS
%   structure created using OPTIMSET to the QUADPROG function when using
%   the 'QP' method. See help optimset for more details.
%
%   SVMTRAIN(...,'SMO_OPTS',SMO_OPTIONS) allows you to set options for the
%   'SMO' method. SMO_OPTIONS should be created using the function
%   SVMSMOSET.
%
%   SVMTRAIN(...,'BOXCONSTRAINT',C) allows you to set the box constraint C
%   for the soft margin. The default value is 1.  C can be a scalar or a
%   vector of the same length as the training data. Note that in older
%   versions of Bioinformatics Toolbox the default value for C was
%   1/sqrt(eps) which will only classify separable data.
%
%   SVMTRAIN(...,'AUTOSCALE', AUTOSCALEVAL) allows you to specify whether
%   or not to automatically shift and scale the data points before
%   training. Default is true.
%
%   SVMTRAIN(...,'SHOWPLOT',true), when used with two-dimensional data,
%   creates a plot of the grouped data and plots the separating line for
%   the classifier.
%
%   Example:
%       % Load the data and select features for classification
%       load fisheriris
%       data = [meas(:,1), meas(:,2)];
%       % Extract the Setosa class
%       groups = ismember(species,'setosa');
%       % Randomly select training and test sets
%       [train, test] = crossvalind('holdOut',groups);
%       cp = classperf(groups);
%       % Use a linear support vector machine classifier
%       svmStruct = svmtrain(data(train,:),groups(train),'showplot',true);
%       % Add a title to the plot
%       title(sprintf('Kernel Function: %s',...
%             func2str(svmStruct.KernelFunction)),...
%             'interpreter','none');
%       % Classify the test set using svmclassify
%       classes = svmclassify(svmStruct,data(test,:),'showplot',true);
%       % See how well the classifier performed
%       classperf(cp,classes,test);
%       cp.CorrectRate
%
%   See also CLASSIFY, CLASSPERF, CROSSVALIND, KNNCLASSIFY, QUADPROG, 
%   SVMCLASSIFY, SVMSMOSET.

%   Copyright 2004-2008 The MathWorks, Inc.
%   $Revision: 1.1.12.13 $  $Date: 2008/06/16 16:32:47 $

%   References:
%
%     [1] Cristianini, N., Shawe-Taylor, J An Introduction to Support
%         Vector Machines, Cambridge University Press, Cambridge, UK. 2000.
%         http://www.support-vector.net
%     [2] Kecman, V, Learning and Soft Computing,
%         MIT Press, Cambridge, MA. 2001.
%     [3] Suykens, J.A.K., Van Gestel, T., De Brabanter, J., De Moor, B.,
%         Vandewalle, J., Least Squares Support Vector Machines,
%         World Scientific, Singapore, 2002.
%     [4] J.C. Platt: A Fast Algorithm for Training  Support Vector
%         Machines,  Advances in Kernel Methods - Support Vector Learning,
%         B. Sch?lkopf, C. Burges, and A. Smola, eds., MIT Press, 1998. 
%     [5] J.C. Platt: Fast Training of Support Vector Machines using
%         Sequential Minimal Optimization Microsoft Research Technical
%         Report MSR-TR-98-14, 1998.
%     [6] http://www.kernel-machines.org/papers/tr-30-1998.ps.gz
%
%   SVMTRAIN(...,'KFUNARGS',ARGS) allows you to pass additional
%   arguments to kernel functions.


% check inputs
bioinfochecknargin(nargin,2,mfilename)

% set defaults


plotflag = false;
% The large scale solver cannot handle this type of problem, so turn it
% off.

qp_opts = optimset('LargeScale','Off','display','off');
smo_opts = svmsmoset;
kfunargs = {};
setPoly = false; usePoly = false;
setMLP = false; useMLP = false;
setSigma = false; useSigma = false;

autoScale = true;
% default optimization method
if ~isempty(which('quadprog'))
    optimMethod = 'QP';
else
    optimMethod = 'SMO';
end

% set default kernel function
kfun = @linear_kernel;




numoptargs = nargin -2;
optargs = varargin;

% check group is a vector -- though char input is special...
if ~isvector(groupnames) && ~ischar(groupnames)
    error('Bioinfo:svmtrain:GroupNotVector',...
        'Group must be a vector.');
end

% grp2idx sorts a numeric grouping var ascending, and a string grouping
% var by order of first occurrence

[groupIndex, groupString] = grp2idx(groupnames);

% make sure that the data are correctly oriented.
if size(groupnames,1) == 1
    groupnames = groupnames';
end
% make sure data is the right size
if size(training,1) ~= length(groupnames)
    if size(training,2) == length(groupnames)
        training = training';
    else
        error('Bioinfo:svmtrain:DataGroupSizeMismatch',...
            'GROUP and TRAINING must have the same number of rows.')
    end
end

% check for NaN in data matrix:
if any(isnan(training(:)))
    error('Bioinfo:svmtrain:NaNinDataMatrix', ...
        'TRAINING data must not contain missing values');
end

% NaNs are treated as unknown classes and are removed from the training
% data
nans = isnan(groupIndex);
if any(nans)
    training(nans,:) = [];
    groupIndex(nans) = [];
end
ngroups = length(groupString);
nPoints = length(groupIndex);
% set default value of box constraint

boxconstraint = ones(nPoints, 1);

if ngroups > 2
    error('Bioinfo:svmtrain:TooManyGroups',...
        'SVMTRAIN only supports classification into two groups.\nGROUP contains %d different groups.',ngroups)
end
% convert to 1, -1.
groupIndex = 1 - (2* (groupIndex-1));

% handle optional arguments

if  numoptargs >= 1
    if rem(numoptargs,2)== 1
        error('Bioinfo:svmtrain:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'kernel_function','method','showplot','kfunargs',...
        'quadprog_opts','polyorder','mlp_params',...
        'boxconstraint','rbf_sigma','autoscale', 'smo_opts'};
    for j=1:2:numoptargs
        pname = optargs{j};
        pval = optargs{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:svmtrain:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:svmtrain:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % kernel_function
                    if ischar(pval)
                        okfuns = {'linear','quadratic',...
                            'radial','rbf','polynomial','mlp'};
                        funNum = strmatch(lower(pval), okfuns);
                        if isempty(funNum)
                            funNum = 0;
                        end
                        switch funNum  %maybe make this less strict in the future
                            case 1
                                kfun = @linear_kernel;
                            case 2
                                kfun = @quadratic_kernel;
                            case {3,4}
                                kfun = @rbf_kernel;
                                useSigma = true;
                            case 5
                                kfun = @poly_kernel;
                                usePoly = true;
                            case 6
                                kfun = @mlp_kernel;
                                useMLP = true;
                            otherwise
                                error('Bioinfo:svmtrain:UnknownKernelFunction',...
                                    'Unknown Kernel Function %s.',pval);
                        end
                    elseif isa (pval, 'function_handle')
                        kfun = pval;
                    else
                        error('Bioinfo:svmtrain:BadKernelFunction',...
                            'The kernel function input does not appear to be a function handle\nor valid function name.')
                    end
                case 2  % method
                    if strncmpi(pval,'qp',2)
                        if isempty(which('quadprog'))
                            warning('Bioinfo:svmtrain:NoOptim',...
                                'The Optimization Toolbox is required to use the quadratic programming method.')
                            optimMethod = 'SMO';
                        else
                            optimMethod = 'QP';
                        end
                    elseif strncmpi(pval,'ls',numel(pval))
                        optimMethod = 'LS';
                    elseif strncmpi(pval, 'smo', numel(pval))
                        optimMethod = 'SMO';
                    else
                        error('Bioinfo:svmtrain:UnknownMethod',...
                            'Unknown method option %s. Valid methods are ''SMO'', ''QP'' and ''LS''',pval);

                    end
                case 3  % display
                    plotflag = opttf(pval,okargs{k},mfilename);
                    if plotflag == true
                        if size(training,2) == 2
                            plotflag = true;
                        else
                            plotflag = false;
                            warning('Bioinfo:svmtrain:OnlyPlot2D',...
                                'The display option can only plot 2D training data.')
                        end

                    end
                case 4 % kfunargs
                    if iscell(pval)
                        kfunargs = pval;
                    else
                        kfunargs = {pval};
                    end
                case 5 % quadprog_opts
                    if isstruct(pval)
                        qp_opts = optimset(qp_opts,pval);
                    elseif iscell(pval)
                        qp_opts = optimset(qp_opts,pval{:});
                    else
                        error('Bioinfo:svmtrain:BadQuadprogOpts',...
                            'QUADPROG_OPTS must be an opts structure.');
                    end
                case 6 % polyorder
                    if ~isscalar(pval) || ~isnumeric(pval)
                        error('Bioinfo:svmtrain:BadPolyOrder',...
                            'POLYORDER must be a scalar value.');
                    end
                    if pval ~=floor(pval) || pval < 1
                        error('Bioinfo:svmtrain:PolyOrderNotInt',...
                            'The order of the polynomial kernel must be a positive integer.')
                    end
                    kfunargs = {pval};
                    setPoly = true;

                case 7 % mlpparams
                    if ~isnumeric(pval)||numel(pval)~=2
                        error('Bioinfo:svmtrain:BadMLPParams',...
                            'MLP_PARAMS must be a two element numeric array.');
                    end
                    if pval(1) <= 0
                        error('Bioinfo:svmtrain:MLPWeightNotPositive',...
                            'The weight for the MLP kernel should be positive.')
                    end
                    if pval(2) > 0
                        warning('Bioinfo:svmtrain:MLPBiasNotNegative',...
                            'The bias for MLP kernel should be negative.')
                    end
                    kfunargs = {pval(1),pval(2)};
                    setMLP = true;
                case 8  % box constraint: it can be a positive numeric scalar
                    % or a numeric vector of the same length as there are
                    % data points
                    if isscalar(pval) && isnumeric(pval) && pval > 0
                        % scalar input: adjust to group size and transform into vector
                        n1 = length(find(groupIndex==1));
                        n2 = length(find(groupIndex==-1));
                        c1 = 0.5 * pval * nPoints / n1;
                        c2 = 0.5 * pval * nPoints / n2;
                        boxconstraint(groupIndex==1) = c1;
                        boxconstraint(groupIndex==-1) = c2;
                    elseif isvector(pval) && isnumeric(pval) && all(pval > 0)
                        % vector input
                        if length(pval) ~= nPoints
                            error('Bioinfo:svmtrain:BoxConstraintNotScalar',...
                                'If box constraint is passed as vector, its size must equal the number of training points');
                        end
                        boxconstraint = pval;
                    else
                        error('Bioinfo:svmtrain:BoxConstraintNotScalar',...
                            'The box constraint must be a numeric scalar or vector > 0.');
                    end
                    % If boxconstraint == Inf then convergence will not
                    % happen so fix the value to 1/sqrt(eps).
                    boxconstraint = min(boxconstraint,repmat(1/sqrt(eps(class(boxconstraint))),...
                        size(boxconstraint)));
                case 9  % rbf sigma
                    if isscalar(pval) && isnumeric(pval)
                        kfunargs = {pval};
                        setSigma = true;
                    else
                        error('Bioinfo:svmtrain:RBFSigmaNotScalar',...
                            'Sigma must be a numeric scalar.');
                    end
                case 10 % autoscale
                    autoScale = opttf(pval,okargs{k},mfilename);
                case 11  % smo_opts
                    smo_opts = pval;
            end
        end
    end
end
if setPoly && ~usePoly
    warning('Bioinfo:svmtrain:PolyOrderNotPolyKernel',...
        'You specified a polynomial order but not a polynomial kernel');
end
if setMLP && ~useMLP
    warning('Bioinfo:svmtrain:MLPParamNotMLPKernel',...
        'You specified MLP parameters but not an MLP kernel');
end
if setSigma && ~useSigma
    warning('Bioinfo:svmtrain:RBFParamNotRBFKernel',...
        'You specified radial basis function parameters but not a radial basis function kernel');
end


% plot the data if requested
if plotflag
    [hAxis,hLines] = svmplotdata(training,groupIndex);
    legend(hLines,cellstr(groupString));
end

% autoscale data if required, we can't use the zscore function here,
% because we need the shift and scale values.
scaleData = [];
if autoScale
    scaleData.shift = - mean(training);
    stdVals = std(training);
    scaleData.scaleFactor = 1./stdVals;
    % leave zero-variance data unscaled:
    scaleData.scaleFactor(~isfinite(scaleData.scaleFactor)) = 1;

    % shift and scale columns of data matrix:
    for c = 1:size(training, 2)
        training(:,c) = scaleData.scaleFactor(c) * ...
            (training(:,c) +  scaleData.shift(c));
    end
end


if strcmpi(optimMethod, 'SMO')

    % if we have a kernel that takes extra arguments we must define a new
    % kernel function handle to be passed to seqminopt
    if ~isempty(kfunargs)
        tmp_kfun = @(x,y) feval(kfun, x,y, kfunargs{:});
    else
        tmp_kfun = kfun;
    end

    [alpha bias] = seqminopt(training, groupIndex, ...
        boxconstraint, tmp_kfun, smo_opts);

    svIndex = find(alpha > sqrt(eps));
    sv = training(svIndex,:);
    alphaHat = groupIndex(svIndex).*alpha(svIndex);

else % QP and LS both need the kernel matrix:

    % calculate kernel function and add additional term required
    % for two-norm soft margin
    try
        kx = feval(kfun,training,training,kfunargs{:});
        % ensure function is symmetric
        kx = (kx+kx')/2 + diag(1./boxconstraint);
    catch theException
        error('Bioinfo:svmtrain:KernelFunctionError',...
            'Error calculating the kernel function:\n%s\n', theException.message);
    end

    % create Hessian
    H =((groupIndex * groupIndex').*kx);

    if strcmpi(optimMethod, 'QP')

        % X=QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,opts)
        [alpha, fval, exitflag, ...
            output, lambda] = quadprog(H,-ones(nPoints,1),[],[],...
            groupIndex',0,zeros(nPoints,1), ...
            Inf *ones(nPoints,1), ones(nPoints,1), ...
            qp_opts); %#ok

        if exitflag <= 0
            error('Bioinfo:svmtrain:UnsolvableOptimization',...
                'Unable to solve the optimization problem:\n%s\n', output.message);
        end

        % The support vectors are the non-zeros of alpha.
        % We could also use the zero values of the Lagrangian (fifth output of
        % quadprog) though the method below seems to be good enough.
        svIndex = find(alpha > sqrt(eps));
        sv = training(svIndex,:);

        % calculate the parameters of the separating line from the support
        % vectors.
        alphaHat = groupIndex(svIndex).*alpha(svIndex);

        % Calculate the bias by applying the indicator function to the support
        % vector with largest alpha.
        [maxAlpha,maxPos] = max(alpha); 
        bias = groupIndex(maxPos) - sum(alphaHat.*kx(svIndex,maxPos));
        % an alternative method is to average the values over all support vectors
        % bias = mean(groupIndex(sv)' - sum(alphaHat(:,ones(1,numSVs)).*kx(sv,sv)));

        % An alternative way to calculate support vectors is to look for zeros of
        % the Lagrangians (fifth output from QUADPROG).
        %
        % [alpha,fval,output,exitflag,t] = quadprog(H,-ones(nPoints,1),[],[],...
        %             groupIndex',0,zeros(nPoints,1),inf *ones(nPoints,1),zeros(nPoints,1),opts);
        %
        % sv = t.lower < sqrt(eps) & t.upper < sqrt(eps);
    else  % Least-Squares
        % now build up compound matrix for solver
        A = [0 groupIndex';groupIndex,H];
        b = [0;ones(size(groupIndex))];
        x = A\b;

        % calculate the parameters of the separating line from the support
        % vectors.
        sv = training;
        bias = x(1);
        alphaHat = groupIndex.*x(2:end);
        svIndex = (1:nPoints)';
    end
end
svm_struct.SupportVectors = sv;
svm_struct.Alpha = alphaHat;
svm_struct.Bias = bias;
svm_struct.KernelFunction = kfun;
svm_struct.KernelFunctionArgs = kfunargs;
svm_struct.GroupNames = groupnames;
svm_struct.SupportVectorIndices = svIndex;
svm_struct.ScaleData = scaleData;
svm_struct.FigureHandles = [];
if plotflag
    hSV = svmplotsvs(hAxis,hLines,groupString,svm_struct);
    svm_struct.FigureHandles = {hAxis,hLines,hSV};
end



