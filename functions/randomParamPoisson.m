function [param] = randomParamPoisson(mu, sigma, integer, min, max)
%randomParam generates random parameter according to gaussian
%distribution
%   Parameters:
%       mu - mean
%       sigma - standard deviation
%       integer - true if returned value must be an integer value, false
%         otherwise
%       min - minimum value that can be generated (optional)
%       max - maximum value that can be generated (optional)
%   Output:
%       param - generated random value

rng('shuffle')

if nargin < 3
    error('too few arguments')
end


param = random('Poisson', mu);

if (nargin>=4)
    while (param > max) || (param < min)
        param = random('Poisson', mu);
    end
end
    
if (nargin >= 3) && (integer == true)
    param = round(param);
end

end

