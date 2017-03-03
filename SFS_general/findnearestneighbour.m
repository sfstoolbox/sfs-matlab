function [idx,weights] = findnearestneighbour(A,b,N)
%FINDNEARESTNEIGHBOUR finds the N nearest neighbours and weights (for N<=2)
%
%   Usage: [idx,weights] = findnearestneighbour(A,b,N)
%
%   Input parameters:
%       A        - points in R^3 / m [nx3]
%       b        - desired point in R^3 / m [1x3]
%       N        - number of nearest neighbours to find
%
%   Output parameters:
%       idx      - row indices of N points in A [Nx1]
%       weights  - weights [1x1] or [2x1]
%
%   FINDNEARESTNEIGHBOUR(A,b,N) returns N indices idx for the nearest neighbour
%   points from A to point b. For N<=2, weights for linear interpolation are also
%   returned.
%
%   See also: find, findrows, get_ir

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    N = 1;
end
if nargout==2 && N>=3
    error(['%s: calculation of weights for 3 or more nearest neighbours ', ...
        'is not implemented.'],upper(mfilename));
end
% Ensure row vector
if size(b,1)>1
    b=b.';
end


%% ===== Computation =====================================================
% Calculate distance between points
distance = vector_norm(bsxfun(@minus,A.',b.'),1);
% Sort the distances in order to find the n lowest once
[~,idx] = sort(distance);
idx = idx(1:min(N,length(idx))).';

% Determine weights for linear 1D interpolation over angle
cos_alpha = sum(bsxfun(@times,A(idx,:).',b.'))./vector_norm(A(idx,:),2).'/norm(b);
% Ensure range -1...1
cos_alpha = min(cos_alpha,1);
cos_alpha = max(cos_alpha,-1);
angles = acos(cos_alpha);
if angles == 0
    weights = 1;
else
    angles = fliplr(angles);
    weights = angles.'/sum(angles);
end
assert(all(weights >= 0), '%s: Negative weights. Shall never happen.', ...
    upper(mfilename))
