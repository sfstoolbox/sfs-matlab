function [node_indeces,weights] = findconvexcone(x0,xs)
%FINDCONVEXCONE finds three points that span a cone
%
%   Usage: [node_indeces,weights] = findconvexcone(x0,xs);
%
%   Input parameters:
%       x0                    - points in R3 (N x 3) 
%       xs                    - a point in R3 (1 x 3)
%
%   output parameters:
%       node_indeces          - row indeces of 3 points in x0 (3 x 1)
%       weights               - weights (3 x 1)
%
%   FINDCONVEXCONE(x0,xs) returns node_indeces and non-negative weights 
%   [w1;w2;w3] such that
%      
%   xs = w1*x1 + w2*x2 + w3*x3 , 
%   where [x1; x2; x3] = x0(node_indeces,:) .
%
%
%   See also: findnearestneighbor
%
%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Computation =====================================================

% find x0 with smallest angle to xs
xs_norm = repmat(xs./norm(xs,2),size(x0,1),1);
x0_norm = x0./repmat(vector_norm(x0,2),[1,3]);
[~,closest_node] = max(vector_product(x0_norm,xs_norm,2));

% from the triangulation of the convex hull:
% get all triangles at "closest_node"
triangles = convhull(x0);
mask = logical(sum(triangles==closest_node,2));
triangles = triangles(mask,:);
assert(~isempty(triangles),'something went wrong: x0 non-convex?')

% one of these cones contains xs
linear_factors = zeros(size(triangles,1),3);
for n = 1:size(triangles,1);
    A = x0(triangles(n,:),:);
    linear_factors(n,:) = xs/A;
end
% there's one conic combination (i.e. non-negative coefficients)
mask = logical(~sum(linear_factors < 0,2));
if ~sum(mask) 
    error('xs is not in convex cone of x0.');
end
assert(sum(mask) == 1, 'Should not happen!');

weights = linear_factors(mask,:).';
node_indeces = triangles(mask,:).';

[weights,order] = sort(weights,'descend');
node_indeces = node_indeces(order);

