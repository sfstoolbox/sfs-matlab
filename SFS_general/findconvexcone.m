function [x0_indeces,weights] = findconvexcone(x0,xs)
%FINDCONVEXCONE finds three points from x0 with xs in their conic span
%
%   Usage: [x0_indeces,weights] = findconvexcone(x0,xs);
%
%   Input parameters:
%       x0                    - point cloud in R3 (N x 3)
%       xs                    - point in R3 (1 x 3)
%
%   output parameters:
%       x0_indeces            - row indeces of 3 points in x0 (3 x 1)
%       weights               - weights (3 x 1)
%
%   FINDCONVEXCONE(x0,xs) returns three row indeces into x0 and
%   non-negative weights [w1;w2;w3] such that xs lies in the convex cone
%   with minimum solid angle.
%       xs = w1*x1 + w2*x2 + w3*x3 , 
%       where [x1; x2; x3] = x0(x0_indeces,:) .
%     
%   (If all x0 and xs have unit norm this is VBAP.)  
%   
%   This may fail when 
%     a) x0 is not convex, or
%     b) The convex hull of x0 does not contain the origin.
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
xs_normed = repmat(xs./norm(xs,2),size(x0,1),1);
x0_normed = x0./repmat(vector_norm(x0,2),[1,3]);
[~,most_aligned_point] = max(vector_product(x0_normed,xs_normed,2));

% Delaunay triangulation of convex hull
triangles = convhull(x0);

% get all triangles at "most aligned point"
mask = logical(sum(triangles==most_aligned_point,2));
triangles = triangles(mask,:);
if isempty(triangles)
    error('x0 contains a point in the interior of its convex hull');
end

% one of the triangles span a convex cone that contains xs
for n = 1:size(triangles,1);
    A = x0(triangles(n,:),:);
    weights = xs/A;
    if ~sum(weights < 0) % non-negative weights == conic combination
        x0_indeces = triangles(n,:);
        break;
    end
end
if sum(weights < 0)
    error('xs is not in convex cone of x0.');
end

[weights,order] = sort(weights.','descend');
x0_indeces = x0_indeces(order).';

