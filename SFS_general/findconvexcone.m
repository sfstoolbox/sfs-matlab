function [idx,weights] = findconvexcone(x0,xs)
%FINDCONVEXCONE selects up to 3 points from x0 with xs in their conic span
%
%   Usage: [idx,weights] = findconvexcone(x0,xs)
%
%   Input parameters:
%       x0          - point cloud on a sphere around the origin / m [nx3]
%       xs          - desired direction as point in R^3 / m [1x3]
%
%   Output parameters:
%       idx         - row indices of N points in x0 [Nx1]
%                     where N is 1,2 or 3
%       weights     - weights [Nx1]
%
%   FINDCONVEXCONE(x0,xs) returns 1,2 or 3 row indices into x0 and non-negative
%   weights w1, ..., w3 such that w1*x1 + w2*x2 + w3*x3 with
%   [x1; x2; x3] == x0(idx,:) composes the point inside the triangle spanned
%   by x1, x2, x3.
%
%   x1...x3 are selected from the convex hull in R3.
%   Various precautions are taken to make this well-behaved in most cases.
%
%   (If all x0 and xs have unit norm this is VBAP.)
%
%   See also: findnearestneighbour, test_interpolation_point_idx

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
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Prepare Grid (see local functions below) ========================
% Normalise x0 and xs to lie on unit sphere as only direction is relevant
xs = xs./norm(xs,2);
radii = vector_norm(x0,2);
if abs(max(radii) - min(radii)) >1e-3
     warning('%s: Grid is apparently not a sphere.', upper(mfilename))
end
x0 = x0./repmat(radii,[1,size(x0,2)]);

% Rotate to principal axes to enable 2D arrays
[x0, xs] = rotate_to_principal_axes(x0, xs);

% Calculate dummy points to enable "partial" arrays
dummy_points = augment_bounding_box(x0);
if ~isempty(dummy_points)
    dummy_indices = (1:size(dummy_points,1)) + size(x0,1);
    x0 = [x0; dummy_points];
end


%% ===== Computation =====================================================
% Delaunay triangulation of convex hull
simplices = convhulln(x0);

% Find x0 with smallest angle to xs
[~,most_aligned_point] = ...
    max(vector_product(x0,repmat(xs,size(x0,1),1),2));

% The simplices at "most aligned point" are the most likely candidates,
% put them at the beginning of the list
mask = logical(sum(simplices==most_aligned_point,2));
simplices = [simplices(mask,:); simplices(~mask,:) ];

% One of these simplices contains xs
for n = 1:size(simplices,1);
    A = x0(simplices(n,:),:);
    weights = xs/A;
    weights(abs(weights)<1e-10) = 0;
    if all(weights >= 0) % non-negative weights == convex combination
        idx = simplices(n,:);
        break;
    end
end
assert(all(weights >= 0), '%s: Negative weights. Shall never happen.', ...
    upper(mfilename))

if ~isempty(dummy_points)
    % Remove possible dummies from selected points
    dummy_mask = (dummy_indices == idx);
    if any(dummy_mask)
        idx(dummy_mask) = [];
        if ~weights(dummy_mask)==0
            warning('%s: Requested point lies outside grid.', upper(mfilename))
        end
        weights(dummy_mask) = [];
    end
end

% Normalise weights
weights = weights/sum(weights);

[weights,order] = sort(weights.','descend');
idx = idx(order).';
end

% =========================================================================

function [x0, xs] = rotate_to_principal_axes(x0, xs, gamma)
%ROTATE_TO_PRINCIPAL_AXES rotates x0 and xs to x0's principal axes.
% If the ratio of second-to-smallest singular values is < gamma,
% the third dimension is discarded.
%
%   Input parameters:
%       x0          - point cloud in R^3
%       xs          - point in R^3
%       gamma       - scalar in 0 < gamma << 1
%
%   Output parameters:
%       x0          - point cloud in R^3 or R^2
%       xs          - point in R^3 or R^2
if nargin < 3
    gamma = 0.1; % inverse of aspect ratio of principal axes
end

[~,S,V] = svd(x0);
x0 = x0*V;
xs = xs*V;
S = diag(S);
if S(end)/S(end-1) < gamma
    x0(:,3) = [];
    xs(3) = [];
    warning('SFS:findconvexcone','%s: Grid is apparently two-dimensional. ', ...
        upper(mfilename));
end
end

% =========================================================================

function dummy_points = augment_bounding_box(x0)
%AUGMENT_BOUNDING_BOX yields dummy points such that the origin is
% contained in the cartesian bounding box of x0.
dummy_points = -diag(sign(max(x0)) + sign(min(x0)));
dummy_points(~any(dummy_points,2),:) = [];
end
