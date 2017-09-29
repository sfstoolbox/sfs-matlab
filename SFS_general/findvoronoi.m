function [idx,weights] = findvoronoi(x0,xs)
%FINDVORONOI finds the corresponding voronoi regions to the points x0
% surrounding the desired point xs. The weights are derived from the voronoi
% region surface area differences with and without xs.
%
%   Usage: [idx,weights] = findvoronoi(x0,xs)
%
%   Input parameters:
%       x0          - point cloud on a sphere around the origin / m [nx3]
%       xs          - desired point in R^3 / m [1x3]
%
%   Output parameters:
%       idx         - row indices of N points in x0 [Nx1]
%       weights     - weights [Nx1]
%
%   
%
%   See also: findnearestneighbour, findconvexcone,
%             test_interpolation_point_selection

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

%% ===== Prepare Grid ====================================================
% Normalise x0 and xs to lie on unit sphere
xs = xs./norm(xs,2);
radii = vector_norm(x0,2);
if abs(max(radii) - min(radii)) >1e-3
     warning('%s: Grid is apparently not a sphere.', upper(mfilename))
end
x0 = x0./repmat(radii,[1,size(x0,2)]);

% Center of the Sphere
center = [0 0 0];

% In case xs is colinear with or equals 1 point there is no interpolation needed
eq_idx = findrows(x0, xs);
if ~isempty(eq_idx) && size(eq_idx,1)<2
    weights = 1;
    idx = eq_idx;
    return
elseif ~isempty(eq_idx) && ~size(eq_idx,1)<2
    warning('%s: Grid is apparently colinear through origin.', upper(mfilename))
end

% Rotate to principal axes to enable 2D arrays
[x0, xs, status_2d] = rotate_to_principal_axes(x0, xs);

% Calculate dummy points to enable "partial" arrays
dummy_points = augment_bounding_box(x0);
if status_2d % Add dummy points in 2D case
    dummy_points = cat(1,dummy_points,[[0 0 1];[0 0 -1]]);
end
if ~isempty(dummy_points)
    dummy_indices = (1:size(dummy_points,1)) + size(x0,1);
    x0 = [x0; dummy_points];
end

%% ===== Computation =====================================================
% Delaunay triangulation of convex hull with and without xs
simplices_old = convhulln(x0);
simplices_new = convhulln([x0; xs]);

% Extract all neighbors of x0 sharing a triangle with xs, denoted as x0_s 
xs_idx = size([x0;xs], 1);
[row, ~] = find(simplices_new == xs_idx);
xs_tri = simplices_new(row, :);
idx = unique(xs_tri(xs_tri ~= xs_idx));
% x0_s = x0(idx, :);

% Extract all triangles from the simplices with at least one x0_s as a vertex
[simplices_new_s, simplices_old_s] = deal([]);
for n = 1:size(idx)
    [row, ~] = find(simplices_new == idx(n));
    [row_old, ~] = find(simplices_old == idx(n));
    simplices_new_s = cat(1, simplices_new_s, simplices_new(row, :));
    simplices_old_s = cat(1, simplices_old_s, simplices_old(row_old, :));
end
simplices_new_s = unique(simplices_new_s, 'rows');
simplices_old_s = unique(simplices_old_s, 'rows');

% Compute spherical voronoi_regions for each x0_s
[regions_new_s, vertices_new_s] = calc_voronoi_regions([x0;xs], center, ...
    simplices_new_s);
[regions_old_s, vertices_old_s] = calc_voronoi_regions([x0;xs], center, ...
    simplices_old_s);

% Prepare regions for calculation of the surface area
% Sorting the regions
[regions_new_s] = sort_voronoi_vertices_of_regions(simplices_new_s, ... 
    regions_new_s);
[regions_old_s] = sort_voronoi_vertices_of_regions(simplices_old_s, ...
    regions_old_s);

% Surface area calculation
[area_new, area_old] = deal(zeros(size(idx,1),1));
for i=1:size(idx,1)
    area_new(i) = calc_surface_area( ...
        vertices_new_s(regions_new_s{idx(i)},:), 1);
    area_old(i) = calc_surface_area( ...
        vertices_old_s(regions_old_s{idx(i)},:), 1);    
end

% Calculate weights
weights = (area_old - area_new) ./ sum(area_old - area_new);

if ~isempty(dummy_points)
    % Remove possible dummies from selected points
    dummy_mask = ismember(idx,dummy_indices);
    if any(dummy_mask)
        idx(dummy_mask) = [];
        if ~weights(dummy_mask)==0
            warning('%s: Requested point lies outside grid.', upper(mfilename))
        end
        weights(dummy_mask) = [];
    end
end

% Normalise weights
weights = weights./sum(weights);

[weights,order] = sort(weights,'descend');
idx = idx(order);
end

% =========================================================================

function [regions, vertices] = calc_voronoi_regions(x0, center, simplices)
%CALC_VORONOI_REGIONS calculates the voronoi vertices and regions of the given
% points x0. In case specific simplices are provided as input argument, only
% those will be included in the calculation.
%
%   Input parameters:
%       x0             - point cloud on the surface of a unit sphere in R^3 [nx3]
%       center         - center of the sphere in R^3 [1x3]
%       simplices      - pre-calculated simplices, each consisting of 
%                        3 indices of the points given in x0 [Nx3]
%
%   Output parameters:
%       regions        - cell array consisting of the indices of the vertices
%                        to belong to the region of each point of x0 [nx3]
%       vertices       - voronoi vertices in R^3 [Mx3]
%
%   Code based on scipy.spatial._spherical_voronoi.SphericalVoronoi

if nargin < 3
    simplices = convhulln(x0); % In case simplices are not provided as input arg
end

% Tetrahedrons from Delaunay triangulation with shape [4x3x2n-4]
% Add center of sphere to each of the simplices
tri = x0(simplices.', 1:end);
tetrahedrons=[];
for n=1:3:size(tri,1);
    tetrahedrons=cat(3,tetrahedrons,cat(1,tri(n:n+2,:),center));
end

% Calculate circumcenters of tetrahedrons
circumcenters=calc_circumcenters(tetrahedrons);

% Project circumcenters of the tetrahedrons to the surface of the unit
% sphere and thereby get the voronoi vertices with shape [2n-4x3]
vertices = bsxfun(@rdivide, circumcenters, vector_norm(circumcenters,2));

% Calculate regions from triangulation
% simplex_indices have shape [1x2n-4]
simplex_indices = 1:1:size(simplices,1);
% tri_indices have shape [1x6n-12]
tri_indices = sort( repmat(simplex_indices, 1,3) );
% point_indices have shape [1x6n-12]
point_indices = reshape(simplices.',1, []);
% array_associations have shape [6n-12x2]
array_associations = cat(1, point_indices, tri_indices).';
array_associations = sortrows(array_associations);
array_associations = cast(array_associations,'int32');

% Group by generator indices to produce unsorted regions in a cell array
regions = accumarray(array_associations(:,1),array_associations(:,2),[], ...
    @(x){x.'},{});
end

% =========================================================================

function [circumcenters] = calc_circumcenters(tetrahedrons)
%CALC_CIRCUMCENTERS calculates the circumcenters of the circumspheres of 
% tetrahedrons.
%
%   Input parameters:
%       tetrahedrons   - tetrahedrons defined by 4 points in R^3 [4x2xN]
%
%   Output parameters:
%       circumcenters  - circumcenters in R^3 [Nx3]
%
%   Code based on scipy.spatial._spherical_voronoi.calc_circumcenters
%   An implementation based on http://mathworld.wolfram.com/Circumsphere.html

num = size(tetrahedrons,3);
a = cat(2, tetrahedrons, ones(4,1,num));

sums = sum(tetrahedrons.^2, 2);
d = cat(2, sums, a);

[dx1, dy1, dz1] = deal(d);
dx1(:, 2, :) = [];
dy1(:, 3, :) = [];
dz1(:, 4, :) = [];

% Calculating det in 3D array; inefficient
[dx, dy, dz, ad] = deal(zeros(1,num));
for i = 1:num
    dx(i) = det(dx1(:,:,i));
    dy(i) = -det(dy1(:,:,i));
    dz(i) = det(dz1(:,:,i));
    ad(i) = det(a(:,:,i));
end

nominator = cat(1, dx, dy, dz);
denominator = 2.*ad;

circumcenters = bsxfun(@rdivide, nominator, denominator).';
end

% =========================================================================

function [ regions_sorted ] = sort_voronoi_vertices_of_regions( simplices, regions )
%SORT_VORONOI_VERTICES_OF_REGIONS sorts the indices of the voronoi vertices for
% each region such that the resulting points are in a clockwise or
% counterclockwise order around the generator point.
%
%   Input parameters:
%       simplices      - array of shape (Nx3) consisting of N points on 
%                        the surface of a unit sphere in R^3 
%       regions        - center of the sphere in R^3 [1x3]
%
%   Output parameters:
%       regions        - cell array consisting of the indices of the vertices to
%                        belong to the region of each point of x0 [Nx3]
%
%   Code based on scipy.spatial._spherical_voronoi.SphericalVoronoi

ARRAY_FILLER = -2;
sorted_vertices = zeros(1, max(cellfun('length', regions)), 'int32');
regions_sorted = cell(size(regions,1),1);

for n=1:size(regions,1)
    remaining_count = 1;
    remaining = regions{n};
    remaining_size = size(remaining,2);
    sorted_vertices(:) = ARRAY_FILLER;
    
    % In case the region (cell) is empty
    if remaining_size<1
        regions_sorted{n} = 0;
    else
        current_simplex = remaining(1);
        for i=1:3
            k = simplices(current_simplex,i);
            if k ~= n
                current_vertex = k;
                break
            end
        end
        sorted_vertices(remaining_count) = current_simplex;
        remaining_count = remaining_count+1;
        
        % remaining_filter()
        for l=1:size(remaining,2)
            if remaining(l) == current_simplex
                remaining(l) = ARRAY_FILLER;
            end
        end
        
        while remaining_size >= remaining_count
            cs_identified = 0;
            for i=1:remaining_size
                if remaining(i) == ARRAY_FILLER
                    continue
                end
                s = remaining(i);
                for j=1:3
                    if current_vertex == simplices(s, j)
                        current_simplex = remaining(i);
                        cs_identified = cs_identified +1;
                        break
                    end
                end
                if cs_identified > 0
                    break
                end
            end
            
            for i=1:3
                s = simplices(current_simplex, i);
                if s~=n && s~=current_vertex
                    current_vertex = s;
                    break
                end
            end
            sorted_vertices(remaining_count) = current_simplex;
            remaining_count = remaining_count +1;
            % remaining_filter()
            for l=1:size(remaining,2)
                if remaining(l) == current_simplex
                    remaining(l) = ARRAY_FILLER;
                end
            end
        end
        regions_sorted{n} = sorted_vertices(sorted_vertices > ARRAY_FILLER);
    end
end
end

% =========================================================================

function [ surface_area ] = calc_surface_area( vertices, radius )
%CALC_SURFACE_AREA calculates the surface area of a polygon on the surface of a
% sphere. The input vertices need to be sorted around their generator point.
%
%   Input parameters:
%       vertices       - sorted vertices [Nx3]
%       radius         - radius of the sphere
%
%   Output parameters:
%       surface_area   - calculated surface area of the given polygon
%
%   Implementation based on: http://mathworld.wolfram.com/LHuiliersTheorem.html
%   Code based on: voronoi_utility.py from https://doi.org/10.5281/zenodo.13688

n = size(vertices, 1);
root_point = vertices(1,:);
surface_area = 0;

b_point = vertices(2,:);
root_b_dist = calc_haversine_dist(root_point, b_point);
for i=2:(n-1)
    a_point = b_point;
    b_point = vertices(i+1,:);
    root_a_dist = root_b_dist;
    root_b_dist = calc_haversine_dist(root_point, b_point);
    a_b_dist = calc_haversine_dist(a_point, b_point);
    s = (root_a_dist + root_b_dist + a_b_dist) ./ 2;
    surface_area = surface_area + 4.* atan( sqrt( ...
        tan(0.5 .* s) .* ...
        tan(0.5 .* (s-root_a_dist)) .* ...
        tan(0.5 .* (s-root_b_dist)) .* ...
        tan(0.5 .* (s-a_b_dist)) ));
end
% Avoid complex numbers, which occur possibly due to rounding errors
surface_area = real(surface_area * radius^2);
end

% =========================================================================

function [ sph_dist ] = calc_haversine_dist( x1, x2 )
%CALC_HAVERSINE_DISTANCE calculates the distance between two points on the 
% surface of a unit sphere based on the haversine formula provided here: 
% https://en.wikipedia.org/wiki/Haversine_formula
%
%   Code based on: voronoi_utility.py from https://doi.org/10.5281/zenodo.13688

[az1, el1, ~] = cart2sph(x1(:,1), x1(:,2), x1(:,3));
[az2, el2, ~] = cart2sph(x2(:,1), x2(:,2), x2(:,3));

el1 = -(el1-pi/2);
el2 = -(el2-pi/2);

% long/lat is not the same as spherical coordinates - phi differs by pi/4
sph_dist = 2 * asin( sqrt(( (1 - cos(el2 - el1))./2) ...
                          + sin(el1) .* sin(el2) .* ((1-cos(az2 - az1))./2)) );
end

% =========================================================================

function [x0, xs, status_2d] = rotate_to_principal_axes(x0, xs, gamma)
%ROTATE_TO_PRINCIPAL_AXES rotates x0 and xs to x0's principal axes.
% If the ratio of second-to-smallest singular values is < gamma,
% the third dimension of x0 is set to zero.
%
%   Input parameters:
%       x0          - point cloud in R^3 [nx3]
%       xs          - point in R^3 [1x3]
%       gamma       - scalar in 0 < gamma << 1
%
%   Output parameters:
%       x0          - point cloud in R^3
%       xs          - point in R^3
%       status_2d   - true or false
if nargin < 3
    gamma = 0.1; % inverse of aspect ratio of principal axes
end
status_2d = 0;

[~,S,V] = svd(x0);
x0 = x0*V;
xs = xs*V;
S = diag(S);
if S(end)/S(end-1) < gamma
     x0(:,3) = 0;
     %xs(3) = 0;
    status_2d = 1;
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
