function [idx,weights] = findvoronoi(x0,xs)
%FINDVORONOI finds the corresponding voronoi regions to the points x0
%surrounding the desired point xs. The weights are derived from the Voronoi
%region surface area differences with and without xs.
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
%   See also: findnearestneighbour, findconvexcone,
%             test_interpolation_point_selection

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2018 SFS Toolbox Developers                             *
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

% Check dimensionality and rotate to principal axes
[x0,xs,dim,eq_idx] = check_dimensionality(x0,xs);

% In 1D case (xs is colinear with or equals one x0) no interpolation is needed
if dim==1
    idx = eq_idx;
    weights = ones(size(idx,1),1) / size(idx,1);
    return
end

% In 2D case linear interpolation using neighboring x0 is sufficient
if dim==2
   [idx, weights] = findnearestneighbour(x0, xs, 2);
   return
end

% In 2.5D case order x0 with respect to azimuth angle
if dim == 2.5
    [~,idx_sorted,az] = sort_azimuth(x0);
end

% If the grid and the center of the unit sphere are coplanar in an otherwise
% 2.5D case, dummy points are used as a workaround, preventing degenerate
% tetrahedra as part of the calculation of the new voronoi vertices.
% This case is denoted as 2.55D and is currently handled like a 3D case with
% dummy points added to the grid.
dummy_points = [];
if dim == 2.55
    dummy_points = [[0 0 -1];[0 0 1]]; % Add north and south pole
    dummy_indices = (1:size(dummy_points,1)) + size(x0,1);
    x0 = [x0;dummy_points];
end

%% ===== Computation =====================================================

% Delaunay triangulation of convex hull with xs (new)
simplices_new = convhulln([x0; xs]);

% Extract all neighbors of x0 sharing a triangle with xs, denoted as x0_s
xs_idx = size([x0;xs], 1);
[row, ~] = find(simplices_new == xs_idx);
xs_tri = simplices_new(row, :);
idx = unique(xs_tri(xs_tri ~= xs_idx));

% Extract all triangles from the simplices with at least one x0_s as a vertex
simplices_new_s = [];
for n = 1:size(idx)
    [row, ~] = find(simplices_new == idx(n));
    simplices_new_s = cat(1,simplices_new_s,simplices_new(row,:));
end
simplices_new_s = unique(simplices_new_s,'rows');

% Compute new spherical voronoi_regions for each x0_s
[regions_new_s, vertices_new_s] = calc_voronoi_regions([x0;xs],center, ...
    simplices_new_s);

% Prepare new regions for calculation of the surface area
% Sorting the regions
[regions_new_s] = sort_voronoi_vertices_of_regions(simplices_new_s, ... 
    regions_new_s);

% Special 2.5D case handling: skip computation of old voronoi regions
if dim ~= 2.5
    % Delaunay triangulation of convex hull without xs (old)
    simplices_old = convhulln(x0);
    
    % Extract all triangles from the simplices with at least one x0_s as vertex
    simplices_old_s = [];
    for n = 1:size(idx)
        [row_old, ~] = find(simplices_old == idx(n));
        simplices_old_s = cat(1,simplices_old_s,simplices_old(row_old,:));
    end
    simplices_old_s = unique(simplices_old_s,'rows');
    
    % Compute old spherical voronoi_regions for each x0_s
    [regions_old_s, vertices_old_s] = calc_voronoi_regions([x0;xs],center, ...
        simplices_old_s);
    
    % Prepare old regions for calculation of the surface area
    % Sorting the regions
    [regions_old_s] = sort_voronoi_vertices_of_regions(simplices_old_s, ...
        regions_old_s);
    
end

% Surface area calculation
[area_new, area_old] = deal(zeros(size(idx,1),1));
for ii=1:size(idx,1)
    area_new(ii) = calc_surface_area( ...
        vertices_new_s(regions_new_s{idx(ii)},:),1);
    % Special 2.5D case handling: alternative calculation of old area
    if dim ~= 2.5
        area_old(ii) = calc_surface_area( ...
            vertices_old_s(regions_old_s{idx(ii)},:),1);  
    else
        % Old area calculation based on area of spherical lune
        az_diff = zeros(size(idx,1),1);
        m = [size(idx,1) 1:size(idx,1) 1];
        for n = 1:size(idx,1)
            az_diff(n) = (az(m(n+2)) - az(m(n))) / 2;
            if az_diff(n)<0
                az_diff(n) = (az(m(n+2)) - az(m(n)) + 360) / 2;
            end
        end
        area_old(:) = pi/90 * az_diff(:);
        area_old = area_old(circshift(idx_sorted,2));
    end
end

% Calculate weights
weights = (area_old - area_new) ./ sum(area_old - area_new);

% Remove possible dummies from selected points
if ~isempty(dummy_points)
    dummy_mask = ismember(idx,dummy_indices);
    if any(dummy_mask)
        idx(dummy_mask) = [];
        if ~weights(dummy_mask)==0
            warning('%s: Requested point lies outside grid.', upper(mfilename))
        end
        weights(dummy_mask) = [];
    end
end

[weights,order] = sort(weights,'descend');
idx = idx(order);
end


%% ===== Functions =======================================================

function [regions,vertices] = calc_voronoi_regions(x0,center,simplices)
%CALC_VORONOI_REGIONS calculates the Voronoi vertices and regions of the given
%points x0. In case specific simplices are provided as input argument, only
%those will be included in the calculation.
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
%       vertices       - Voronoi vertices in R^3 [Mx3]
%
%   Code based on scipy.spatial._spherical_voronoi.SphericalVoronoi

if nargin < 3
    simplices = convhulln(x0); % in case simplices are not provided as input arg
end

% Tetrahedrons from Delaunay triangulation with shape [4x3x2n-4].
% Add center of sphere to each of the simplices
tri = x0(simplices.', 1:end);
tetrahedrons = [];
n = 1:3:size(tri,1);
for m=n
    tetrahedrons = cat(3,tetrahedrons,cat(1,tri(m:m+2,:),center));
end

% Calculate surface normal of each triangle via cross product of triangle edges
N = bsxfun(@cross,tri(n+1,:)-tri(n,:),tri(n+2,:)-tri(n,:));
% Determine direction of projection into correct hemisphere
project_dir = sign(vector_product(direction_vector(tri(n,:),[0 0 0]),N(:,:),2));

% Calculate circumcenters of tetrahedrons
circumcenters = calc_circumcenters(tetrahedrons);

% Project circumcenters of the tetrahedrons to the surface of the unit
% sphere and thereby get the voronoi vertices with shape [2n-4x3]
% Consider the surface normal direction of each triangle for projection into the
% correct hemisphere
circumcenters = bsxfun(@times, circumcenters, project_dir);
vertices = bsxfun(@rdivide,circumcenters,vector_norm(circumcenters,2));

% Calculate regions from triangulation.
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
%tetrahedrons.
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
for ii=1:num
    dx(ii) = det(dx1(:,:,ii));
    dy(ii) = -det(dy1(:,:,ii));
    dz(ii) = det(dz1(:,:,ii));
    ad(ii) = det(a(:,:,ii));
end

nominator = cat(1, dx, dy, dz);
denominator = 2.*ad;

circumcenters = bsxfun(@rdivide, nominator, denominator).';
end

% =========================================================================

function [regions_sorted] = sort_voronoi_vertices_of_regions(simplices,regions)
%SORT_VORONOI_VERTICES_OF_REGIONS sorts the indices of the Voronoi vertices for
%each region such that the resulting points are in a clockwise or
%counterclockwise order around the generator point.
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
        for ii=1:3
            k = simplices(current_simplex,ii);
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
            for ii=1:remaining_size
                if remaining(ii) == ARRAY_FILLER
                    continue
                end
                s = remaining(ii);
                for jj=1:3
                    if current_vertex == simplices(s, jj)
                        current_simplex = remaining(ii);
                        cs_identified = cs_identified +1;
                        break
                    end
                end
                if cs_identified > 0
                    break
                end
            end
            
            for ii=1:3
                s = simplices(current_simplex, ii);
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

function [surface_area] = calc_surface_area(vertices,radius)
%CALC_SURFACE_AREA calculates the surface area of a polygon on the surface of a
%sphere. The input vertices need to be sorted around their generator point.
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
for ii=2:(n-1)
    a_point = b_point;
    b_point = vertices(ii+1,:);
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

function [sph_dist] = calc_haversine_dist(x1,x2)
%CALC_HAVERSINE_DISTANCE calculates the distance between two points on the 
% surface of a unit sphere based on the Haversine formula provided here: 
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

function [x0,idx,az] = sort_azimuth( x0 )
%SORT_AZIMUTH sorts x0 based on azimuth angle

% Compute azimuth
az = atan2d(x0(:,2), x0(:,1));
% Sort azimuth
[az, idx] = sort(az);
% Reorder
x0 = [x0(idx,1) x0(idx,2) x0(idx,3)];
end