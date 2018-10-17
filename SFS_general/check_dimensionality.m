function [x0, xs, dim, eq_idx] = check_dimensionality(x0, xs, tol, gamma)
%CHECK_DIMENSIONALITY checks dimensionality and rotates the grid to its
%principal axes.
%
%   Usage: [x0,xs,dim,eq_idx] = check_dimensionality(x0,xs,tol)
%
%   Input parameters:
%       x0          - point cloud on a sphere around the origin / m [nx3]
%       xs          - desired direction as point in R^3 / m [1x3]
%       tol         - tolerance for equality check
%       gamma       - scalar in 0 < gamma << 1
%
%   Output parameters:
%       x0          - original point cloud rotated to principal axes / m [nx3]
%       xs          - original query point rotated to principal axes / m [1x3]
%       dim         - dimensionality, either 1, 2, 2.5 or 3
%       eq_idx      - indice of point from point cloud equal to query point 
%                     within tolerance specified in tol or -1 otherwise
%
%   1D:    xs is colinear with or equal to one x0
%   2D:    all x0 and xs are in one plane
%   2.5D:  all x0 are in one plane, but xs is not
%   2.55D: 2.5D case with coplanarity of grid and the center of the sphere
%   3D:    otherwise
%
%   See also: findvoronoi, test_interpolation_point_selection,
%   rotate_to_principal_axes
%
%

if nargin < 3
    tol = 1e-6; % in case no tolerance is provided as input arg
    gamma = 0.1; % inverse of aspect ratio of principal axes
end

% In case no 1D or 2D case is detected
dim = 3;

% Check for 1D case (equality\collinearity within tolerance) and return idx
eq = vector_norm(bsxfun(@minus,x0,xs),2);
eq_idx = find(eq<=tol);
if ~isempty(eq_idx)
    dim = 1;
    warning('SFS:check_dimensionality',...
        '%s: Query point is apparently colinear with or equal to one grid point.', ...
        upper(mfilename))
    return
else
    eq_idx = -1;
end

% Check for 2D or 2.5D case and rotate to principal axes
x0mean = mean(x0,1);
[~,S,V] = svd(bsxfun(@minus,x0,x0mean));
x0mean = x0mean*V;
x0 = x0*V;
xs = xs*V;
S = diag(S);
if S(end)/S(end-1) < gamma
    dim = 2.5;
    warning('SFS:check_dimensionality',...
        '%s: Grid is apparently two-dimensional. ',upper(mfilename));
    % Check for 2D case
    if abs(xs(3) - x0mean(3)) < tol
        dim = 2;
        x0(:,3) = x0(:,3) - x0mean(3);
        xs(3) = xs(3) - x0mean(3);
    % Check for coplanarity of grid and center of the sphere, denoted as 2.55D
    elseif abs(x0mean(3)) < tol
        dim = 2.55;
    end
end

end