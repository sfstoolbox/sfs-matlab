function x0 = secondary_source_positions(conf)
%SECONDARY_SOURCE_POSITIONS generates the positions, directions, and weights of
%   the secondary sources
%
%   Usage: x0 = secondary_source_positions(conf)
%
%   Input options:
%       conf   - configuration struct (see SFS_config)
%
%   Output options:
%       x0     - secondary source positions, directions and weights
%                [x0 y0 z0 nx0 ny0 nz0 w] / m
%
%   SECONDARY_SOURCES_POSITIONS(conf) generates the positions and directions
%   x0 of secondary sources for a given geometry
%   (conf.secondary_sources.geometry) and array size
%   (conf.secondary_sources.size). Alternatively, if
%   conf.secondary_sources.geomrtry is set to 'custom' the field
%   conf.secondary_sources.x0 is used to supply x0. It can be a [n 7] matrix
%   consiting of n sources or it can be a SOFA file/struct from the source
%   positions are extracted.
%
%   The direction of the sources is given as their unit vectors pointing in the
%   given direction. For a linear array the secondary sources are pointing
%   towards the negative y-direction. If you create a linear array with default
%   position conf.secondary_sources.center = [0 0 0], your listening area is in
%   the area y<0, which means the y value of conf.xref should also be <0!
%
% See also: secondary_source_selection, secondary_source_tapering

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

% NOTE: If you wanted to add a new type of loudspeaker array, do it in a way,
% that the loudspeakers are ordered in a way, that one can go around for closed
% arrays. Otherwise the tapering window function will not work properly.


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargstruct(conf);


%% ===== Configuration ===================================================
% Given secondary sources are used in the 'custom' section
%conf.secondary_sources.x0;
% Array type
geometry = conf.secondary_sources.geometry;
if ~strcmp('custom',geometry)
    % Center of the array
    X0 = conf.secondary_sources.center;
    % Number of secondary sources
    nls = conf.secondary_sources.number;
    x0 = zeros(nls,7);
    % Diameter/length of array
    L = conf.secondary_sources.size;
end


%% ===== Main ============================================================
if strcmp('line',geometry) || strcmp('linear',geometry)
    % === Linear array ===
    %
    %                     y-axis
    %                       ^
    %                       |
    %                       |  secondary sources
    %                       |        |
    %                       |        v
    %  ------x--x--x--x--x--x--x--x--x--x--x-------> x-axis
    %        |  |  |  |  |  |  |  |  |  |  | <- secondary source direction
    %                       |
    %                       |
    %
    %% Positions of the secondary sources
    x0(:,1) = X0(1) + linspace(-L/2,L/2,nls)';
    x0(:,2) = X0(2) * ones(nls,1);
    x0(:,3) = X0(3) * ones(nls,1);
    % Direction of the secondary sources pointing to the -y direction
    x0(:,4:6) = direction_vector(x0(:,1:3),x0(:,1:3)+repmat([0 -1 0],nls,1));
    % Weight each secondary source by the inter-loudspeaker distance
    x0(:,7) = L./(nls-1);
elseif strcmp('circle',geometry) || strcmp('circular',geometry)
    % === Circular array ===
    %
    %                  y-axis
    %                    ^
    %                    |
    %                    x
    %               x    |     x
    %               \    |     /
    %          x_        |         _x
    %            -       |        -
    %       x-_          |          _-x
    %                    |
    %  ----x---------------------------x------> x-axis
    %         _          |          _
    %       x-           |           -x
    %          _-        |        -_
    %         x          |          x
    %              /     |     \
    %              x     |     x
    %                    x
    %                    |
    %
    % 'circle' is special case of 'rounded-box' with fully rounded corners
    t = (0:nls-1)/nls;
    [x0(:,1:3), x0(:,4:6), x0(:,7)] = rounded_box(t,1.0);  % 1.0 for circle
    % Scale unit circle and shift center to X0
    x0(:,1:3) = bsxfun(@plus, x0(:,1:3).*L/2, X0);
    % Scale weights
    x0(:,7) = x0(:,7).*L/2;
elseif strcmp('box',geometry)
    % === Boxed loudspeaker array ===
    %
    %                  y-axis
    %                    ^
    %                    |
    %        x   x   x   x   x   x   x
    %        |   |   |   |   |   |   |
    %     x--            |            --x
    %                    |
    %     x--            |            --x
    %                    |
    %     x--            |            --x
    %                    |
    %  ---x-----------------------------x-----> x-axis
    %                    |
    %     x--            |            --x
    %                    |
    %     x--            |            --x
    %                    |
    %     x--            |            --x
    %        |   |   |   |   |   |   |
    %        x   x   x   x   x   x   x
    %                    |
    %
    % 'box' is special case of 'rounded-box' where there is no rounding
    % and the sources in the corners are skipped
    %
    % Number of secondary sources per linear array
    % ensures that nls/4 is always an integer.
    if rem(nls,4)~=0
        error(['%s: conf.secondary_sources.number has to be a multiple of' ...
            ' 4.'],upper(mfilename));
    else
        nbox = nls/4;
    end
    % Distance between secondary sources
    dx0 = L/(nbox-1);
    % Length of one edge of the rectangular bounding box
    Lbound = L + 2*dx0;
    % Index t for the positions on the boundary
    t = linspace(-L/2,L/2,nbox)./Lbound;  % this skips the corners
    t = [t, t+1, t+2, t+3]*0.25;  % repeat and shift to get all 4 edges
    % 'box' is special case of 'rounded-box' where there is no rounding
    [x0(:,1:3), x0(:,4:6), x0(:,7)] = rounded_box(t,0.0);  % 0.0 for square
    % Scale "unit" box and shift center to X0
    x0(:,1:3) = bsxfun(@plus, x0(:,1:3).*Lbound/2, X0);
    % Scale integration weights
    x0(:,7) = x0(:,7).*Lbound/2;
    % Correct weights of loudspeakers near corners
    corners = [1,nbox,nbox+1,2*nbox,2*nbox+1,3*nbox,3*nbox+1,4*nbox];
    x0(corners,7) = (1 + sqrt(2)) *dx0/2;  % instead of 3/2 * dx0
elseif strcmp('rounded-box', geometry)
    % Ratio for rounding the edges
    ratio = 2*conf.secondary_sources.corner_radius./L;
    t = (0:nls-1)/nls;
    [x0(:,1:3), x0(:,4:6), x0(:,7)] = rounded_box(t, ratio);
    % Scale "unit" rounded-box and shift center to X0
    x0(:,1:3) = bsxfun(@plus, x0(:,1:3).*L/2, X0);
    % Scale integration weights
    x0(:,7) = x0(:,7).*L/2;
elseif strcmp('spherical',geometry) || strcmp('sphere',geometry)
    % Get spherical grid + weights
    [points,weights] = get_spherical_grid(nls,conf);
    % Secondary source positions
    x0(:,1:3) = L/2 * points + repmat(X0,nls,1);
    % Secondary source directions
    x0(:,4:6) = direction_vector(x0(:,1:3),repmat(X0,nls,1));
    % Secondary source weights + distance scaling
    x0(:,7) = weights .* L^2/4;
elseif strcmp('custom',geometry)
    % Custom geometry defined by conf.secondary_sources.x0.
    % This could be in the form of a n x 7 matrix, where n is the number of
    % secondary sources or as a SOFA file/struct.
    if ischar(conf.secondary_sources.x0) || isstruct(conf.secondary_sources.x0)
        x0 = sofa_get_secondary_sources(conf.secondary_sources.x0);
    else
        x0 = conf.secondary_sources.x0;
    end
    isargsecondarysource(x0);
else
    error('%s: %s is not a valid array geometry.',upper(mfilename),geometry);
end
