function x0 = secondary_source_positions(conf)
%SECONDARY_SOURCE_POSITIONS Generates the positions and directions of the
%   secondary sources
%
%   Usage: x0 = secondary_source_positions([conf])
%
%   Input options:
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       x0          - secondary source positions, directions and weights / m
%
%   SECONDARY_SOURCES_POSITIONS(conf) generates the positions and directions
%   x0 of secondary sources for a given geometry
%   (conf.secondary_sources.geometry) and array size
%   (conf.secondary_sources.size). Alternatively, if conf.secondary_sources.x0
%   is set, it returns the positions and directions specified there.
%   The direction of the sources is given as their unit vectors pointing in the
%   given direction. For a linear array the secondary sources are pointing
%   towards the negative y-direction. If you create a linear array with default
%   position conf.secondary_sources.center = [0 0 0], your listening area is in
%   the area y<0, which means the y value of conf.xref should also be <0!
%
%   Default geometry for a linear array:
%
%                                y-axis
%                                   ^
%                                   |
%                                   |  secondary sources
%                                   |        |
%                                   |        v
%       -------------x--x--x--x--x--x--x--x--x--x--x------------> x-axis
%                    |  |  |  |  |  |  |  |  |  |  | <- secondary source direction
%                                   |              
%                                   |
%                                   |
%
%   Default geometry for a circular/spherical array:
%
%                                y-axis
%                                   ^
%                                   |
%                                   x
%                              x    |     x
%                              \    |     /
%                         x_        |         _x
%                           -       |        -
%                      x-_          |          _-x
%                                   |         
%       --------------x---------------------------x------------------> x-axis
%                        _          |          _
%                      x-           |           -x
%                         _-        |        -_
%                        x          |          x
%                             /     |     \
%                             x     |     x
%                                   x
%                                   |
%
%
%   Default geometry for a box-shape array:
%
%                                y-axis
%                                   ^
%                                   |
%                       x   x   x   x   x   x   x  
%                       |   |   |   |   |   |   |            
%                    x--            |            --x
%                                   |         
%                    x--            |            --x
%                                   |
%                    x--            |            --x
%                                   |
%       -------------x-----------------------------x-----------------> x-axis
%                                   |
%                    x--            |            --x
%                                   |
%                    x--            |            --x
%                                   |
%                    x--            |            --x
%                       |   |   |   |   |   |   |
%                       x   x   x   x   x   x   x
%                                   |
%
% see also: secondary_source_selection, secondary_source_tapering 

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************

% NOTE: If you wanted to add a new type of loudspeaker array, do it in a way,
% that the loudspeakers are ordered in a way, that one can go around for closed
% arrays. Otherwise the tapering window function will not work properly.


%% ===== Checking of input  parameters ===================================
nargmin = 0;
nargmax = 1;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
% Given secondary sources
x0 = conf.secondary_sources.x0;
% Check if we have already predefined secondary sources
if ~isempty(x0)
    isargsecondarysource(x0);
    % If we have predefined secondary sources return at this point
    return
end
% Array type
geometry = conf.secondary_sources.geometry;
% Center of the array
X0 = conf.secondary_sources.center;
% Number of secondary sources
nls = conf.secondary_sources.number;
% Diameter/length of array
L = conf.secondary_sources.size;


%% ===== Main ============================================================
x0 = zeros(nls,7);
if strcmp('line',geometry) || strcmp('linear',geometry)
    % === Linear array ===
    % Positions of the secondary sources
    x0(:,1) = X0(1) + linspace(-L/2,L/2,nls)';
    x0(:,2) = X0(2) * ones(nls,1);
    x0(:,3) = X0(3) * ones(nls,1);
    % Direction of the secondary sources pointing to the -y direction
    x0(:,4:6) = direction_vector(x0(:,1:3),x0(:,1:3)+repmat([0 -1 0],nls,1));
    % equal weights for all sources
    x0(:,7) = ones(nls,1);
elseif strcmp('circle',geometry) || strcmp('circular',geometry)
    % === Circular array ===
    % Azimuth angles
    phi = linspace(0,(2-2/nls)*pi,nls)'; % 0..2pi
    %phi = linspace(pi/2,(5/2-2/nls)*pi,nls)'; % pi/2..5/2pi, Room Pinta
    % Elevation angles
    theta = zeros(nls,1);
    % Positions of the secondary sources
    [cx,cy,cz] = sph2cart(phi,theta,L/2);
    x0(:,1:3) = [cx,cy,cz] + repmat(X0,nls,1);
    % Direction of the secondary sources
    x0(:,4:6) = direction_vector(x0(:,1:3),repmat(X0,nls,1).*ones(nls,3));  
    % equal weights for all sources
    x0(:,7) = ones(nls,1);
elseif strcmp('box',geometry)
    % === Boxed loudspeaker array ===
    % Number of secondary sources per linear array
    % ensures that nls/4 is always an integer.
    nbox = round(nls/4);
    % distance between secondary sources
    dx0 = L/(nbox-1);
    % Position and direction of the loudspeakers
    % top
    x0(1:nbox,1) = X0(1) + linspace(-L/2,L/2,nbox)';
    x0(1:nbox,2) = X0(2) + ones(nbox,1) * L/2 + dx0;
    x0(1:nbox,3) = X0(3) + zeros(nbox,1);
    x0(1:nbox,4:6) = direction_vector(x0(1:nbox,1:3), ...
        x0(1:nbox,1:3)+repmat([0 -1 0],nbox,1));
    % right
    x0(nbox+1:2*nbox,1) = X0(1) + ones(nbox,1) * L/2 + dx0;
    x0(nbox+1:2*nbox,2) = X0(2) + linspace(L/2,-L/2,nbox)';
    x0(nbox+1:2*nbox,3) = X0(3) + zeros(nbox,1);
    x0(nbox+1:2*nbox,4:6) = direction_vector(x0(nbox+1:2*nbox,1:3), ...
        x0(nbox+1:2*nbox,1:3)+repmat([-1 0 0],nbox,1));
    % bottom
    x0(2*nbox+1:3*nbox,1) = X0(1) + linspace(L/2,-L/2,nbox)';
    x0(2*nbox+1:3*nbox,2) = X0(2) - ones(nbox,1) * L/2 - dx0;
    x0(2*nbox+1:3*nbox,3) = X0(3) + zeros(nbox,1);
    x0(2*nbox+1:3*nbox,4:6) = direction_vector(x0(2*nbox+1:3*nbox,1:3), ...
        x0(2*nbox+1:3*nbox,1:3)+repmat([0 1 0],nbox,1));
    % left
    x0(3*nbox+1:nls,1) = X0(1) - ones(nbox,1) * L/2 - dx0;
    x0(3*nbox+1:nls,2) = X0(2) + linspace(-L/2,L/2,nbox)';
    x0(3*nbox+1:nls,3) = X0(3) + zeros(nbox,1);
    x0(3*nbox+1:nls,4:6) = direction_vector(x0(3*nbox+1:nls,1:3), ...
        x0(3*nbox+1:nls,1:3)+repmat([1 0 0],nbox,1));
    % equal weights for all sources
    x0(:,7) = ones(nls,1);
elseif strcmp('spherical',geometry) || strcmp('sphere',geometry)
    % get spherical grid + weights
    [points,weights] = get_spherical_grid(nls,conf);
    % secondary source positions
    x0(:,1:3) = L/2 * points + repmat(X0,nls,1);
    % secondary source directions
    x0(:,4:6) = direction_vector(x0(:,1:3),repmat(X0,nls,1));
    % secondary source weights
    x0(:,7) = weights;
    % add integration weights (because we integrate over a sphere) to the grid
    % weights
    [~,theta] = cart2sph(x0(:,1),x0(:,2),x0(:,3)); % get elevation
    x0(:,7) = x0(:,7) .* cos(theta);
    
else
    error('%s: %s is not a valid array geometry.',upper(mfilename),geometry);
end
