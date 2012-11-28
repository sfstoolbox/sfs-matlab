function x0 = secondary_source_positions(L,conf)
%SECONDARY_SOURCE_POSITIONS Generates the positions and directions of the
%   secondary sources
%
%   Usage: x0 = secondary_source_positions(L,[conf])
%
%   Input options:
%       L           - the size (m) of the array (length for a linear array,
%                     diameter for a circle or a box)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output options:
%       x0          - secondary source positions and directions (m)
%
%   SECONDARY_SOURCES_POSITIONS(L,conf) generates the positions and directions
%   x0 of secondary sources for a given geometry (conf.array) and array size
%   (L). Alternatively, if conf.x0 is set, it returns the positions and
%   directions specified there.
%   The direction of the sources is given as their unit vectors pointing in the
%   given driection.
%
%   Geometry (for a linear array):
%
%                                y-axis
%                                   ^
%                                   |
%                                   |
%                                   |
%          v--v--v--v--v--v--v--v--v|-v--v
%          |              X0        |
%    (Loudspeaker)  (Array center)  |
%                                   |
%       --------------------------------------------------------> x-axis
%
% see also: secondary_source_selection, secondary_source_number, tapering_window

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************

% NOTE: If you wanted to add a new type of loudspeaker array, do it in a way,
% that the loudspeakers are ordered in a way, that one can go around for closed
% arrays. Otherwise the tapering window function will not work properly.


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(L);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
% Array type
array = conf.array;
% Center of the array
X0 = position_vector(conf.X0);
% Distance between secondary sources
dx0 = conf.dx0;
% Given secondary sources
x0 = conf.x0;


%% ===== Main ============================================================
% Check if we have already predefined secondary sources
if ~isempty(x0)
    isargsecondarysource(x0);
    % If we have predefined secondary sources return at this point
    return
end

% Get the number of secondary sources
[nls, L] = secondary_source_number(L,conf);

x0 = zeros(nls,6);
if strcmp('linear',array)
    % === Linear array ===
    % Positions of the secondary sources
    x0(:,1) = X0(1) + linspace(-L/2,L/2,nls)';
    x0(:,2) = X0(2) * ones(nls,1);
    x0(:,3) = X0(3) * ones(nls,1);
    % Direction of the secondary sources
    x0(:,4:6) = direction_vector(x0(:,1:3),x0(:,1:3)+repmat([0 1 0],nls,1));
elseif strcmp('circle',array) || strcmp('circular',array)
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
elseif strcmp('box',array)
    % === Boxed loudspeaker array ===
    % Number of secondary sources per linear array
    % Note, that the call to the secondary_source_number function above
    % ensures that nls/4 is always an integer.
    nbox = nls/4;
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
elseif strcmp('U',array)
    to_be_implemented(mfilename);
else
    error('%s: %s is not a valid array type.',upper(mfilename),array);
end
