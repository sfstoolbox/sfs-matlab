function [points,weights] = get_spherical_grid(number,conf)
%GET_SPHERICAL_GRID retruns grid points and weights
%
%   Usage: [points,weights] = get_spherical_grid(number,[conf])
%
%   Input parameters:
%       number  - number of grid points
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       points  - grid points
%       weights - integration weights for the grid points
%
%   GET_SPHERICAL_GRID(number) returns the points and weights for a grid on a
%   sphere. The type of grid is specified by conf.secondary_sources.grid.
%   For available grids, have a look at http://github.com/sfstoolbox/data.
%   It expects the grid files at SFS_basepath/data/spherical_grids. If the
%   desired file is not available on the hard disk, the function tries to
%   download it directly from github.
%   For conf.secondary_sources.grid='gauss' the grid positions are calculated
%   after Ahrens (2012), p. 121
%
%   References:
%       J. Ahrens (2012) - "Analytic Methods of Sound Field Synthesis", Springer.
%
%   see also: secondary_source_positions,
%       weights_for_points_on_a_sphere_rectangle

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


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(number);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
spherical_grid = conf.secondary_sources.grid;


%% ===== Main ============================================================
filename = sprintf('%06.0fpoints.mat',number);
basepath = get_sfs_path();

if strcmp('equally_spaced_points',spherical_grid)
    % check if we have a squared number of points, because only for those values
    % are equally spaced points grids available.
    if mod(number,sqrt(number))~=0
        error('%s: number has to be a squared number.',upper(mfilename));
    end
    file = [basepath '/data/spherical_grids/equally_spaced_points/' filename];
    url = ['https://dev.qu.tu-berlin.de/projects/data/repository/revisions/' ...
        'master/raw/spherical_grids/equally_spaced_points/' filename];
    % download file if not present
    if ~exist(file,'file')
        download_file(url,file);
    end
    tmp = load(file,'-ascii');
    points = tmp(:,1:3);
    weights = tmp(:,4);
elseif strcmp('fabian',spherical_grid)
    % here we have only one number of secondary sources available
    if number~=11345
        error('%s: this grid is only available for 11345 sources.', ...
            upper(mfilename));
    end
    file = [basepath '/data/spherical_grids/fabian/' filename];
    url = ['https://dev.qu.tu-berlin.de/projects/data/repository/revisions/' ...
        'master/raw/spherical_grids/fabian/' filename];
    % download file if not present
    if ~exist(file,'file')
        download_file(url,file);
    end
    tmp = load(file,'-ascii');
    points = tmp(:,1:3);
    weights = tmp(:,4);
elseif strcmp('gauss',spherical_grid)
    % the number of secondary sources needs to be 2,8,18,32, ... , 
    % see Ahrens (2012)
    if mod(number,sqrt(number/2))~=0
        error(['%s: the number of secondary sources needs to be ', ...
            '2*n^2 for a gauss grid.'],upper(mfilename));
    end
    number = sqrt(number/2);
    % get gauss points and weights
    [p,w] = legpts(number);
    % sampling points along azimuth
    PHI = linspace(0,2*pi,2*number+1);
    % remove the last one, because phi=0 and phi=2pi are the same
    PHI = PHI(1:end-1);
    % sampling points along elevation
    THETA = acos(p)-pi/2;
    % get grid points
    [phi,theta] = meshgrid(PHI,THETA);
    [~,weights] = meshgrid(PHI,w);
    r = ones(size(phi));
    % convert to cartesian
    [points(:,1) points(:,2) points(:,3)] = sph2cart(phi(:),theta(:),r(:));
    % incorporating integration weights
    weights = weights(:).*cos(theta(:));
else
    error(['%s: the given spherical grid is not available, have a look at ' ...
        'http://github.com/sfstoolbox/data for avialable grids.'], ...
        upper(mfilename));
end
