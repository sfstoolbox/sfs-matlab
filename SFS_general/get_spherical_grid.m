function [points,weights] = get_spherical_grid(number,conf)
%GET_SPHERICAL_GRID returns grid points and weights
%
%   Usage: [points,weights] = get_spherical_grid(number,conf)
%
%   Input parameters:
%       number  - number of grid points
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       points  - grid points
%       weights - integration weights for the grid points
%
%   GET_SPHERICAL_GRID(number,conf) returns the points and weights for a grid on
%   a sphere. The type of grid is specified by conf.secondary_sources.grid.
%   For available grids, have a look at http://github.com/sfstoolbox/data.
%   It expects the grid files at SFS_basepath/data/spherical_grids. If the
%   desired file is not available on the hard disk, the function tries to
%   download it directly from github.
%   For conf.secondary_sources.grid='gauss' the grid positions are calculated
%   after Ahrens (2012), p. 121 (see also: Rafaely (2015), p. 64)
%
%   References:
%       J. Ahrens (2012) - "Analytic Methods of Sound Field Synthesis", Springer.
%       B. Rafaely (2015) - "Fundamentals of Spherical Array Processing", Springer.
%
%   See also: secondary_source_positions,
%       weights_for_points_on_a_sphere_rectangle

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(number);
isargstruct(conf);


%% ===== Configuration ===================================================
spherical_grid = conf.secondary_sources.grid;


%% ===== Main ============================================================
filename = sprintf('%06.0fpoints.mat',number);
basepath = get_sfs_path();

if strcmp('equally_spaced_points',spherical_grid)
    % Check if we have a squared number of points, because only for those values
    % are equally spaced points grids available.
    if mod(number,sqrt(number))~=0
        error('%s: number has to be a squared number.',upper(mfilename));
    end
    file = [basepath '/data/spherical_grids/equally_spaced_points/' filename];
    url = ['https://raw.githubusercontent.com/sfstoolbox/data/master/' ...
           'spherical_grids/equally_spaced_points/' filename];
    % Download file if not present
    if ~exist(file,'file')
        download_file(url,file);
    end
    tmp = load(file,'-ascii');
    points = tmp(:,1:3);
    weights = tmp(:,4);
elseif strcmp('fabian',spherical_grid)
    % Here we have only one number of secondary sources available
    if number~=11345
        error('%s: this grid is only available for 11345 sources.', ...
            upper(mfilename));
    end
    file = [basepath '/data/spherical_grids/fabian/' filename];
    url = ['https://raw.githubusercontent.com/sfstoolbox/data/master/' ...
           'spherical_grids/fabian/' filename];
    % Download file if not present
    if ~exist(file,'file')
        download_file(url,file);
    end
    tmp = load(file,'-ascii');
    points = tmp(:,1:3);
    weights = tmp(:,4);
elseif strcmp('gauss',spherical_grid)
    % The number of secondary sources needs to be 2,8,18,32, ... ,
    % see Ahrens (2012)
    if mod(number,sqrt(number/2))~=0
        error(['%s: the number of secondary sources needs to be ', ...
            '2*n^2 for a gauss grid.'],upper(mfilename));
    end
    number = sqrt(number/2);
    % Get gauss points and weights
    [p,w] = legpts(number);
    % Sampling points along azimuth
    PHI = linspace(0,2*pi,2*number+1);
    % Remove the last one, because phi=0 and phi=2pi are the same
    PHI = PHI(1:end-1);
    % Sampling points along elevation
    THETA = acos(p)-pi/2;
    % Get grid points
    [phi,theta] = meshgrid(PHI,THETA);
    [~,weights] = meshgrid(PHI,w);
    r = ones(size(phi));
    % Convert to cartesian
    [points(:,1) points(:,2) points(:,3)] = sph2cart(phi(:),theta(:),r(:));
    % Normalize integration weights
    weights = weights(:)*pi/number;
else
    error(['%s: the given spherical grid is not available, have a look at ' ...
        'http://github.com/sfstoolbox/data for avialable grids.'], ...
        upper(mfilename));
end
