function [xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z,conf)
%XYZ_GRID returns a xyz-grid for the listening area
%
%   Usage: [xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z)
%
%   Input parameters:
%       X        - x-axis / m; single value or [xmin,xmax]
%       Y        - y-axis / m; single value or [ymin,ymax]
%       Z        - z-axis / m; single value or [zmin,zmax]
%       conf     - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       xx,yy,zz - matrices representing the xy-grid / m
%       x,y,z    - x-, y-, z-axis / m
%
%   XYZ_GRID(X,Y,Z) creates a xyz-grid to avoid a loop in the sound field
%   calculation for the whole listening area. It returns also the x-, y-, z-axis
%   for the listening area, defined by the points given with X,Y,Z.
%
%   see also: xyz_axes, xyz_axes_selection, sound_field_mono

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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Computation =====================================================
% creating x-, y-axis
[x,y,z] = xyz_axes(X,Y,Z,conf);
% check which dimensions will be non singleton
dimensions = xyz_axes_selection(x,y,z);
% create xyz-grid
if all(dimensions)
    % create a 3D grid => size(xx)==[resolution resolution resolution]
    [xx,yy,zz] = meshgrid(x,y,z);
elseif dimensions(1) && dimensions(2)
    % create a 2D grid => size(xx)==[resolution resolution]
    [xx,yy] = meshgrid(x,y);
    zz = meshgrid(z,y);
elseif dimensions(1) && dimensions(3)
    [xx,zz] = meshgrid(x,z);
    yy = meshgrid(y,z);
elseif dimensions(2) && dimensions(3)
    [yy,zz] = meshgrid(y,z);
    xx = meshgrid(x,z);
elseif any(dimensions)
    % create a 1D grid => size(xx)==[resolution 1]
    xx = x;
    yy = y;
    zz = z;
else
    % create a 0D grid => size(xx)==[1 1]
    xx = x(1);
    yy = y(1);
    zz = z(1);
end
