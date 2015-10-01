function [x,y,z] = xyz_axes(X,Y,Z,conf)
%XYZ_AXES returns the x-, y-, and z-axis for the listening area
%
%   Usage: [x,y,z] = xyz_aex(X,Y,Z,[conf])
%
%   Input parameters:
%       X        - x-axis / m; single value or [xmin,xmax]
%       Y        - y-axis / m; single value or [ymin,ymax]
%       Z        - z-axis / m; single value or [zmin,zmax]
%       conf     - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x,y,z    - x-, y-, z-axis / m
%
%   XYZ_AXES(X,Y,Z,conf) creates the x-, y-, and -z-axis for the listening area.
%
%   See also: xyz_grid, xyz_axes_selection

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


% ===== Configuration ====================================================
resolution = conf.resolution;


%% ===== Computation =====================================================
x = X; y = Y; z = Z;
dimensions = xyz_axes_selection(X,Y,Z);
if dimensions(1), x = linspace(x(1),x(2),resolution)'; end
if dimensions(2), y = linspace(y(1),y(2),resolution)'; end
if dimensions(3), z = linspace(z(1),z(2),resolution)'; end
