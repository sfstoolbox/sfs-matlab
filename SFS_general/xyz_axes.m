function [x,y,z] = xyz_axes(X,Y,Z,conf)
%XYZ_AXES returns the x-, y-, and z-axis for the listening area
%
%   Usage: [x,y,z] = xyz_aex(X,Y,Z,[conf])
%
%   Input parameters:
%       X        - [xmin,xmax]
%       Y        - [ymin,ymax]
%       Z        - [zmin,zmax]
%       conf     - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x,y,z    - x-, y-, z-axis
%
%   XYZ_AXES(X,Y,Z,conf) creates the x-, y-, and -z-axis for the listening area.
%
%   see also: xyz_grid, xyz_axes_selection

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


%% ===== Checking input parameters =======================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
[X,Y,Z] = axis_vector(X,Y,Z);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


% ===== Configuration ====================================================
xysamples = conf.xysamples;


%% ===== Computation =====================================================
% creating x-, y-, and z-axis
x = linspace(X(1),X(2),xysamples);
y = linspace(Y(1),Y(2),xysamples);
z = linspace(Z(1),Z(2),xysamples);
