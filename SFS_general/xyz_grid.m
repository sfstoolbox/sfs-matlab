function [xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z,conf)
%XY_GRID returns a xy-grid for the listening area
%
%   Usage: [xx,yy,zz,x,y,z] = xyz_grid(X,Y,Z,[conf])
%
%   Input parameters:
%       X        - [xmin,xmax]
%       Y        - [ymin,ymax]
%       Z        - [zmin,zmax]
%       conf     - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       xx,yy,zz - matrices representing the xy-grid
%       x,y,z    - x-, y-, z-axis
%
%   XYZ_GRID(X,Y,Z,conf) creates a xyz-grid to avoid a loop in the wave field
%   calculation for the whole listening area. It returns also the x-, y-, z-axis
%   for the listening area, defined by the points given with X,Y,Z.
%
%   see also: wave_field_mono_wfs_25d

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
isargvector(X,Y);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


% ===== Configuration ====================================================
xysamples = conf.xysamples;
zreferenceaxis = conf.zreferenceaxis;


%% ===== Computation =====================================================
% creating x-, y-axis
x = linspace(X(1),X(2),xysamples);
y = linspace(Y(1),Y(2),xysamples);
z = linspace(Z(1),Z(2),xysamples);
% check which dimensions will be non singleton
dimensions = active_dimensions(x,y,z);
% create xyz-grid
if dimensions(1) && dimensions(2)
    [xx,yy] = meshgrid(x,y);
    % create the z grid regarding its reference axis.
    % this means that z changes its values along this particular axis.
    if strcmp('y',zreferenceaxis)
        [~,zz] = meshgrid(x,z);
    elseif strcmp('x',zreferenceaxis)
        zz = meshgrid(z,y);
    else
        error('%s: zreferenceaxis has to be ''y'' or ''x''.',upper(mfilename));
    end
else
    [yy,zz] = meshgrid(y,z);
    xx = meshgrid(x,y);
end

