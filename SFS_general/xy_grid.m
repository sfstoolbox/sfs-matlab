function [xx,yy,x,y] = xy_grid(X,Y,conf)
%XY_GRID returns a xy-grid for the listening area
%
%   Usage: [xx,yy,x,y] = xy_grid(X,Y,[conf])
%
%   Input parameters:
%       X       - x-axis / m; single value or [xmin,xmax]
%       Y       - y-axis / m; single value or [ymin,ymax]
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       xx,yy   - matrices representing the xy-grid / m
%       x,y     - x-, y-axis / m
%
%   XYGRID(X,Y,conf) creates a xy-grid to avoid a loop in the sound field
%   calculation for the whole listening area. It returns also the x-, y-axis for
%   the listening area, defined by the points given with X,Y.
%
%   see also: sound_field_mono_wfs

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
isargvector(X,Y);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


% ===== Configuration ====================================================
resolution = conf.resolution;


%% ===== Computation =====================================================
% creating x-, y-axis
x = linspace(X(1),X(2),resolution);
y = linspace(Y(1),Y(2),resolution);
% create xy-grid
[xx,yy] = meshgrid(x,y);
