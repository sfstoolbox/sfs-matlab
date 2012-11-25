function [x1,x2,xref] = active_dimensions(x,y,z,conf)
%ACTIVE_DIMENSIONS returns the first two active dimensions and the
%corresponding xref values
%
%   Usage: [x1,x2,xref] = active_dimensions(x,y,z,[conf])
%
%   Input options:
%       x,y,z   - vectors conatining the x-, y- and z-axis values
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options:
%       x1      - first active dimension (this could be the x or y axis)
%       x2      - second active dimension (this could be the y or z axis)
%       xref    - corresponding xref values (dim: 1x2)
%
%   ACTIVE_DIMENSIONS(x,y,z,conf) returns the first two active dimensions.
%   An active dimension is any dimension of x,y,z for which the first and
%   the last value of the axis vector is not the same. In addition the
%   coresponding xref values are returned. For example, if the first two 
%   active dimensions are y,z then xref = [conf.xref(2) conf.xref(3)].
%
%   see also: norm_wavefield, plot_wavefield

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


%% ===== Checking of input parameters ====================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
isargvector(x,y,z);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
xref = position_vector(conf.xref);


%% ===== Computation =====================================================
% Check if we have any inactive dimensions
if x(1)==x(end)
    x1 = y;
    x2 = z;
    xref = [xref(2) xref(3)];
elseif y(1)==y(end)
    x1 = x;
    x2 = z;
    xref = [xref(1) xref(3)];
else
    x1 = x;
    x2 = y;
    xref = [xref(1) xref(2)];
end
