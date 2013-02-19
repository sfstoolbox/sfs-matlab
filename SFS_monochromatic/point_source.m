function S = point_source(x,y,xs,f,conf)
%POINT_SOURCE returns the Green's function for a point source
%
%   Usage: S = point_source(x,y,xs,f,[conf])
%
%   Input options:
%       x,y     - x,y points for which the Green's function should be calculated
%       xs      - position of the point source
%       f       - frequency of the point source
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       S       - Wave field of a point source located at x0,y0
%
%   POINT_SOURCE(x,y,xs,f) calculates the wave field of a point source
%   located at xs for the given points x,y and the frequency f. The wave
%   field is calculated by the Greens function.
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: line_source

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


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargmatrix(x,y);
isargposition(xs);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
omega = 2*pi*f;
% Source model for a point source: 3D Green's function.
%
%              1  e^(-i w/c |x-xs|)
% G(x-xs,w) = --- -----------------
%             4pi      |x-xs|
%
% see: Williams1999, p. 198
%
S = 1/(4*pi) * exp(-1i*omega/c.*sqrt((x-xs(1)).^2+(y-xs(2)).^2)) ./ ...
        sqrt((x-xs(1)).^2+(y-xs(2)).^2);
