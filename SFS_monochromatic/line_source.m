function S = line_source(x,y,xs,f,conf)
%LINE_SOURCE returns the Green's function for a line source
%
%   Usage: S = line_source(x,y,xs,f,[conf])
%
%   Input options:
%       x,y      - x,y points for which the Green's function should be calculated
%       xs       - position of the line source
%       f        - frequency of the line source
%       conf     - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       S        - Wave field of a line source located at xs
%
%   LINE_SOURCE(x,y,xs,f) calculates the wave field of a line source
%   located at xs for the given points x,y and the frequency f. The wave
%   field is calculated by the Greens function.
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: point_source

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
xs = position_vector(xs);
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
% Source model for a line source: 2D Green's function.
%
%              i   (2)/ w        \
% G(x-xs,w) =  -  H0  | - |x-xs| |
%              4      \ c        /
%
% see: Williams1999, p. 266
%
S = 1i/4 * besselh(0,2,omega/c* ...
    sqrt( (x-xs(1)).^2 + (y-xs(2)).^2 ));
