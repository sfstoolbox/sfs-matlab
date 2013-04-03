function g = greens_function_imp(x,y,xs,src,conf)
%GREENS_FUNCTION_IMP returns a Green's function in the time domain
%
%   Usage: g = greens_function_imp(x,y,xs,src,[conf])
%
%   Input options:
%       x,y     - x,y points for which the Green's function should be calculated
%       xs      - position of the source
%       src     - source model of the Green's function. Valid models are:
%                   'ps' - point source
%                   'ls' - line source
%                   'pw' - plane wave
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       g       - Green's function evaluated at the points x,y
%
%   GREENS_FUNCTION_IMP(x,y,xs,src) calculates the Green's function for the
%   given source model located at xs for the given points x,y. Note, that the
%   delta function for the time t is not performed and the result is independent
%   of t.
%
%   References:
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: greens_function_mono, wave_field_mono

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


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargmatrix(x,y);
isargposition(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Computation =====================================================
% calculate Green's function for the given source model
if strcmp('ps',src)
    % Source model for a point source: 3D Green's function.
    %
    %              1  delta(t - |x-xs|/c)
    % g(x-xs,t) = --- -------------------
    %             4pi       |x-xs|
    %
    % see: Williams1999, p. FIXME: ??
    %
    g = 1./(4*pi) ./ sqrt((x-xs(1)).^2+(y-xs(2)).^2);

elseif strcmp('ls',src)
    % Source model for a line source: 2D Green's function.
    %
    %              
    % g(x-xs,t) =  
    %              
    %
    % see: Williams1999, p. FIXME
    %
    to_be_implemented;

elseif strcmp('pw',src)
    % Source model for a plane wave:
    %
    % g(x,t) = delta(t - nx/c)
    %
    % see: Williams1999, p. FIXME
    %
    % direction of plane wave
    g = 1;

else
    error('%s: %s is not a valid source model for the Green''s function', ...
        upper(mfilename),src);
end
