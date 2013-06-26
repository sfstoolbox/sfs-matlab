function N = nfchoa_order(nls)
%NFCHOA_ORDER returns the maximum order for the spherical harmonics for the
%given number of secondary sources
%
%   Usage: N = nfchoa_order(nls)
%
%   Input parameters:
%       nls     - number of secondary sources
%
%   Output parameters:
%       N       - spherical harmonics order
%
%   NFCHOA_ORDER(nls) returns the maximum order of spherical harmonics for the
%   given number of loudspeakers. This is calculated after ... with N = nls/2.
%
%   References:
%       FIXME: add one
%
%   see also: driving_function_imp_nfchoa, driving_function_mono_nfchoa

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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargpositivescalar(nls)


%% ===== Computation =====================================================
% get maximum order of spherical harmonics
if isodd(nls)
    N = (nls-1)/2;
else
    N = floor((nls+1)/2);
end
