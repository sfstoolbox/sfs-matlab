function out = sphbesselh(nu,k,z)
% SPHBESSELH spherical hankel function of order nu, kind k, and argument z
%
%   Usage: out = sphbesselh(nu,k,z)
%
%   Input parameters:
%       nu  - order of hankel function
%       k   - kind of hankel function (1 ^= first kind, 2 ^= second kind)
%       z   - argument of hankel function
%
%   Output parameters:
%       out - value of hankel function at point z
%
%   SPHBESSELH(nu,k,z) spherical hankel function of order nu, kind k, and
%   argument z
%
%   see also: sphbesselj, sphbessely

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
nargmax = 3;
narginchk(nargmin,nargmax);
isargscalar(nu)
if (k==1)
    sign = 1;
elseif (k==2)
    sign = -1;
else
    error(['%s: Invalid kind of Hankel function is asked ',...
           '(k has to be 1 or 2).'],upper(mfilename));
end
isargnumeric(z)


%% ===== Computation =====================================================
out = sphbesselj(nu, z) + 1j .* sign .* sphbessely(nu, z);
