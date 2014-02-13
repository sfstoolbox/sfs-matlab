function out = sphbesselj(nu,z)
% SPHBESSELJ spherical bessel function of first kind of order nu, and argument z
%
%   Usage: out = sphbesselj(nu,z)
%
%   Input parameters:
%       nu  - order of bessel function
%       z   - argument of bessel function
%
%   Output parameters:
%       out - value of bessel function at point z
%
%   SPHBESSELJ(nu,z) spherical bessel function of order nu, frist type, and
%   argument z
%
%   see also: sphbesselh, sphbessely

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargscalar(nu)
isargnumeric(z)


%% ===== Computation =====================================================
out = zeros(size(z));

% avoid division by "0"
if (nu==0)
    out(z==0) = 1;
elseif (nu~=0)
    out(z==0) = 0;
end

% finally evaluate for z~=0
out(z~=0) = sqrt(pi./(2.*z(z~=0))) .* besselj(nu+0.5, z(z~=0));
