function out = besselh_derived(nu,k,z)
% BESSELH_DERIVED derivative of cylindrical hankel function of k kind of order nu, and argument z
%
%   Usage: out = besselh_derived(nu,k,z)
%
%   Input parameters:
%       nu  - order of bessel function
%       z   - argument of bessel function
%
%   Output parameters:
%       out - value of bessel function at point z
%
%   BESSELH_DERIVED(nu,z) derivation of cylindrical hankel function of 
%   order nu, k kind, and argument z
%
%   References:
%       (4.1-51) in Ziomek (1995) - "Fundamentals of acoustic field theory 
%                                   and space-time signal processing"
%
%   see also: besselh

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
isargscalar(k)
isargnumeric(z)


%% ===== Computation =====================================================
out = 0.5*(besselh(nu-1,k,z) - besselh(nu+1,k,z)); 