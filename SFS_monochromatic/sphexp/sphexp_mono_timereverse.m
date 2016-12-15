function ABnm = sphexp_mono_timereverse(ABnm)
%SPHEXP_MONO_TIMEREVERSE compute coefficients of a time reversed sound field
%
%   Usage: ABnm = sphexp_mono_timereverse(ABnm)
%
%   Input parameters:
%       ABnm          - spherical expansion coefficients of sound field
%
%   Output parameters:
%       ABnm          - spherical expansion coefficients of time reversed
%                       sound field

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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

%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargvector(ABnm);
L = length(ABnm);
isargsquaredinteger(L);

%% ===== Computation ====================================================

for n=0:(sqrt(L)-1)
    v = sphexp_index(-n:n,n);
    ABnm(v) = conj(sphexp_access(ABnm,n:-1:-n,n));
end
