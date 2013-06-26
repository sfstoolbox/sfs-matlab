function [z,p] = sphbesselh_zeros(order)
%SPHBESSELH_ZEROS finds zeros of spherical hankel function
%
%   Usage: [z,p] = sphbesselh_zeros(order)
%
%   Input parameters:
%       order       - order of hankel function
%
%   Output parameters:
%       FIXME: correct the following documentation
%       z       - zeros ...
%       p       - zeros ...
%
%   SPHBESSELH_ZEROS(order) finds zeros for a spherical hankel functino of
%   the specified order.
%
%   see also: sphbesselh, driving_function_imp_nfchoa

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
isargpositivescalar(order);


%% ===== Computation =====================================================
% compute normalized roots/zeros of spherical Hankel function
% FIXME: this works only until a order of 85. For higher orders factorial will
% return Inf. Hence, for higher orders we have to find another way to calculate
% the zeros of the Hankel funtion.
B = zeros(1,order+2);
A = B;
for n=0:order
    B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B = B(end:-1:1);
% find zeros/roots
z = roots(B);
% find roots (for another kind of hankel function ??? This is needed for a plane
% wave in NFCHOA)
A(2) = 1;
p = roots(A);
