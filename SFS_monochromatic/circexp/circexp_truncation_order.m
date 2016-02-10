function Nce = circexp_truncation_order(r, f, nmse, conf)
%CIRCEXP_TRUNCATION_ORDER computes truncation order for circular expansion 
%coefficients of an arbitrary sound field
%
%   Usage: Nce = circexp_truncation_order(r, f, nmse, conf)
%
%   Input parameters:
%       r           - max 2D distance from expansion center / m
%       f           - frequency / Hz
%       nmse        - maximum bound for normalized mean squared error
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Nce         - Maximum order for circular expansion
%
%   CIRCEXP_TRUNCATION_ORDER(r, f, epsilon, conf) yields the order up to which 
%   a the circular expansion coefficients of an arbitrary sound field have
%   be summed up. For a given frequency and maximum radius the normalized 
%   truncation mean squared error is below the specified error bound (nmse).
%
%   References:
%       Kennedy et al. (2007) - "Intrinsic Limits of Dimensionality and 
%                               Richness in Random Multipath Fields",
%                               IEEE Transactions on Signal Processing

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargpositivescalar(f,r,nmse);

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Computation ====================================================
% See Kennedy et al. (eq. 36)
lambda = c/f;  % wave length
delta = max(0, ceil(0.5*log(0.0093/nmse)));
Nce = ceil(pi*r*exp(1)/lambda) + delta;

end

