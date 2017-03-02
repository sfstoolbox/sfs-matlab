function Nse = sphexp_truncation_order(r,f,nmse,conf)
%SPHEXP_TRUNCATION_ORDER yields the bound of summation for a spherical expansion
%of an arbitrary sound field
%
%   Usage: Nse = sphexp_truncation_order(r,f,nmse,conf)
%
%   Input parameters:
%       r           - max 3D distance from expansion center / m
%       f           - frequency / Hz
%       nmse        - maximum bound for normalized mean squared error
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Nse         - maximum order for spherical expansion
%
%   SPHEXP_TRUNCATION_ORDER(r,f,nmse,conf) yields the order up to which
%   a the spherical expansion coefficients of an arbitrary sound field have
%   be summed up. For a given frequency (f) the normalized truncation mean
%   squared error is below the specified error bound (nmse) at any point inside
%   a spherical volume around the expansion center with a radius r. Note, that
%   this approximation holds for any sound field and might heavily over-estimate
%   the needed order in specific cases.
%
%   References:
%       Kennedy et al. (2007) - "Intrinsic Limits of Dimensionality and
%                               Richness in Random Multipath Fields",
%                               IEEE Transactions on Signal Processing
%
%   see also: sphexp_truncation

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
isargstruct(conf);

%% ===== Configuration ==================================================
c = conf.c;

%% ===== Computation ====================================================
% See Kennedy et al. (eq. 42/43)
lambda = c/f;  % wave length
delta = max(0,ceil(log(0.67848/nmse)));  % nmse dependent term
Nse = ceil(pi*r*exp(1)/lambda) + delta;
