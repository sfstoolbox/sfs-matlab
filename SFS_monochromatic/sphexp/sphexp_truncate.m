function Anm = sphexp_truncate(Pnm, N, M, Mshift)
% Truncates spherical expansion by setting remaining coefficients to zero
%
%   Usage: Anm = sphexp_truncate(Pnm, N, M, Mshift)
%
%   Input parameters:
%       Pnm         - 1D array of spherical expansion coefficients [n x Nf]
%       N           - maximum degree of spherical expansion
%       M           - maximum order of spherical expansion
%       Mshift      - shift for asymmetric trunction with respect to order
%
%   Output parameters:
%       Anm         - 1D array of bandlimited spherical expansion
%                     coefficients [N x Nf]
%
%   SPHEXP_TRUNCATE(Pnm, N, M, Mshift) sets coefficients belonging to an degree
%   n higher than N or an order m exceeding: -M+Mshift to +M+Mshift to zero.
%
%   see also: sphexp_truncation_order

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

%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(Pnm);
isargpositivescalar(N);
switch nargin
  case 2
    M = N;
    Mshift = 0;
  case 3
    isargpositivescalar(M);
    Mshift = 0;
  case 4
    isargscalar(Mshift);
end
%% ===== Variable =======================================================
Nse = min(sqrt(size(Pnm, 1))-1, N);

%% ===== Computation ====================================================
Anm = zeros(size(Pnm));

for m=max(-M+Mshift,-Nse):min(M+Mshift,Nse)
  v = sphexp_index(m,abs(m):Nse);
  Anm(v,:) = Pnm(v,:);
end

end

