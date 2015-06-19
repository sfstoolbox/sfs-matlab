function Anm = sphexp_truncate(Pnm, N, Nshift)
% Truncates spherical expansion by setting remaining coefficients to zero
%
%   Usage: Anm = sphexp_truncate(Pnm, N, Nshift)
%
%   Input parameters:
%       Pnm         - 1D array of spherical expansion coefficients [n x Nf]
%       N           - maximum order of spherical expansion
%       Nshift      - 
%
%   Output parameters:
%       Anm         - 1D array of bandlimited spherical expansion
%                     coefficients [N x Nf]
%
%   SPHEXP_TRUNCATE(Pnm, N, Nshift) sets coefficients belonging to an order
%   higher than N to zero.
%
%   see also: sphexp_access sphexp_truncation_order

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
nargmax = 3;
narginchk(nargmin,nargmax);
isargmatrix(Pnm);
isargpositivescalar(N);
if nargin == nargmin
  Nshift = 0;
else
  isargscalar(Nshift);
end
%% ===== Variable =======================================================
Nse = sqrt(size(Pnm, 1))-1; 

%% ===== Computation ====================================================
Anm = zeros(size(Pnm));

for m=max(-N+Nshift,-Nse):min(N+Nshift,Nse)
  v = sphexp_index(m,abs(m):Nse);
  Anm(v,:) = Pnm(v,:);
end

end

