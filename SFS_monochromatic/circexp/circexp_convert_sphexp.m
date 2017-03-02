function Am = circexp_convert_sphexp(Anm)
%CIRCEXP_CONVERT_SPHEXP converts regular spherical expansion coefficients into
%regular circular expansion coefficients
%
%   Usage: Am = circexp_convert_sphexp(Anm)
%
%   Input parameters:
%       Anm           - regular spherical expansion coefficients
%
%   Output parameters:
%       Am            - regular circular expansion coefficients
%
%   References:
%       Hahn, Winter, Spors (2016) -
%           "Local Wave Field Synthesis by Spatial Band-limitation in the
%            Circular/Spherical Harmonics Domain", 140th AES Convention

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

%% ===== Checking of input parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargvector(Anm);
Nse = sqrt(size(Anm,1))-1;

%% ===== Computation ====================================================

% Implementation of Hahn2016, Eq. (32)
Am = zeros(2*Nse+1,size(Anm,2));
for m=-Nse:Nse
    v = sphexp_index(m);  % (n,m) = (abs(m),m);
    Am(m+Nse+1,:) = Anm(v,:)./(4*pi.*1j.^(m-abs(m))* ...
                     sphharmonics(abs(m),-m,0,0));
end
