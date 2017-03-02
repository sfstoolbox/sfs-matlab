function Anm = sphexp_convert_circexp(Am)
%SPHEXP_CONVERT_CIRCEXP compute converts regular circular expansion into
%spherical expansion coefficients
%
%   Usage: Anm = sphexp_convert_circexp(Am)
%
%   Input parameters:
%       Am            - regular circular expansion coefficients
%
%   Output parameters:
%       Anm           - regular spherical expansion coefficients
%
%   References:
%       Hahn, Spors (2015) - "Sound Field Synthesis of Virtual Cylindrical
%                            Waves using Circular and Spherical Loudspeaker
%                            Arrays", 138th AES Convention

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
isargvector(Am);
Nce = (length(Am)-1)/2;

%% ===== Computation ====================================================

% Implementation of Hahn2015, Eq. (14)
Anm = zeros((Nce+1).^2,1);
for m=-Nce:Nce
    % for theta=0 the legendre polynom is zero if n+m is odd
    for n = abs(m):2:Nce
        v = sphexp_index(m,n);
        Anm(v) = 4*pi.*1j.^(m-n).*sphharmonics(n,-m,0,0).*Am(m+Nce+1);
    end
end
