function [ Ynm ] = sphharmonics(n,m,beta,alpha);
% SPHHARMONICS spherical harmonics function
%
%   Usage: Ynm = sphbesselh(n,m,beta,alpha)
%
%   Input parameters:
%       n     - spherical harmonic degree
%       m     - spherical harmonic order
%       beta  - colatitude to be calculated
%       alpha - azimuth to be calculated
%
%   Output parameters:
%       Ynm   - values of spherical harmonics function
%
%   SPHHARMONICS(n,m,alpha,beta) spherical harmonics function of degree n and
%   order m for the angles alpha, beta.
%   alpha and beta can be arrays but have to be of same size or one of them
%   has to be a scalar.
%
%   see also: sphbesselj, sphbessely, sphbesselh, asslegendre

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
% Copyright (c) 2012      Jens Ahrens                                        *
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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargpositivescalar(n)
isargscalar(m)
isargnumeric(alpha,beta)
if n<abs(m)
    warning( 'Absolute value of order m must be less than or equal to the degree n.' ); 
    Ynm = zeros(size(alpha));
    return;
end


%% ===== Computation =====================================================
Lnm = asslegendre( n, abs(m), cos( beta ) );
factor_1 = ( 2*n + 1 ) / ( 4*pi );
factor_2 = factorial( n - abs(m) ) ./ factorial( n + abs(m) );
Ynm = (-1).^m .* sqrt( factor_1 .* factor_2 ) .* Lnm .* exp( 1i .* m .* alpha );
