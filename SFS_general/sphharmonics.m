function [ Ynm ] = sphharmonics(n,m,theta,phi);
% SPHHARMONICS spherical harmonics function
%
%   Usage: Ynm = sphharmonics(n,m,theta,phi)
%
%   Input parameters:
%       n     - spherical harmonic degree
%       m     - spherical harmonic order
%       theta - elevation angle
%       phi   - azimuth angle
%
%   Output parameters:
%       Ynm   - values of spherical harmonics function
%
%   SPHHARMONICS(n,m,theta,phi) spherical harmonics function of degree n and
%   order m for the angles phi, theta.
%   phi and theta can be arrays but have to be of same size or one of them
%   has to be a scalar.
%
%   see also: sphbesselj, sphbessely, sphbesselh, asslegendre

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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking input parameters =======================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargpositivescalar(n)
isargscalar(m)
isargnumeric(phi,theta)
if n<abs(m)
    warning( 'Absolute value of order m must be less than or equal to the degree n.' ); 
    Ynm = zeros(size(phi));
    return;
end
phi = correct_azimuth(phi);
theta = correct_elevation(theta);


%% ===== Computation =====================================================
% NOTE: we use here sin(theta) opposite to cos(theta) as presented in Ahrens (2012),
% because we have another coordinate system convention -- compare Ahrens, Fig.A.1
Lnm = asslegendre(n,abs(m),sin(theta));
factor_1 = (2*n+1) / (4*pi);
factor_2 = factorial(n-abs(m)) ./ factorial(n+abs(m));
Ynm = (-1).^m .* sqrt(factor_1.*factor_2) .* Lnm .* exp(1i.*m.*phi);
