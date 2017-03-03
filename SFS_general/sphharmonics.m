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
%   See also: sphbesselj, sphbessely, sphbesselh, asslegendre

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
