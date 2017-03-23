function out = sphbesselj(nu,z)
%SPHBESSELJ spherical bessel function of first kind of order nu, and argument z
%
%   Usage: out = sphbesselj(nu,z)
%
%   Input parameters:
%       nu  - order of bessel function
%       z   - argument of bessel function
%
%   Output parameters:
%       out - value of bessel function at point z
%
%   SPHBESSELJ(nu,z) spherical bessel function of order nu, frist type, and
%   argument z
%
%   See also: sphbesselh, sphbessely

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargscalar(nu)
isargnumeric(z)


%% ===== Computation =====================================================
out = zeros(size(z));

% Avoid division by "0"
if (nu==0)
    out(z==0) = 1;
elseif (nu~=0)
    out(z==0) = 0;
end

% Finally evaluate for z~=0
out(z~=0) = sqrt(pi./(2.*z(z~=0))) .* besselj(nu+0.5, z(z~=0));
