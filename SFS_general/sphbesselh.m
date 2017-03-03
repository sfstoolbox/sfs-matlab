function out = sphbesselh(nu,k,z)
% SPHBESSELH spherical hankel function of order nu, kind k, and argument z
%
%   Usage: out = sphbesselh(nu,k,z)
%
%   Input parameters:
%       nu  - order of hankel function
%       k   - kind of hankel function (1 ^= first kind, 2 ^= second kind)
%       z   - argument of hankel function
%
%   Output parameters:
%       out - value of hankel function at point z
%
%   SPHBESSELH(nu,k,z) spherical hankel function of order nu, kind k, and
%   argument z
%
%   See also: sphbesselj, sphbessely

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
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargscalar(nu)
if (k==1)
    sign = 1;
elseif (k==2)
    sign = -1;
else
    error(['%s: Invalid kind of Hankel function is asked ',...
           '(k has to be 1 or 2).'],upper(mfilename));
end
isargnumeric(z)


%% ===== Computation =====================================================
out = sphbesselj(nu, z) + 1j .* sign .* sphbessely(nu, z);
