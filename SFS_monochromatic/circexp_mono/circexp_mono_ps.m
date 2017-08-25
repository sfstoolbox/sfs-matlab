function Pm = circexp_mono_ps(xs,Nce,f,xq,conf)
%CIRCEXP_MONO_PS calculates the circular basis expansion of a point source
%
%   Usage: Pm = circexp_mono_ps(xs,Nce,f,xq,conf)
%
%   Input parameters:
%       xs      - position of point source / m [1 x 3]
%       Nce     - maximum order of circular basis expansion
%       xq      - optional expansion center / m [1 x 3]
%       f       - frequency of the monochromatic source / Hz
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       Pm      - regular circular expansion coefficients
%                 for m = 0:Nce, [1 x Nce+1]
%
%   CIRCEXP_MONO_PS(xs,Nce,f,xq,conf) returns the circular basis expansion of
%   a point source.
%
%   See also: circexp_mono_pw

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


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargcoord(xq,xs);
isargpositivescalar(Nce,f);


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
xs = xs - xq;  % shift coordinates
[phis, rs] = cart2pol(xs(1),xs(2));
k = 2*pi*f/c;

%                       j^(|m|-m)
% P_m = -jk h_|m|(k rs) --------- e^(-im phis)
%                          4pi
Pm = zeros(1,2*Nce+1);
for m=-Nce:Nce
  Pm(m+Nce+1) = ...
      -1i*k*sphbesselh(abs(m),2,k*rs).*1i.^(abs(m)-m).*exp(-1i*m*phis);
end
Pm = Pm./(4*pi);
