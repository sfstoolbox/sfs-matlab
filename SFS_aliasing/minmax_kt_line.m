function [kGtmin, kGtmax] = minmax_kt_line(x0, xl, Ll, alphal)
%

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2018 SFS Toolbox Developers                             *
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

%% ===== Main ============================================================

% intersection of:
%            g1: x = xl + gamma * nl     -Ll/2 <= gamma <= Ll/2
%  g2 (n0-axis): 0 = (x - x0)^T * n0
%
% -> gamma = ((x0 - xl)^T * n0) / (nl^T * n0)

xl = bsxfun(@minus, xl, x0(:,1:3));  % shift xl about x0
nl = [cos(alphal), sin(alphal), 0];  % direction vector of line
n0 = x0(:,4:6);
phin0 = cart2pol(n0(:,1),n0(:,2));

gamma = -(xl(:,1).*n0(:,1) + xl(:,2).*n0(:,2))./(n0*nl.');

select = gamma <= 0 | gamma > Ll/2;
gamma1(select) = +Ll/2;
gamma1(~select) = gamma;

select = gamma < -Ll/2 | gamma > 0;
gamma2(select) = -Ll/2;
gamma2(~select) = gamma;

select = ~isinf(gamma1);
if any(select)
  x1(select,:) = bsxfun(@plus, xl(select,:), gamma1(select).*nl);
elseif any(~select)
  x1(~select,:) = sign(gamma1(~select))*nl;
end

select = ~isinf(gamma2);
if any(select)
  x2(select,:) = bsxfun(@plus, xl(select,:), gamma2(select).*nl);
elseif any(~select)
  x2(~select,:) = sign(gamma2(~select))*nl;
end

if abs(x1) < 1e-10
  x1 = x2;
  degenerated = true;
else
  degenerated = false;
end
if abs(x2) < 1e-10
  x2 = x1;
  if degenerated
    error('degenerated case');
  end
end

phi1 = cart2pol(x1(:,1),x1(:,2));
phi2 = cart2pol(x2(:,1),x2(:,2));

kGt1 = sin(phin0 - phi1);
kGt2 = sin(phin0 - phi2);

kGtmin = min(kGt1,kGt2);
kGtmax = max(kGt1,kGt2);
