function [kGtmin, kGtmax] = minmax_kt_circle(x0, xc, Rc)
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
phin0 = cart2pol(x0(:,4),x0(:,5));

% range for k_Gt(x-x_0) (tangential component of k_G)
xc = bsxfun(@minus, xc, x0(:,1:3));
[phil, rc] = cart2pol(xc(:,1),xc(:,2));
klt = sin(phin0 - phil);
kln = sqrt(1 - klt.^2);

rho = Rc./rc;

kGtmin = zeros(size(xc,1),1);
select = rho > 1 | -sqrt(1-rho.^2) > klt;
kGtmin(select) = -1;
kGtmin(~select) = klt(~select).*sqrt(1-rho(~select).^2) ...
    - kln(~select).*rho(~select);

kGtmax = zeros(size(xc,1),1);
select = rho > 1 | sqrt(1-rho.^2) < klt;
kGtmax(select) = +1;
kGtmax(~select) = klt(~select).*sqrt(1-rho(~select).^2) ...
    + kln(~select).*rho(~select);
