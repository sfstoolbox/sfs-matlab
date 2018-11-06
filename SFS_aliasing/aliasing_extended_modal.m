function f = aliasing_extended_modal(x0, kSx0, x, minmax_kGt_fun, xc, M, conf)
%ALIASING_EXTENDED_MODAL aliasing frequency for an extended listening area for 
%an circular control area at xc with R=M/k where synthesis is focused on.
%
%   Usage: f = aliasing_extended_control(x0, kSx0, x, minmax_kGt_fun, minmax_kSt_fun, conf)
%
%   Input options:
%       x0              - position, direction, and sampling distance of 
%                         secondary sources [N0x7] / m
%       kSx0            - normalised local wavenumber vector of virtual sound 
%                         field [N0x3]
%       x               - position for which aliasing frequency is calculated
%                         [Nx3]
%       minmax_kGt_fun  - function handle to determine the extremal value of 
%                         the tangential component of k_G(x-x0)
%                         [kGtmin, kGtmax] = minmax_kGt_fun(x0)
%       xc              - center of circular control area
%       M               - modal order which defines the radius R=M/k
%       conf            - configuration struct (see SFS_config)
%
%   Output parameters:
%       f   - aliasing frequency [Nx1]
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

phik = cart2pol(kSx0(:,1),kSx0(:,2));  % azimuth angle of kS(x0)
phin0 = cart2pol(x0(:,4),x0(:,5));  % azimuth angle of normal vector n0(x0)

% secondary source selection
select = cos(phin0 - phik) >= 0;  
x0 = x0(select,:);
% k_St(x0) (tangential component of k_S(x0) )
kSt = sin(phin0(select) - phik(select));  

% sampling distance
deltax0 = abs(x0(:,7));

f = inf(size(x,1), 1);
for xdx = 1:size(x,1);
    % mininum and maximum values of k_Gt(x - x_0) 
    % (tangential component of k_G(x-x0))
    [kGtmin, kGtmax] = minmax_kGt_fun(x0,x(xdx,:));
    % aliasing frequency for x0
    f0 = conf.c./(deltax0.*max(abs(kSt-kGtmin),abs(kSt-kGtmax)));
    % radius of local region at f0
    rM = M.*conf.c./(2.*pi.*f0);
    % mininum and maximum values of kSt(x_0) at f0
    [kStmin, kStmax] = minmax_kt_circle(x0, xc, rM);
    select = kSt > kStmin & kSt < kStmax;
    if sum(select) ~= 0
        f(xdx) = min(f0(select));
    end
end
