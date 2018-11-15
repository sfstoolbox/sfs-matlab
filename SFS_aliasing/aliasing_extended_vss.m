function f = aliasing_extended_vss(x0, xv, kSxv,  x, minmax_kGt_fun, conf)
%ALIASING_EXTENDED_CONTROL aliasing frequency for an extended listening area 
%with an defined control area where the sound field synthesis is prioritized
%
%   Usage: f = aliasing_extended_control(x0, kSx0, x, minmax_kGt_fun, minmax_kSt_fun, conf)
%
%   Input options:
%       x0              - position, direction, and sampling distance of 
%                         secondary sources [N0x7] / m
%       kSx0            - normalised local wavenumber vector kS(x0) 
%                         of virtual sound field at x0 [N0x3]
%       x               - position for which aliasing frequency is calculated
%                         [Nx3]
%       minmax_kGt_fun  - function handle to determine the extremal value of 
%                         the tangential component of k_G(x-x0)
%                         [kGtmin, kGtmax] = minmax_kGt_fun(x0,x)
%       minmax_kSt_fun  - function handle to determine the extremal value of 
%                         the tangential component of k_S(x0)
%                         [kStmin, kStmax] = minmax_kSt_fun(x0)
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

phik = cart2pol(kSxv(:,1),kSxv(:,2));  % azimuth angle of kS(xv)
phinv = cart2pol(xv(:,4),xv(:,5));  % azimuth angle of normal vector nv(xv)
phin0 = cart2pol(x0(:,4),x0(:,5));  % azimuth angle of normal vector n0(x0)

% selection of virtual secondary source 
select = cos(phinv - phik) >= 0;  
xv = xv(select,:);
% k_S,tv(xv) (tangential component of kS(xv) at xv )
kStv = sin(phinv(select) - phik(select));  

% sampling distance
deltax0 = abs(x0(:,7));
deltaxv = abs(xv(:,7));

f = inf(size(x,1), 1);

eps = 0.01;
for xdx = 1:size(x,1);
  
    % mininum and maximum values of k_G,tv(x-xv) 
    % (tangential component of k_G(x-xv) at xv)
    [kGtvmin, kGtvmax] = minmax_kGt_fun(xv,x(xdx,:));
    lambda = max(deltaxv.*max(abs(kStv-kGtvmin),abs(kStv-kGtvmax)));
    
    % mininum and maximum values of k_G,t0(x-x0) 
    % (tangential component of k_G(x-x0) at x0)
    [kGt0min, kGt0max] = minmax_kGt_fun(x0,x(xdx,:));
    
    for vdx=1:size(xv,1)  % interate over all active xv
        
        vectorxvx0 = bsxfun(@minus, xv(vdx,1:3), x0(:,1:3));  % vector xv - x0
        
        % secondary source selection for ALL x0 and a single xv
        active = vectorxvx0*xv(vdx,4:6).' >= 0;
      
        % polar angle and radius of vector (xv-x0)
        phixvx0 = cart2pol(vectorxvx0(active,1),vectorxvx0(active,2));
        
        % k_FS,tv(xv-x0) (tangential component of k_FS(xv-x0) at xv)
        kFStv = sin(phinv(vdx) - phixvx0); 
        % k_FS,t0(xv-x0) (tangential component of k_FS(xv-x0) at x0)
        kFSt0 = sin(phin0(active) - phixvx0); 
        
        % range for aliasing wavelength at x0 caused by discrete SSD (eta = 1)
        lambda10_1 = abs(kFSt0 - kGt0min(active)).*deltax0(active);
        lambda10_2 = abs(kFSt0 - kGt0max(active)).*deltax0(active);

        lambda10min = min(lambda10_1,lambda10_2);
        lambda10max = max(lambda10_1,lambda10_2);
        
        % aliasing wavelength at xv caused by discrete virtual SSD (zeta = 1)
        lambda01 = abs((kStv(vdx) - kFStv).*deltaxv(vdx));
        
        % select the x0 where the virtual SSD wavelength is in the range of SSD 
        % wavelength
        select = lambda01 < (1+eps).*lambda10max & lambda10min < (1+eps).*lambda01;
        if any(select)
          lambda = max(lambda, max(lambda01(select)));
        end
        
        select = lambda01 < eps;
        if any(select)
          lambda = max(lambda, max(lambda10max(select)));
        end
    end
        
    f(xdx) = conf.c./lambda;
end
